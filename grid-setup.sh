#!/bin/bash
# Author: Noel Dawe
shopt -s expand_aliases

function print_help() {
    echo "Usage : $0 [clean|local|unpack|build-root-python|build|worker]"
    exit
}

if [[ $# -eq 0 ]]
then
    print_help
fi

BASE=${PWD}
github_deps=deps/packages.github
repo='http://linuxsoft.cern.ch/cern/slc5X/updates/x86_64/RPMS/'

PYTHON_VERS_MAJOR=2.7
PYTHON_VERS=2.7.2
ROOT_VERS=5.34.14
SETUPTOOLS_VERS=2.2

ROOT_VERSION_CVMFS=5.34.14-x86_64-slc5-gcc4.3
PYTHON_VERSION_CVMFS=2.6.5-x86_64-slc5-gcc43

USE_CVMFS=true

# debugging grid issues
echo "Python site imported from:"
python -c "import site; print site.__file__"
# clear PYTHONPATH
unset PYTHONPATH


function download_from_github() {
    cd ${BASE}
    GIT_USER=${1}
    PACKAGE=${2}
    TAG=${3}
    if [[ ! -e ${PACKAGE}.tar.gz ]]
    then
        if ! wget --no-check-certificate -O ${PACKAGE}.tar.gz https://github.com/${GIT_USER}/${PACKAGE}/tarball/${TAG}
        then
            echo "Failed to download package ${PACKAGE} from github"
            exit 1
        fi
    fi
}

function unpack_github_tarball() {
    cd ${BASE}
    GIT_USER=${1}
    PACKAGE=${2}
    if [[ ! -e ${PACKAGE} ]]
    then
        if tar -pzxf ${PACKAGE}.tar.gz
        then
            rm -f ${PACKAGE}.tar.gz
            mv ${GIT_USER}-${PACKAGE}-* ${PACKAGE}
        else
            exit 1
        fi
    fi
}

function install_python_package() {
    cd ${BASE}
    PACKAGE=${1}
    if [[ -d ${PACKAGE} ]]
    then
        echo "Installing ${PACKAGE}..."
        cd ${PACKAGE}
        cp ${BASE}/setuptools-${SETUPTOOLS_VERS}.tar.gz .
        #echo ">>> lib dirs: ${PYTHON_LIB}:${BASE}/rpmroot/usr/lib64"
        #echo ">>> include dirs: ${RPM_INCLUDE}"
        #python setup.py build_ext --library-dirs="${PYTHON_LIB}:${BASE}/rpmroot/usr/lib64" --include-dirs="${RPM_INCLUDE}"
        #python setup.py build -e "/usr/bin/env python"
        #python setup.py setopt --command build_ext --option library-dirs --set-value ${PYTHON_LIB}
        if ! python setup.py install --user
        then
            echo "Failed to install package ${PACKAGE}"
            exit 1
        fi
        cd ..
    fi
}

function install_rpm() {
    cd ${BASE}
    if [[ ! -e rpmroot ]]
    then
        mkdir rpmroot
    fi
    echo "Installing ${2}..."
    cd rpmroot
    if ! wget ${1}${2}
    then
        echo "Failed to download ${1}${2}"
        exit 1
    fi
    rpm2cpio ${2} | cpio -idmv
    rm -rf ${2}
    cd ${BASE}
    export RPM_INCLUDE=${BASE}/rpmroot/usr/include/${3}${RPM_INCLUDE:+:$RPM_INCLUDE}
}

function setup_rpm() {
    export RPM_INCLUDE=${BASE}/rpmroot/usr/include${RPM_INCLUDE:+:$RPM_INCLUDE}
    export LD_LIBRARY_PATH=${BASE}/rpmroot/usr/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
    export PATH=${BASE}/rpmroot/usr/bin${PATH:+:$PATH}
}

function setup_root() {
    source ${BASE}/root/bin/thisroot.sh
}

function setup_python() {
    export PYTHONPATH=${BASE}/python/lib/python${PYTHON_VERS_MAJOR}/site-packages${PYTHONPATH:+:$PYTHONPATH}
    export LD_LIBRARY_PATH=${BASE}/python/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
}

function determine_python() {
    export PYTHON_VERSION=`python -c "import distutils.sysconfig; print distutils.sysconfig.get_python_version()"`
    export PYTHON_LIB=`python -c "import distutils.sysconfig; import os; print os.path.dirname(distutils.sysconfig.get_python_lib(standard_lib=True))"`
    export PYTHON_INCLUDE=`python -c "import distutils.sysconfig; import os; print os.path.dirname(distutils.sysconfig.get_python_inc())"`/python${PYTHON_VERSION}
    echo "Python version is "${PYTHON_VERSION}
    echo "Python lib is located in "${PYTHON_LIB}
    echo "Python include path is "${PYTHON_INCLUDE}
    export PYTHONUSERBASE=${BASE}/user-python
    export PATH=${PYTHONUSERBASE}/bin${PATH:+:$PATH}
    export LD_LIBRARY_PATH=$PYTHON_LIB${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
    export LDFLAGS="-L${PYTHON_LIB} ${LDFLAGS}"
    export CPPFLAGS="-I${PYTHON_INCLUDE} ${CPPFLAGS}"
}

function setup_CVMFS() {
    export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
    source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh
}

function setup_ROOT_CVMFS() {
    source ${ATLAS_LOCAL_ROOT_BASE}/packageSetups/atlasLocalROOTSetup.sh \
        --skipConfirm --rootVersion=${ROOT_VERSION_CVMFS}
}

function setup_python_CVMFS() {
    source ${ATLAS_LOCAL_ROOT_BASE}/packageSetups/atlasLocalPythonSetup.sh \
        --pythonVersion=${PYTHON_VERSION_CVMFS}
}

case "${1}" in
clean)

    echo "Cleaning up..."
    rm -rf user-python 
    if [[ -f ${github_deps} ]]
    then
        while read line
        do
            line=($line)
            package=${line[1]}
            rm -rf ${package}
            rm -f ${package}.tar.gz
        done < ${github_deps}
    fi
    if [[ -f deps/dependencies ]]
    then
        cd deps
        while read line
        do
            tokens=($line)
            rm -rf ${tokens[0]}
            rm -f ${tokens[0]}.tar.gz
        done < dependencies
        cd ${BASE}
    fi
    rm -f ${BASE}/setuptools-${SETUPTOOLS_VERS}.tar.gz
    ;;

local)
    
    if [[ -f deps/dependencies ]]
    then
        cd deps
        while read line
        do
            tokens=($line)
            package=${tokens[0]}
            url=${tokens[1]}
            if [[ ! -e ${package}.tar.gz ]]
            then
                wget --no-check-certificate ${url}
            fi
        done < dependencies
        cd ${BASE}
    fi
    if [[ -f ${github_deps} ]]
    then
        while read line
        do
            download_from_github $line
        done < ${github_deps}
    fi
    # get setuptools
    if [[ ! -e setuptools-${SETUPTOOLS_VERS}.tar.gz ]]
    then
        wget --no-check-certificate https://pypi.python.org/packages/source/s/setuptools/setuptools-${SETUPTOOLS_VERS}.tar.gz
    fi
    ;;

build-root-python)

    # install python and ROOT

    if [[ ! -e python ]]
    then
        if ! wget http://www.python.org/ftp/python/${PYTHON_VERS}/Python-${PYTHON_VERS}.tar.bz2
        then
            echo "Failed to download Python"
            exit 1
        fi
        tar -xjf Python-${PYTHON_VERS}.tar.bz2
        cd Python-${PYTHON_VERS}
        ./configure --enable-shared --prefix=${BASE}/python
        make -j
        make install
        cd ..
        rm -rf Python-${PYTHON_VERS}
        rm -rf Python-${PYTHON_VERS}.tar.bz2
    fi

    setup_python

    if [[ ! -e root ]]
    then
        if ! wget ftp://root.cern.ch/root/root_v${ROOT_VERS}.source.tar.gz
        then
            echo "Failed to download ROOT"
            exit 1
        fi
        gzip -dc root_v${ROOT_VERS}.source.tar.gz | tar -xf -
        rm -rf root_v${ROOT_VERS}.source.tar.gz
        cd root
        ./configure --with-python-libdir=${BASE}/python/lib --with-python-incdir=${BASE}/python/include/python${PYTHON_VERS_MAJOR} \
                    --with-dcap-libdir=/atlas/software/ATLASLocalRootBase/x86_64/gLite/current/d-cache/dcap/lib64 \
                    --with-dcap-incdir=/atlas/software/ATLASLocalRootBase/x86_64/gLite/current/d-cache/dcap/include
        make -j
        cd ..
    fi
    
    setup_root
    ;;

unpack)
    
    if [[ -f ${github_deps} ]]
    then
        while read line
        do 
            unpack_github_tarball $line
        done < ${github_deps}
    fi
    ;;

build)
    
    if $USE_CVMFS
    then
        setup_CVMFS
        setup_python_CVMFS
        setup_ROOT_CVMFS
    fi
    
    determine_python
    export ROOTPY_NO_EXT=1
 
    if [[ ! -e user-python ]]
    then
        echo "Creating user python area"
        mkdir -p user-python/lib/python${PYTHON_VERSION}/site-packages/
        mkdir user-python/bin
    fi

    # install setuptools
    tar -zxvf ${BASE}/setuptools-${SETUPTOOLS_VERS}.tar.gz
    cd setuptools-${SETUPTOOLS_VERS}
    python setup.py install --user
    cd -

    if [[ -f deps/dependencies ]]
    then
        cd deps
        while read line
        do
            tokens=($line)
            package=${tokens[0]}
            tar -zxvf ${package}.tar.gz
            cd ${package}
            python setup.py install --user
            cd ..
        done < dependencies
        cd ${BASE}
    fi
    if [[ -f ${github_deps} ]]
    then
        while read line
        do 
            unpack_github_tarball $line
            line=($line)
            install_python_package ${line[1]}
        done < ${github_deps}
    fi
    # source user build script
    if [[ -f grid.build ]]
    then
        if ! source grid.build
        then
            echo "Failed to execute user build script"
            exit 1
        fi
    fi
    ;;

worker)
    
    if $USE_CVMFS
    then
        setup_CVMFS
        setup_python_CVMFS
        setup_ROOT_CVMFS
    fi

    if [[ -e python ]]
    then
        setup_python
    fi
    if [[ -e root ]]
    then
        setup_root
    fi
    
    determine_python 
    export ROOTPY_GRIDMODE=true

    # source user setup script
    if [[ -f grid.setup ]]
    then
        source grid.setup
    fi
    ;;

*)
    
    print_help
    ;;
esac
