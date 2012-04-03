#!/bin/bash

BASE=${PWD}

if [[ "${1}" = "clean" ]]
then
    echo "Cleaning up..."
    rm -rf rootpy
    rm -rf goodruns
    rm -rf cython
    rm -rf lxml
    rm -rf yaml
    rm -rf atlastools
    rm -rf python
    rm -rf root
    rm -rf rpmroot
    rm -rf tabulartext
    exit
fi

repo='http://linuxsoft.cern.ch/cern/slc5X/updates/x86_64/RPMS/'
GIT_USER=ndawe

# install python and ROOT

PYTHON_VERS_MAJOR=2.7
PYTHON_VERS=2.7.2

use_precompiled_python=true
use_precompiled_root=true

if [[ ! -e python ]]
then
    if $use_precompiled_python
    then
        wget http://hep.phys.sfu.ca/~endw/grid/python.tar.gz
        tar -zxf python.tar.gz
        rm -rf python.tar.gz
    else
        wget http://www.python.org/ftp/python/${PYTHON_VERS}/Python-${PYTHON_VERS}.tar.bz2
        tar -xjf Python-${PYTHON_VERS}.tar.bz2
        cd Python-${PYTHON_VERS}
        ./configure --enable-shared --prefix=${BASE}/python
        make -j
        make install
        cd ..
        rm -rf Python-${PYTHON_VERS}
        rm -rf Python-${PYTHON_VERS}.tar.bz2
    fi
fi

export PATH=${BASE}/python/bin${PATH:+:$PATH}
export PYTHONPATH=${BASE}/python/lib/python${PYTHON_VERS_MAJOR}/site-packages${PYTHONPATH:+:$PYTHONPATH}
export LD_LIBRARY_PATH=${BASE}/python/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}

ROOT_VERS=5.30.00

if [[ ! -e root ]]
then
    if ${use_precompiled_root}
    then
        wget http://hep.phys.sfu.ca/~endw/grid/root.tar.gz
        tar -zxf root.tar.gz
        rm -rf root.tar.gz
    else
        wget ftp://root.cern.ch/root/root_v${ROOT_VERS}.source.tar.gz
        gzip -dc root_v${ROOT_VERS}.source.tar.gz | tar -xf -
        rm -rf root_v${ROOT_VERS}.source.tar.gz
        cd root
        ./configure --with-python-libdir=${BASE}/python/lib --with-python-incdir=${BASE}/python/include/python${PYTHON_VERS_MAJOR} \
                    --with-dcap-libdir=/atlas/software/ATLASLocalRootBase/x86_64/gLite/current/d-cache/dcap/lib64 \
                    --with-dcap-incdir=/atlas/software/ATLASLocalRootBase/x86_64/gLite/current/d-cache/dcap/include
        make -j
        cd ..
    fi
fi

source root/bin/thisroot.sh

python_version=`python -c "import distutils.sysconfig; print distutils.sysconfig.get_python_version()"`
PYTHON_LIB=`python -c "import distutils.sysconfig; import os; print os.path.dirname(distutils.sysconfig.get_python_lib(standard_lib=True))"`
echo "Python version "${python_version}
echo "Python lib located in "${PYTHON_LIB}

function download_python_module_from_github() {
    PACKAGE=${1}
    ROOT=https://raw.github.com/${GIT_USER}/${PACKAGE}/master/
    if [[ ! -e ${PACKAGE} ]]
    then
        mkdir ${PACKAGE}
        cd ${PACKAGE}
        wget --no-check-certificate ${ROOT}MANIFEST
        while read line;
        do
            if [[ ${line} != \#* ]];
            then
                echo "${ROOT}${line}" >> MANIFEST_1
            fi
        done < MANIFEST
        cat MANIFEST_1
        wget -x -nH --cut-dirs=3 --no-check-certificate -i MANIFEST_1
        cd ${BASE}
    fi
}

function install_python_module() {
    cd ${BASE}
    if [[ ! -e ${1} ]]
    then
        echo "Checking out ${1}..."
        if svn checkout ${2} ${1}
        then
            : # do nothing
        else
            echo "Subversion checkout failed."
            echo "wget'ting instead..."
            download_python_module_from_github ${1}
        fi
        if [[ -d ${1} ]]
        then
            echo "Installing ${1}..."
            cd ${1}
            echo ">>> lib dirs: ${PYTHON_LIB}:${BASE}/rpmroot/usr/lib64"
            echo ">>> include dirs: ${RPM_INCLUDE}"
            python setup.py build_ext --library-dirs="${PYTHON_LIB}:${BASE}/rpmroot/usr/lib64" --include-dirs="${RPM_INCLUDE}"
            python setup.py build -e "/usr/bin/env python"
            python setup.py install
            cd ..
        fi
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
    wget ${1}${2}
    rpm2cpio ${2} | cpio -idmv
    rm -rf ${2}
    cd ${BASE}
    export RPM_INCLUDE=${BASE}/rpmroot/usr/include/${3}${RPM_INCLUDE:+:$RPM_INCLUDE}
}

export RPM_INCLUDE=${BASE}/rpmroot/usr/include${RPM_INCLUDE:+:$RPM_INCLUDE}

if [[ ! -e rpmroot ]]
then
    install_rpm ${repo} libxml2-2.6.26-2.1.12.x86_64.rpm libxml2
    install_rpm ${repo} libxml2-devel-2.6.26-2.1.12.x86_64.rpm libxml2
    install_rpm ${repo} libxslt-1.1.17-2.el5_2.2.x86_64.rpm libxslt
    install_rpm ${repo} libxslt-devel-1.1.17-2.el5_2.2.x86_64.rpm libxslt
fi

export LD_LIBRARY_PATH=${BASE}/rpmroot/usr/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export PATH=${BASE}/rpmroot/usr/bin${PATH:+:$PATH}

if ! ${use_precompiled_python}
then
    install_python_module cython http://svn.github.com/cython/cython.git
    install_python_module lxml http://svn.github.com/lxml/lxml.git
    install_python_module yaml http://svn.pyyaml.org/pyyaml/tags/3.10/
fi

install_python_module rootpy http://svn.github.com/ndawe/rootpy.git
install_python_module goodruns http://svn.github.com/ndawe/goodruns.git
install_python_module atlastools http://svn.github.com/ndawe/atlastools.git
install_python_module tabulartext http://svn.github.com/ndawe/tabulartext.git
