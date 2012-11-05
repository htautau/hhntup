# This script will work in either bash or zsh.

# deterine path to this script
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in
SOURCE_HIGGSTAUTAU_SETUP="${BASH_SOURCE[0]:-$0}"

DIR_HIGGSTAUTAU_SETUP="$( dirname "$SOURCE_HIGGSTAUTAU_SETUP" )"
while [ -h "$SOURCE_HIGGSTAUTAU_SETUP" ]
do 
  SOURCE_HIGGSTAUTAU_SETUP="$(readlink "$SOURCE_HIGGSTAUTAU_SETUP")"
  [[ $SOURCE_HIGGSTAUTAU_SETUP != /* ]] && SOURCE_HIGGSTAUTAU_SETUP="$DIR_HIGGSTAUTAU_SETUP/$SOURCE_HIGGSTAUTAU_SETUP"
  DIR_HIGGSTAUTAU_SETUP="$( cd -P "$( dirname "$SOURCE_HIGGSTAUTAU_SETUP"  )" && pwd )"
done
DIR_HIGGSTAUTAU_SETUP="$( cd -P "$( dirname "$SOURCE_HIGGSTAUTAU_SETUP" )" && pwd )"

echo "sourcing ${SOURCE_HIGGSTAUTAU_SETUP}..."

export PATH=${DIR_HIGGSTAUTAU_SETUP}${PATH:+:$PATH}
export PYTHONPATH=${DIR_HIGGSTAUTAU_SETUP}${PYTHONPATH:+:$PYTHONPATH}

if [ -f ${DIR_HIGGSTAUTAU_SETUP}/externaltools/setup.sh ]
then
    source ${DIR_HIGGSTAUTAU_SETUP}/externaltools/setup.sh
fi
