#!/bin/bash
# Travis install script for AmBiVErT non docker installs
# Basic for now but not done inline as this is likely to get more elaborate

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then

    brew install python3
    brew upgrade python
    pip3 install .

elif [[ ($USE_PYTHON_VERSION == 'py37') && ($TRAVIS_OS_NAME == 'linux') ]]; then

    sudo apt-get install python3.7 python3.7-pip python3.7-dev
    sudo pip3.7 install .
    #export PATH=$PATH:$(dirname $(which python3.7))
    echo $(which python3.7)
    echo $(readlink $(which python3.7))
   

else

    sudo apt-get install python3 python3-pip python3-dev
    sudo pip3 install .

fi

python3 --version
pip3 --version

