#!/bin/bash
# Travis install script for AmBiVErT non docker installs
# Basic for now but not done inline as this is likely to get more elaborate

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then

    brew install python3
    pip3 install .

else

    sudo apt-get install python3 python3-pip python3-dev
    sudo pip3 install .

fi

python3 --version
pip3 --version

