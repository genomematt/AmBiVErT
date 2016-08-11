#!/bin/bash
# Travis install script for AmBiVErT non docker installs
# Basic for now but not done inline as this is likely to get more elaborate

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then

    brew install python3

else

    sudo apt-get install python3.5

fi

python3 --version
pip3 --version

pip3 install .
