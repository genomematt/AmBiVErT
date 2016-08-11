#!/bin/bash
# Travis install script for AmBiVErT non docker installs
# Basic for now but not done inline as this is likely to get more elaborate

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then

    brew install python3

fi

python --version
pip3 --version

pip3 install .
