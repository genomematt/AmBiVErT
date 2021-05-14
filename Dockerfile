# How to run AmBiVErT in a docker container

# Docker File Author / Maintainer
# MAINTAINER Matthew Wakefield <matthew.wakefield@unimelb.edu.au>

# First you will need install docker.  On MacOS you do
# brew install Caskroom/cask/virtualbox
# brew install docker
# brew install boot2docker
# boot2docker download
# boot2docker init
# boot2docker up
# $(boot2docker shellinit)

# You then need to create a docker container.  In the same directory as this file you type

# docker build -t ambivert .  

# This will run the following executable part of this file

FROM ubuntu:trusty-20190425
MAINTAINER Matthew Wakefield <matthew.wakefield@unimelb.edu.au>
RUN apt-get update && apt-get install -y \
	python3-pip \
	git
#RUN pip install --upgrade pip && \
#	pip install --upgrade setuptools
RUN pip3 install git+https://github.com/genomematt/AmBiVErT.git

# congratulations - you now have a docker container with AmBiVErT installed!

# To test it run
# docker run -i ambivert ambivert --version
# (note the first ambivert is the container name, the second the program name)

# To do anything useful with it you will need to mount some data in the container
# The simplest version of this is to mount a directory from the host.
# This must be an absolute path, and something like $HOME/my/data/directory:/data
# will mount your host directory in the container as /data

# docker run -it -v $HOME/repos/AmBiVErT/ambivert/tests/data/:/data/ ambivert ambivert -f /data/testdata_R1.fastq -r /data/testdata_R2.fastq --fasta/data/testdatareferences.fasta

# If you want to use the parallel processing capabilities of ambivert one option is to define a cpu set
# The docker run option --cpuset=0,1 sets the docker container to run on the first and second cores of the host machine
# The ambivert option --threads 2 tells ambivert to use two concurrent threads for mapping

# docker run -it --cpuset=0,1 -v $HOME/repos/AmBiVErT/ambivert/tests/data/:/data/ ambivert ambivert -f /data/testdata_R1.fastq -r /data/testdata_R2.fastq --fasta /data/testdatareferences.fasta --threads 2

