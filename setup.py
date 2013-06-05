#!/usr/bin/env python

from distutils.core import setup

setup(
    name='AmBiVErT',
    version='0.1.0',
    author='Matthew Wakefield',
    author_email='matthew.wakefield@unimelb.edu.au',
    packages=['ambivert'],
    url='https://git@bitbucket.org/genomematt/ambivert.git',
    license='LICENSE.txt',
    description='AmBiVErT - AMplicon BInning Variant caller with ERror Truncation.\
                 For calling variants in amplicon based sequencing experiments',
    long_description=open('README.txt').read(),
    install_requires=[
        "cogent",
    ],
    classifiers=[
          'Development Status :: 5 - Alpha',
          'License :: OSI Approved :: MIT',
          'Operating System :: POSIX',
          'Programming Language :: Python',
    ],

)
