#!/usr/bin/env python

# from distutils.core import setup
from setuptools import setup

setup(
    name='AmBiVErT',
    version='0.2.dev1',
    author='Matthew Wakefield',
    author_email='matthew.wakefield@unimelb.edu.au',
    packages=['ambivert'],
    include_package_data = True,
    url='https://git@bitbucket.org/genomematt/ambivert.git',
    license='GPLv3',
    entry_points={
        'console_scripts': ['ambivert = ambivert.ambivert:main',
                            'truseq_manifest = ambivert.truseq_manifest:main',
                            'simulate_variants = ambivert.simulate_variants:main',
                           ]
    },

    description='AmBiVErT - AMplicon BInning Variant caller with ERror Truncation.\
                 For calling variants in amplicon based sequencing experiments',
    long_description=open('README.txt').read(),
    classifiers=[
          'Development Status :: 3 - Alpha',
          'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
          'Operating System :: POSIX',
          'Programming Language :: Python :: 3.3',
          'Programming Language :: Python :: 3.4',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Private :: Not yet ready for uploading to PyPI',
    ],

)
