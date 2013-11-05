#!/usr/bin/env python

from distutils.core import setup

setup(
    name='AmBiVErT',
    version='0.1.0',
    author='Matthew Wakefield',
    author_email='matthew.wakefield@unimelb.edu.au',
    packages=['ambivert'],
    include_package_data = True,
    url='https://git@bitbucket.org/genomematt/ambivert.git',
    license='GPL',
    entry_points={
        'console_scripts': ['ambivert = ambivert.ambivert:main',
                            'truseq_manifest = ambivert.truseq_manifest:main',
                            'simulate_mutations = ambivert.simulate_mutations:main',
                           ]
    },

    description='AmBiVErT - AMplicon BInning Variant caller with ERror Truncation.\
                 For calling variants in amplicon based sequencing experiments',
    long_description=open('README.txt').read(),
    classifiers=[
          'Development Status :: 3 - Alpha',
          'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
          'Operating System :: POSIX',
          'Programming Language :: Python :: 2.7',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],

)
