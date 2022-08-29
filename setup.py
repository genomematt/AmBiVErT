#!/usr/bin/env python
import sys, types
#be assertive that we want a shared library on linux
try:
  import dl
except ImportError:
  sys.modules['dl'] = types.ModuleType('dl')
  sys.modules['dl'].RTLD_NOW = None

from setuptools import setup

try:
    from setuptools.extension import Library
except ImportError:
    from setuptools import Library

import glob

try:
    from platform import linux_distribution
    from os import environ 
    if linux_distribution()[0] == 'CentOS':
        environ['CFLAGS'] = "-std=gnu99"
        print('setting environ["CFLAGS"] = "-std=gnu99"')
except:
    pass


setup(
    name='ambivert',
    version='0.5.2',
    author='Matthew Wakefield',
    author_email='matthew.wakefield@unimelb.edu.au',
    install_requires = [
      'setuptools >= 1.1.6',
    ],
    packages=['ambivert',
              'ambivert.align',
              'ambivert.tests',
              'ambivert.tests.data',
              'lib.align',
              ],
    zip_safe = False,
    ext_modules = [
      Library(
        'ambivert.align.align_c',
        sources = glob.glob('lib/align/*.c'),
        depends = glob.glob('lib/align/*.h')
      ),
    ],
    include_package_data = True,
    url='https://git@bitbucket.org/genomematt/ambivert.git',
    license='GPLv3',
    entry_points={
        'console_scripts': ['ambivert = ambivert.ambivert:main',
                            'ambivert_truseq_manifest = ambivert.truseq_manifest:main',
                            'ambivert_simulate_variants = ambivert.simulate_variants:main',
                            'ambivert_insert_variants = ambivert.insert_variants:main',
                           ]
    },
    test_suite = "ambivert.tests.test_ambivert",

    description='AmBiVErT - AMplicon BInning Variant caller with ERror Truncation.\
                 For calling variants in amplicon based sequencing experiments',
    long_description='AmBiVErT - AMplicon BInning Variant caller with ERror Truncation.\
                 For calling variants in amplicon based sequencing experiments',
    classifiers=[
          'Development Status :: 5 - Production/Stable',
          'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
          'Operating System :: POSIX',
          'Programming Language :: Python :: 3.3',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Programming Language :: Python :: 3.9',
          'Programming Language :: Python :: 3.10',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],

)
