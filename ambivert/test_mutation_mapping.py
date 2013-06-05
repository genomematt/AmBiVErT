#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Matthew Wakefield on 2013-05-03.
Copyright (c) 2013  Matthew Wakefield and The Walter and Eliza Hall Institute. All rights reserved.
"""
from __future__ import print_function
import sys
import os
#from sequence_utilities import *

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2013,  Matthew Wakefield and The Walter and Eliza Hall Institute"
__credits__ = ["Matthew Wakefield",]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Development"


def point_mutate_sequence(sequence,chromosome=None,one_based_start=None,one_based_site=1):
    site = one_based_site - 1
    reference = sequence[site]
    for base in ['A','C','G','T']:
        if base == reference:
            continue
        mdz_tag = '{pre}{mut}{post}'.format(pre = site if site else '',
                                            mut = base,
                                            post = (len(sequence)-1)-site if (len(sequence)-1)-site else '')
        name = '{chromosome}_{start}_{cigar_len}M_{mdz_tag}'.format(chromosome=chromosome,
                                                                    start=one_based_start,
                                                                    cigar_len=len(sequence),
                                                                     mdz_tag=mdz_tag)
        yield name,sequence[:site]+base+sequence[site+1:]


if __name__ == '__main__':
    [print(x) for x in point_mutate_sequence('GATC')]
    [print(x) for x in point_mutate_sequence('GATC',chromosome='test',one_based_start='1230789',one_based_site=2)]
    [print(x) for x in point_mutate_sequence('GATC',chromosome='test',one_based_start='1230789',one_based_site=4)]

