#!/usr/bin/env python3
# encoding: utf-8
"""
insilicoPCR_to_amplicons.py

A script for formating in silico PCR fasta files from UCSC as fasta suitable as use for a reference
in ambivert analysis.

https://genome.ucsc.edu/cgi-bin/hgPcr

Created by Matthew Wakefield.
Copyright (c) 2013-2016  Matthew Wakefield and The University of Melbourne. All rights reserved.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

"""
#from __future__ import print_function, division, unicode_literals
import sys, os
import argparse
from ambivert.sequence_utilities import *

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2013-2016,  Matthew Wakefield and The University of Melbourne"
__credits__ = ["Matthew Wakefield"]
__license__ = "GPL"
__version__ = "0.5.1"
__maintainer__ = "Matthew Wakefield"
__email__ = "matthew.wakefield@unimelb.edu.au"
__status__ = "Development"


def command_line_interface(*args,**kw): #pragma no cover
    parser = argparse.ArgumentParser(description='A script for formating in silico PCR fasta files from UCSC \
                                              as fasta suitable as use for a reference in ambivert analysis.')
    parser.add_argument('--fasta',
                        type=argparse.FileType('rt'),
                        default=sys.stdin,
                        help='UCSC in silico PCR fasta. Default: stdin')
    parser.add_argument('--output',
                        type=argparse.FileType('wt'),
                        default=sys.stdout,
                        help='a multi fasta output file of sequence targets. Default: stdout')
    parser.add_argument('--softmask_probes',
                        action="store_true",
                        help='softmask primer sequences in lower case')
    return parser.parse_args(*args,**kw)
    
def process_isp_fasta(name, seq, target_name='UNKNOWN',softmask_probes=True):
    location, size, primer1, primer2 = name.split(' ')
    #location is of format 'chr22:31000551+31001000'
    chrom = location.split(':')[0]
    start = location.split(':')[1].split('+')[0].split('-')[0]
    end = location.split(':')[1].split('+')[-1].split('-')[-1]
    strand = '+' if location.split(':')[1].count('+') else '-'
    without_primer = seq[len(primer1):-len(primer2)].upper()
    if softmask_probes:
        primer1 = primer1.lower()
        primer2 = primer2.lower()
    return target_name, chrom, start, end, strand, primer1+without_primer+primer2
        

def make_fasta(infile=sys.stdin, output=sys.stdout, **kw):
    with output as outfile, infile as isp_fasta:
        for name,seq in parse_fasta(isp_fasta):
            target_name, chrom, start, end, strand, sequence = process_isp_fasta(name, seq, **kw)
            description = '{0}_{1}_{2}_{3}_{4} {1} {2} {3} {4}'.format(target_name, chrom, start, end, strand)
            #print(description,sequence)
            outfile.write(format_fasta(description,sequence))

def main(): #pragma no cover
    args = command_line_interface()
    make_fasta(infile=args.infile, output=args.outfile, softmask_probes=args.softmask_probes)
    pass

if __name__ == '__main__': #pragma no cover
    main()
    

