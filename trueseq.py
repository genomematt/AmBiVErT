#!/usr/bin/env python
# encoding: utf-8
"""
trueseq.py

Utilities for handling Illumina TruSeq amplicon manifest files

Created by Matthew Wakefield on 2013-05-02.
Copyright (c) 2013  Matthew Wakefield and The University of Melbourne. All rights reserved.
"""
from __future__ import print_function
import sys, os
import argparse
from sequence_utilities import *
from collections import namedtuple

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2013,  Matthew Wakefield and The University of Melbourne"
__credits__ = ["Matthew Wakefield","Graham Taylor"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Matthew Wakefield"
__email__ = "matthew.wakefield@unimelb.edu.au"
__status__ = "Development"


def parse_trueseq_manifest(inputfile):
    def parse_header(infile):
        line = infile.readline().strip('\n').split('\t')
        header = {}
        while line[0] != '[Probes]':
            if line[0] and line[0][0] != '[':
                header[line[0]] = line[1]
            line = infile.readline().strip('\n').split('\t')
            #print('#'*5, repr(line))
        return header
    
    def parse_probes(infile):
        probes = []
        line = infile.readline().strip('\n').split('\t')
        while line[0] != '[Targets]':
            line = infile.readline().strip('\n').split('\t')
        return probes
    
    def parse_targets(infile):
        targets = []
        column_names = ['TargetA', 'TargetB', 'Target Number', 'Chromosome', 'Start Position', 'End Position', 'Probe Strand', 'Sequence', 'Species', 'Build ID']
        TargetRecord = namedtuple('TargetRecord', [ x.replace(' ','_') for x in column_names])
        line = infile.readline().strip('\n').split('\t')
        assert line[:10] == column_names #ignore extra columns that may have been added by excel
        if len(line) > 10:
            print('WARNING Extra columns in mainifest file being ignored:',line[10:],file=sys.stderr)
        line = infile.readline().strip('\n').split('\t')[:10]
        while line[0] != '':
            targets.append(TargetRecord._make(line))
            line = infile.readline().strip('\n').split('\t')[:10]
        return targets
    
    with inputfile as infile:
        header = parse_header(infile)
        probes = parse_probes(infile)
        targets = parse_targets(infile)
    return header, probes, targets

def command_line_interface(*args,**kw):
    parser = argparse.ArgumentParser(description='A script for converting Illumina TrueSeq Amplicon manifest files to fasta target files')
    parser.add_argument('--manifest',
                        type=argparse.FileType('U'),
                        default=sys.stdin,
                        help='an Illumina TrueSeq Amplicon manifest file. Default: stdin')
    parser.add_argument('--output',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='a multi fasta output file of sequence targets. Default: stdout')    
    return parser.parse_args(*args,**kw)

if __name__ == '__main__':
    args = command_line_interface()
    header, probes, targets = parse_trueseq_manifest(args.manifest)
    with args.output as outfile:
        for target in targets:
            name = '{0} {1} {2}-{3}'.format(target.TargetA.split()[0],target.Chromosome,target.Start_Position,target.End_Position)
            outfile.write(format_fasta(name,target.Sequence))
