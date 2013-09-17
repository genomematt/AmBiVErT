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
        column_names = ['Target Region Name', 'Target Region ID', 'Target ID', 'Species', 'Build ID', 'Chromosome', 'Start Position', 'End Position',
                        'Submitted Target Region Strand', 'ULSO Sequence', 'ULSO Genomic Hits', 'DLSO Sequence', 'DLSO Genomic Hits', 'Probe Strand',
                        'Designer', 'Design Score', 'Expected Amplifed Region Size', 'SNP Masking', 'Labels']
        ProbeRecord = namedtuple('TargetRecord', [ x.replace(' ','_') for x in column_names])
        line = infile.readline().strip('\n').split('\t')
        assert line == column_names
        while line[0] != '[Targets]':
            line = infile.readline().strip('\n').split('\t')
            probes.append(ProbeRecord._make(line))
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
    parser = argparse.ArgumentParser(description='A script for converting Illumina TrueSeq Amplicon manifest files to fasta files\
                                                    Produces either a fasta file of target sequences without primers or a file\
                                                    of primer sequences suitable for use by a trimming program (eg Nesoni clip)')
    parser.add_argument('--manifest',
                        type=argparse.FileType('U'),
                        default=sys.stdin,
                        help='an Illumina TrueSeq Amplicon manifest file. Default: stdin')
    parser.add_argument('--output',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='a multi fasta output file of sequence targets. Default: stdout')
    parser.add_argument('--probes',
                        action="store_true",
                        help='output only the primer sequences')
    parser.add_argument('--adaptors',
                        action="store_true",
                        help='append Illumina adaptor sequences to the primer sequences')
    return parser.parse_args(*args,**kw)

def main():
    args = command_line_interface()
    header, probes, targets = parse_trueseq_manifest(args.manifest)
    if args.probes:
        print('in probes')
        if args.adaptors:
            #Trueseq custom amplicon is P7-index1-adaptor-ULSO-target-DLSO-index2-P5
            #Oligonucleotide sequences copyright 2007-2012 Illumina Inc.  All rights reserved
            ULSOadaptor = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT' # P7 end
            #DLSOadaptor = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'
            DLSOadaptorRC = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'# Reverse complemented P5 end
            #The read 1 sequencing primer reads into the reverse complement of the DLSO
            #To generate a clipping sequence we RC the seq primer and append at end of DLSO
            #Stop at the index sequence as it is unlikely to be included and is variable
            #end copyrighted sequences
        else:
            ULSOadaptor = ''
            DLSOadaptor = ''
        with args.output as outfile:
            print(probes)
            for probe in probes:
                if probe.Target_ID:
                    outfile.write(format_fasta(probe.Target_ID+' ULSO',ULSOadaptor+probe.ULSO_Sequence))
                    outfile.write(format_fasta(probe.Target_ID+' DLSO',probe.DLSO_Sequence+DLSOadaptorRC))
    else:
        with args.output as outfile:
            for target in targets:
                name = '{0} {1} {2}-{3}'.format(target.TargetA.split()[0],target.Chromosome,target.Start_Position,target.End_Position)
                outfile.write(format_fasta(name,target.Sequence))
    
    
if __name__ == '__main__':
    main()
