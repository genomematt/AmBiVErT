#!/usr/bin/env python3
# encoding: utf-8
"""
truseq_manifest.py

Utilities for handling Illumina TruSeq amplicon manifest files

Created by Matthew Wakefield.
Copyright (c) 2013-2015  Matthew Wakefield and The University of Melbourne. All rights reserved.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

"""
#from __future__ import print_function, division, unicode_literals
import sys, os
import argparse
from ambivert.sequence_utilities import *
from collections import namedtuple

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2013-2015,  Matthew Wakefield and The University of Melbourne"
__credits__ = ["Matthew Wakefield","Graham Taylor"]
__license__ = "GPL"
__version__ = "0.5.2"
__maintainer__ = "Matthew Wakefield"
__email__ = "matthew.wakefield@unimelb.edu.au"
__status__ = "Development"

def parse_truseq_manifest(inputfile):
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
            if len(line) == 19 and line[0] != '[Targets]':
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
    parser = argparse.ArgumentParser(description='A script for converting Illumina TruSeq Amplicon manifest files to fasta files\
                                                    Produces either a fasta file of target sequences without primers or a file\
                                                    of primer sequences suitable for use by a trimming program (eg Nesoni clip)')
    parser.add_argument('--manifest',
                        type=argparse.FileType('rt'),
                        default=sys.stdin,
                        help='an Illumina TruSeq Amplicon manifest file. Default: stdin')
    parser.add_argument('--output',
                        type=argparse.FileType('wt'),
                        default=sys.stdout,
                        help='a multi fasta output file of sequence targets. Default: stdout')
    parser.add_argument('--probes',
                        action="store_true",
                        help='output only the ULSO and DLSO primer sequences')
    parser.add_argument('--adaptors',
                        action="store_true",
                        help='append Illumina adaptor sequences to the primer sequences')
    parser.add_argument('--with_probes',
                        action="store_true",
                        help='append the ULSO and DLSO sequences to the fasta target sequences')
    parser.add_argument('--softmask_probes',
                        action="store_true",
                        help='append the ULSO and DLSO sequences to the fasta target sequences')
    parser.add_argument('--all_plus',
                        action="store_true",
                        help='reorient target sequences so they are all presented on the plus strand')
    return parser.parse_args(*args,**kw)

def make_probes(header, probes, targets, adaptors=False, output=sys.stdout):
    if adaptors:
        #Truseq custom amplicon is P7-index1-adaptor-ULSO-target-DLSO-index2-P5
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
        DLSOadaptorRC = ''
    with output as outfile:
        #print(probes)
        for probe in probes:
            if probe.Target_ID:
                outfile.write(format_fasta(probe.Target_ID.replace(' ','_')+'_ULSO',ULSOadaptor+probe.ULSO_Sequence))
                outfile.write(format_fasta(probe.Target_ID.replace(' ','_')+'_DLSO',probe.DLSO_Sequence+DLSOadaptorRC))
    pass

def make_fasta(header, probes, targets, output=sys.stdout, **kw):
    with output as outfile:
        for name,sequence in make_sequences(header, probes, targets, **kw):
            #print(name,sequence)
            outfile.write(format_fasta(name,sequence))
            
    
def make_sequences(header, probes, targets,  with_probes=False, softmask_probes=False, all_plus=True):
    probes_by_targetid = {probe.Target_ID:(probe.ULSO_Sequence,probe.DLSO_Sequence) for probe in probes}
    for target in targets:
        if with_probes:
            name = '{0}_{1}_{2}_{3} {1} {2} {3} {4}'.format(target.TargetB.replace(' ','_'), target.Chromosome, target.Start_Position, target.End_Position, target.Probe_Strand)
        else:
            if all_plus and (target.Probe_Strand == '-'):
                start = int(target.Start_Position) + len(probes_by_targetid[target.TargetA][1]) # add ULSO length to start
                end = int(target.End_Position) - len(probes_by_targetid[target.TargetA][0]) # subtract DLSO length from end
            else:
                start = int(target.Start_Position) + len(probes_by_targetid[target.TargetA][0]) # add DLSO length to start
                end = int(target.End_Position) - len(probes_by_targetid[target.TargetA][1]) # subtract ULSO length from end
            name = '{0}_{1}_{2}_{3} {1} {2} {3} {4}'.format(target.TargetB.replace(' ','_'), target.Chromosome, start, end, target.Probe_Strand)
        
        if with_probes and softmask_probes:
            # For on target sequences Submitted Target Region Strand == probe Probe Strand != target Probe Strand
            # So for the order should be ULSO - Reverse complement of target - DLSO
            # RC(DLSO)-target-RC(ULSO) is equivalent and puts the sequence in read_1 format
            seq = reverse_complement(probes_by_targetid[target.TargetA][1].lower())+target.Sequence+reverse_complement(probes_by_targetid[target.TargetA][0].lower())
        elif with_probes:
            seq = reverse_complement(probes_by_targetid[target.TargetA][1])+target.Sequence+reverse_complement(probes_by_targetid[target.TargetA][0])
        else:
            seq = target.Sequence
        if all_plus and (target.Probe_Strand == '-'):
            if name[-1] == '-':
                name = name[:-1]+'+'
            #print('minus strand - rc ing',name,target.Start_Position, target.End_Position, len(probes_by_targetid[target.TargetA][0]), len(probes_by_targetid[target.TargetA][1]))
            yield(name,reverse_complement(seq))
        else:
            #print('plus strand',name,target.Start_Position, target.End_Position, len(probes_by_targetid[target.TargetA][0]), len(probes_by_targetid[target.TargetA][1]))
            yield(name,seq)

def main():
    args = command_line_interface()
    # Probably should refactor these so functions just yeild (name,sequence) and format output here.
    if args.probes:
        make_probes(*parse_truseq_manifest(args.manifest), adaptors=args.adaptors, output=args.output)
    else:
        make_fasta(*parse_truseq_manifest(args.manifest), output=args.output, with_probes=args.with_probes, softmask_probes=args.softmask_probes, all_plus=args.all_plus)
    pass

if __name__ == '__main__':
    main()
