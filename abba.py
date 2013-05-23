#!/usr/bin/env python
# encoding: utf-8
"""
abba.py

Amplicon Batching Before Alignment - ABBA

This program is designed for the processing of amplicon based resequencing data
generated on second generation sequencing platforms (eg Illumina HiSeq/MiSeq).

The approach used is to batch reads derived from each amplicon together in a
clustering step prior to aligning the sequence to a reference genome.

Steps:
    1)  Cluster similar reads
    2)  Align for overlap
    3a) Align to reference target sequences
    3b) Use external program to align to genome
    4)  Consolidate calls
    5)  Output VCF format calls


Created by Matthew Wakefield and Graham Taylor.
Copyright (c) 2013  Matthew Wakefield and The University of Melbourne. All rights reserved.
"""
from __future__ import print_function
import sys
import os
import itertools
from collections import defaultdict
from cogent.align.algorithm import nw_align, sw_align
from sequence_utilities import parse_fastq, parse_fasta, reverse_complement, flatten_paired_alignment, format_alignment

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2013,  Matthew Wakefield and The University of Melbourne"
__credits__ = ["Matthew Wakefield","Graham Taylor"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Matthew Wakefield"
__email__ = "matthew.wakefield@unimelb.edu.au"
__status__ = "Development"


def smith_waterman(seq1,seq2):
    align_seq1, align_seq2 = sw_align(seq1,seq2)
    start_seq1 = seq1.find(align_seq1.replace('-',''))
    start_seq2 = seq2.find(align_seq2.replace('-',''))
    return align_seq1,align_seq2,start_seq1,start_seq2
    


class AmpliconData(object):
    """A Class for holding read data from amplicon experiments
    Reads are indexed by hashes of both the forward and reverse reads
    Arguments:
        - trim5 : number of bases to trim from 5' default = None
        - trim3   : number of bases to trim from 3' default = None
    """
    def __init__(self, trim5=None, trim3=None):
        self.data = defaultdict(list)
        self.merged = {} # a dictionary of overlapped merged sequences
        self.aligned = {} # a dictionary of aligned sequence tuples with same key as self.data
        self.location = {} # dictionary of chromosome
        self.trim5 = trim5
        self.trim3 = None if trim3==None else -1*abs(trim3)
        pass
    
    def __str__(self):
        return self.data.__str__()
    
    def add_reads(self, f_name, f_seq, f_qual,r_name, r_seq, r_qual):
        amplicon_key = f_seq[self.trim5:self.trim3]+r_seq[self.trim5:self.trim3]
        self.data[amplicon_key].append(((f_name, f_seq, f_qual),(r_name, r_seq, r_qual)))
        pass
    
    def process_twofile_readpairs(self, forward_file, reverse_file, parser=parse_fastq):
        for (f_name, f_seq, f_qual),(r_name, r_seq, r_qual) in itertools.izip(parser(forward_file), parser(reverse_file)):
            assert f_name.split()[0] == r_name.split()[0]
            self.add_reads(f_name, f_seq, f_qual,r_name, r_seq, r_qual)
        pass
    
    #def process_interleaved_readpairs(parser=parse_fastq, filename):
    #    pass
    
    def get_above_threshold(self, threshold=1):
        return [x for x in self.data if len(self.data[x]) > threshold]
    
    
        


# Detect overlaps
# We expect overlaps like these when amplicon < read length 
#          ---------->
#          |||||
#    <----------
#
# We expect mainly overlaps like these
#     ---------->
#          |||||
#         <----------
#
# Non overlaps may still overlap but by less than min_overlap
# Need to take care when aligning to genome.
# Do each end then reconstruct
#    ---------->
#         
#               <----------


def match(a, b, min_match = 10, max_mismatch = 1):
    # Written by Toby Sargeant
    """Overlap match detection"""
    m = numpy.zeros((len(a),len(b)))
    for i in xrange(len(a)):
        m[i,0] = +1 if a[i] == b[0] else -1
    for j in xrange(len(b)):
        m[0,j] = +1 if a[0] == b[j] else -1
    for i in xrange(1, len(a)):
        for j in xrange(1, len(b)):
            m[i,j] = m[i-1,j-1] + (+1 if a[i] == b[j] else -1)
    print(repr(m))
    best = (-1, None)
    for j in xrange(min_match, len(b)):
        if m[-1,j] >= j + 1 - max_mismatch:
            if m[-1,j] > best[0]:
                best = (m[-1,j], len(a) - 1 - j)
    for i in xrange(min_match, len(a)):
        if m[i,-1] >= i + 1 - max_mismatch:
            if m[i,-1] > best[0]:
                best = (m[i,-1], -len(b) + 1 + i)
    return best


def process_commandline_args():
    pass

def main(forward_file, reverse_file, threshold=1):
    amplicons = AmpliconData()
    amplicons.process_twofile_readpairs(forward_file, reverse_file)
    #print(amplicons)
    print('Reads that occur more than {0} times'.format(threshold))
    for key in amplicons.get_above_threshold(threshold):
        [print(x[0][1]) for x in amplicons.data[key]]
        print()
        
    for key in amplicons.get_above_threshold():
        print('Smith Waterman Aligned Pairs:')
        fwd = amplicons.data[key][0][0][1]
        rev = reverse_complement(amplicons.data[key][0][1][1])
        fwd_aligned, rev_aligned, fwd_start, rev_start = smith_waterman(fwd,rev)
        print('Forward ',' '*max((rev_start-fwd_start),0)+fwd[:fwd_start]+fwd_aligned+fwd[fwd_start+len(fwd_aligned.replace('-','')):])
        print('Reverse ',' '*max((fwd_start-rev_start),0)+rev[:rev_start]+rev_aligned+rev[rev_start+len(rev_aligned.replace('-','')):])
        overlapped = flatten_paired_alignment(fwd[:fwd_start]+fwd_aligned, rev_aligned+rev[rev_start+len(rev_aligned.replace('-','')):])
        print('Merged  ',' '*max((rev_start-fwd_start),0)+overlapped)
        print()
        
        print('Aligned to Reference:')
        ref_name,ref_seq = parse_fasta('/Users/wakefield/Software/Allocate/brca1.fasta').next()
        ref_name = ' '.join(ref_name.split( )[:2])
        ref_aligned, query_aligned, ref_start, query_start = smith_waterman(ref_seq,overlapped)
        
        print(format_alignment([(ref_name,ref_seq[:ref_start]+ref_aligned+ref_seq[ref_start+len(ref_aligned.replace('-','')):]),('query',' '*ref_start+query_aligned+' '*(len(ref_seq)-(ref_start+len(ref_aligned.replace('-','')))))]))
        
        
        print()
    
    
    pass

#test data around region BRCA1 c.3481_3491del AGTATCTTCCT/-
#
#>17 dna:chromosome chromosome:GRCh37:17:41244007:41244117:1 
#CGCTTTTGCTAAAAACAGCAGAACTTTCCTTAATGTCATTTTCAGCAAAACTagtatcttcctTTATTTCACCATCATCTAACAGGTCATCAGGTGTCTCAGAACAAACCT
#
## Test data from HaloSeq FFPE run.  First 20 reads including BRCA1 c.3481 and their reverse

if __name__ == '__main__':
    args = process_commandline_args()
    args = {'forward_file':'/Users/wakefield/Software/Allocate/testdata_halo_FFPE_BRCA1_3481_R1.fastq',
            'reverse_file':'/Users/wakefield/Software/Allocate/testdata_halo_FFPE_BRCA1_3481_R2.fastq'}
    main(**args)

