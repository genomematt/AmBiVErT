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
import difflib
from collections import defaultdict, Counter
from cogent.align.algorithm import nw_align, sw_align
from sequence_utilities import parse_fastq, parse_fasta, reverse_complement, flatten_paired_alignment, format_alignment
from trueseq import parse_trueseq_manifest

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
        These arguments are for use if technical variation occurs at ends
        - trim5 : number of bases to trim from 5' for key default = None
        - trim3   : number of bases to trim from 3' for key default = None
    """
    def __init__(self, trim5=None, trim3=None):
        self.data = defaultdict(list)
        self.merged = {} # a dictionary of overlapped merged sequences with same key as self.data
        self.unmergable = [] # a list of keys (to self.data) for which merging failed
        self.aligned = {} # a dictionary of aligned (ref,sample) sequence tuples with same key as self.data
        self.location = {} # dictionary of (chromosome, start, end) with same key as self.data
        self.reference = {} # dictionary of reference sequence keys with same key as self.data
        self.reference_sequences = {} # dictionary of reference sequences
        self.potential_variants = [] # list of keys (to self.data) for which putative variants were identified 
        self.trim5 = trim5
        self.trim3 = None if trim3==None else -1*abs(trim3)
        self.threshold = 0
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
    
    def merge_overlaps(self, threshold=0, minimum_overlap=10):
        self.threshold = threshold #record in object data
        total = len(self.get_above_threshold(threshold))
        completed = 0
        print('Merging overlaps for {0} Pairs'.format(total),file=sys.stderr)
        for key in self.get_above_threshold(threshold):
            print(completed,end=' ',file=sys.stderr)
            fwd = self.data[key][0][0][1]
            rev = reverse_complement(self.data[key][0][1][1])
            fwd_aligned, rev_aligned, fwd_start, rev_start = smith_waterman(fwd,rev)
            #print(fwd,rev,fwd_aligned, rev_aligned, fwd_start, rev_start, sep='\n')
            if len(fwd_aligned) < minimum_overlap:
                self.unmergable.append(key)
            else:
                self.merged[key] = fwd[:fwd_start]+flatten_paired_alignment(fwd_aligned,rev_aligned)+rev[rev_start+len(rev_aligned.replace('-','')):]
                #print(self.merged[key])
            completed += 1
        print("\nSuccessfully merged {0} reads. {1} reads could not be merged".format(len(self.merged),len(self.unmergable)),file=sys.stderr)
        pass
    
    def add_references_from_fasta(self, fastafile):
        self.reference_sequences = {(x[0].split()[0],'unknown',1):x[1] for x in parse_fasta(fastafile)}
        pass
    
    def add_references_from_manifest(self, manifestfile):
        header, probes, targets = parse_trueseq_manifest(manifestfile)
        for target in targets:
            key = (target.TargetA.split()[0],target.Chromosome,target.Start_Position,target.End_Position)
            self.reference_sequences[key] = target.Sequence
        pass
    
    def match_to_reference(self, min_score = 0.8):
        print('Matching read bins to references',file=sys.stderr)
        for merged_key in self.merged:
            best_score = 0
            best_hit = ''
            for ref_key in self.reference_sequences:
                score = difflib.SequenceMatcher(None, self.merged[merged_key], self.reference_sequences[ref_key]).ratio()
                #print(ref_key, score, best_hit, best_score)
                if score > best_score and score >= min_score:
                    best_score = score
                    best_hit = ref_key
            if best_hit:
                self.reference[merged_key]= best_hit
                print('.',end='',file=sys.stderr)
                #print("Matched", self.merged[merged_key], 'to', best_hit)
                #print(self.reference_sequences[best_hit])
            else:
                print("\nWARNING NO MATCH FOR ", self.merged[merged_key], file=sys.stderr)
        print(file=sys.stderr)        
        pass
    
    def align_to_reference(self):
        for merged_key in self.reference:
            aligned_ref_seq, aligned_sample_seq = sw_align(self.reference_sequences[self.reference[merged_key]],self.merged[merged_key])
            print(self.reference[merged_key], file=sys.stderr)
            if aligned_ref_seq != aligned_sample_seq:
                self.potential_variants.append(merged_key)
                print(aligned_ref_seq, file=sys.stderr)
                print(aligned_sample_seq, file=sys.stderr)
            self.aligned[merged_key] = (aligned_ref_seq, aligned_sample_seq)
        pass
    
    def get_amplicon_counts(self):
        amplicon_counts = {key:0 for key in self.reference_sequences}
        if not self.reference:
            self.align_to_reference()
        for merged_key in self.reference:
            amplicon_counts[self.reference[merged_key]] += len(self.data[merged_key])
        return amplicon_counts
            
    

def process_commandline_args():
    pass

def process_amplicon_data(forward_file, reverse_file, threshold=1000):
    amplicons = AmpliconData()
    amplicons.process_twofile_readpairs(forward_file, reverse_file)
    amplicons.merge_overlaps(minimum_overlap=20, threshold=threshold)
    #amplicons.add_references_from_fasta(fastafile='/Users/wakefield/Software/Allocate/brca1.fasta')
    amplicons.add_references_from_manifest(open('/Users/wakefield/Software/abba/TruSeq_CAT_Manifest_Allocate-CAT.txt','U'))
    #amplicons.add_references_from_manifest(open('/Users/wakefield/Software/abba/TruSeq_CAT_Manifest_TC0019069-CAT.txt','U'))
    #amplicons.add_references_from_manifest(open('/Users/wakefield/Software/abba/TruSeq_CAT_Manifest_TC0019072-CAT.txt','U'))
    amplicons.match_to_reference()
    amplicons.align_to_reference()
    #amplicons.to_mpileup(filename=sys.stdout)
    print(amplicons.get_amplicon_counts())
    
    #print(amplicons)
    return amplicons
    

#test data around region BRCA1 c.3481_3491del AGTATCTTCCT/-
#
#>17 dna:chromosome chromosome:GRCh37:17:41244007:41244117:1 
#CGCTTTTGCTAAAAACAGCAGAACTTTCCTTAATGTCATTTTCAGCAAAACTagtatcttcctTTATTTCACCATCATCTAACAGGTCATCAGGTGTCTCAGAACAAACCT
#
## Test data from HaloSeq FFPE run.  First 20 reads including BRCA1 c.3481 and their reverse

if __name__ == '__main__':
    args = process_commandline_args()
    #args = {'forward_file':'/Users/wakefield/Software/Allocate/testdata_halo_FFPE_BRCA1_3481_R1.fastq',
    #        'reverse_file':'/Users/wakefield/Software/Allocate/testdata_halo_FFPE_BRCA1_3481_R2.fastq'}
    
    args = {'forward_file':'/Users/wakefield/Software/ABBA/130503_M00267_0038_L001_MS0302_12015B_E0023_TruSeq_R1.fastq',
            'reverse_file':'/Users/wakefield/Software/ABBA/130503_M00267_0038_L001_MS0302_12015B_E0023_TruSeq_R2.fastq'}
    
    #print(repr({(x[0].split()[0],'unknown',1):x[1] for x in parse_fasta('/Users/wakefield/Software/Allocate/brca1.fasta')}))
    amplicons = process_amplicon_data(**args)
    #print(amplicons.reference_sequences)

