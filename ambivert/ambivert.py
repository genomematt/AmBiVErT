#!/usr/bin/env python
# encoding: utf-8
"""
ambivert.py

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
import sys, os
import itertools, difflib, argparse
import hashlib, cPickle
from collections import defaultdict
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

logfile = sys.stderr

def smith_waterman(seq1,seq2):
    #this currently uses the clunky pycogent pure python implementation
    # slow but readable and trustworthy.
    #should be replaced with something else to 
    # - provide better gap extension penalties
    # - faster
    # - return start and end of alignment without this hack.
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
        print('Merging overlaps for {0} Pairs'.format(total),file=logfile)
        for key in self.get_above_threshold(threshold):
            print(completed,end=' ',file=logfile)
            fwd = self.data[key][0][0][1]
            rev = reverse_complement(self.data[key][0][1][1])
            fwd_aligned, rev_aligned, fwd_start, rev_start = smith_waterman(fwd,rev)
            if len(fwd_aligned) < minimum_overlap:
                self.unmergable.append(key)
            else:
                self.merged[key] = fwd[:fwd_start]+flatten_paired_alignment(fwd_aligned,rev_aligned)+rev[rev_start+len(rev_aligned.replace('-','')):]
            completed += 1
        print("\nSuccessfully merged {0} reads. {1} reads could not be merged".format(len(self.merged),len(self.unmergable)),file=logfile)
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
    
    def match_to_reference(self, min_score = 0.1, trim_primers=15):
        def match_by_edit(merged_key,ref_keys):
            best_score = 0
            best_hit = ''
            for ref_key in ref_keys:
                if trim_primers:
                    #if smith waterman was fast enough could use here instead of diflib to reduce gappy misassignments
                    score = difflib.SequenceMatcher(None, self.merged[merged_key][trim_primers:-trim_primers], self.reference_sequences[ref_key]).ratio()
                else:
                    score = difflib.SequenceMatcher(None, self.merged[merged_key], self.reference_sequences[ref_key]).ratio()
                if score > best_score:
                    best_score = score
                    best_hit = ref_key
            return best_hit,score
        
        print('Matching read bins to references',file=logfile)
        for merged_key in self.merged:
            if merged_key in self.reference:
                print(',',end='',file=logfile)
                continue
            best_hit,best_score = match_by_edit(merged_key,self.reference_sequences)
            if best_hit and best_score > min_score:
                self.reference[merged_key]= best_hit
                print('.',end='',file=logfile)
                #print("Matched", self.merged[merged_key], 'to', best_hit)
                #print(self.reference_sequences[best_hit])
            else:
                print("\nWARNING NO MATCH FOR ", self.merged[merged_key], file=logfile)
        print(file=logfile)        
        pass
    
    def align_to_reference(self):
        for merged_key in self.reference: #use reference to restrict to matched merged pairs
            if not merged_key in self.merged:
                continue #occurs when matched in cached file but not present in data set
            aligned_ref_seq, aligned_sample_seq = sw_align(self.reference_sequences[self.reference[merged_key]],self.merged[merged_key])
            print(self.reference[merged_key], file=logfile)
            if aligned_ref_seq != aligned_sample_seq:
                self.potential_variants.append(merged_key)
                print(aligned_ref_seq, file=logfile)
                print(aligned_sample_seq, file=logfile)
            self.aligned[merged_key] = (aligned_ref_seq, aligned_sample_seq)
        pass
    
    def get_amplicon_counts(self):
        amplicon_counts = {key:0 for key in self.reference_sequences}
        if not self.reference:
            self.align_to_reference()
        for merged_key in self.reference:
            counts = len(self.data[merged_key])
            if counts > self.threshold:
                amplicon_counts[self.reference[merged_key]] += len(self.data[merged_key])
        return amplicon_counts
    
    def save_hash_table(self,newhashfile):
        reference_sha224 = hashlib.sha224(repr(self.reference_sequences)).hexdigest()
        with newhashfile as outfile:
            cPickle.dump((reference_sha224,self.reference),outfile)
        pass
    
    def load_hash_table(self,hashfile):
        with hashfile as infile:
            reference_sha224,refdict = cPickle.load(infile)
            if reference_sha224 == hashlib.sha224(repr(self.reference_sequences)).hexdigest():
                self.reference = refdict
            else:
                print('WARNING: loaded read to reference hash library does not match reference sequences\n\
                    I really hope you know what you are doing... Check and if in doubt use --newhashfile\n \
                    without specifying an existing hash file.',file=logfile)
        pass
    
def process_commandline_args(*args,**kw):
    parser = argparse.ArgumentParser(description="""AmBiVErT: A program for binned analysis of amplicon data
        AmBiVErT clusters identical amplicon sequences and thresholds based on read frequency to remove technical errors.
        Due to sequencing errors occuring with a more random distribution than low frequency variants this approach
        reduces the number of amplicon products that must be assigned to target regions & assessed for variant calls.
        AmBiVErT overlaps forward and reverse reads from the same amplicon and preserves local phasing information.
        Typical running time for first use is several hours, which reduces to less than 10 minutes when the
        hash table calculated on a previous run is supplied for analysis of subsequent samples with the same amplicons.""")
    parser.add_argument('-f','--forward',
                        type=str,
                        help='a fastq format file of forward direction amplicon reads. \
                             May be compressed with gzip or bzip2 with the appropriate suffix (.gz/.bz2)')
    parser.add_argument('-r','--reverse',
                        type=str,
                        help='a fastq format file of reverse direction amplicon reads. \
                             May be compressed with gzip or bzip2 with the appropriate suffix (.gz/.bz2)')
    parser.add_argument('-m','--manifest',
                        type=argparse.FileType('U'),
                        help='an Illumina TrueSeq Amplicon manifest file.')
    parser.add_argument('--fasta',
                        type=argparse.FileType('U'),
                        help='an fasta amplicon manifest file. \
                            Sequences should be limited to the regions to be called and exclude primers & adaptors. \
                            This file can be provided in addition to an Illumina manifest to specify additional off target regions.')
    parser.add_argument('--output',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='output of alignments with variants. Default: stdout')    
    parser.add_argument('--countfile',
                        type=argparse.FileType('w'),
                        help='output of occurance counts per amplicon.  Includes all counts that map to the reference amplicon.\
                             This count does not include reads that occured at frequencies below --threshold <default=20>')    
    #parser.add_argument('--mpileup',
    #                    type=argparse.FileType('w'),
    #                    help='output mpileup files with filenames of the format "<name_provided_as_argument>_<amplicon_name>.mpileup"')    
    parser.add_argument('--threshold',
                        type=int,
                        default=20,
                        help='the minimum occurance threshold.  Amplicons that occur fewer than threshold times are ignored.')    
    parser.add_argument('--overlap',
                        type=int,
                        default=20,
                        help='The minimum overlap required between forward and reverse sequences to merge')    
    parser.add_argument('--primer',
                        type=int,
                        default=15,
                        help='The size of the smallest primer.  This number of bases is trimmed from the end of the merged sequences \
                             to reduce the possibility that small amplicons will fail to match due to primer mismatch')    
    parser.add_argument('--hashtable',
                        type=argparse.FileType('U'),
                        help='Filename for a precomputed hash table that matches amplicons to references.  Generate with --savehashtable')    
    parser.add_argument('--savehashtable',
                        type=argparse.FileType('w'),
                        help='Output a precomputed hash table that matches amplicons to references.  Use to speed up matching with --hashtable')    
    return parser.parse_args(*args,**kw)

def process_amplicon_data(forward_file, reverse_file,
                          manifest=None, fasta=None,
                          threshold=50, overlap=20, primer=15, 
                          savehashtable=None, hashtable=None,
                          ):
    amplicons = AmpliconData()
    amplicons.process_twofile_readpairs(forward_file, reverse_file)
    amplicons.merge_overlaps(minimum_overlap=20, threshold=threshold)
    if manifest:
        amplicons.add_references_from_manifest(manifest)
    if fasta:
        amplicons.add_references_from_fasta(fasta)
    if hashtable:
        amplicons.load_hash_table(hashtable)
    amplicons.match_to_reference()
    if savehashtable:
        amplicons.save_hash_table(savehashtable)
    amplicons.align_to_reference()
    return amplicons

if __name__ == '__main__':
    testargs = ['--forward','/Users/wakefield/Software/ABBA/130503_M00267_0038_L001_MS0302_12015B_E0023_TruSeq_R1.fastq',
                '--reverse','/Users/wakefield/Software/ABBA/130503_M00267_0038_L001_MS0302_12015B_E0023_TruSeq_R2.fastq',
                '--manifest','/Users/wakefield/Software/ABBA/TruSeq_CAT_Manifest_Allocate-CAT.txt',
                #'--hashtable','/Users/wakefield/Software/ABBA/cached_read_to_reference_table',
                '--savehashtable','/Users/wakefield/Software/ABBA/new_cached_read_to_reference_table',
                '--threshold','20',
                '--overlap','20',
                '--primer','15',
                '--countfile','/Users/wakefield/Software/ABBA/130503_M00267_0038_L001_MS0302_12015B_E0023_TruSeq_AmpliconCounts',
                ]
    args = process_commandline_args(testargs)
    #args = process_commandline_args()
    amplicons = process_amplicon_data(args.forward,args.reverse,
                                      args.manifest,args.fasta,
                                      args.threshold,args.overlap,args.primer,
                                      args.savehashtable,args.hashtable,
                                      )
    if args.countfile:
        with args.countfile as outfile:
            amplicon_counts = amplicons.get_amplicon_counts()
            for key in sorted(amplicon_counts.keys()):
                outfile.write('{0}\t{1}\n'.format(key,amplicon_counts[key]))
    
    for key in amplicons.potential_variants:
            aligned_ref_seq, aligned_sample_seq = amplicons.aligned[key]
            print(amplicons.reference[key],file=args.output)
            print(aligned_ref_seq,file=args.output)
            print(aligned_sample_seq,file=args.output)
            matches = ''
            for a,b in itertools.izip(aligned_ref_seq, aligned_sample_seq):
                if a == b:
                    matches += '.'
                else:
                    matches += b
            print(matches,file=args.output)
        
