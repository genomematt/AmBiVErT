#!/usr/bin/env python3
# encoding: utf-8
"""
ambivert.py

    ___              ____  _ _    ________    ______
   /   |  ____ ___  / __ )(_) |  / / ____/___/_  __/
  / /| | / __ `__ \/ __  / /| | / / __/ / ___// /   
 / ___ |/ / / / / / /_/ / / | |/ / /___/ /   / /    
/_/  |_/_/ /_/ /_/_____/_/  |___/_____/_/   /_/     
                                                    
AmBiVErT: A program for binned analysis of amplicon data

AMplicon BInning Variant-caller with ERror Truncation

This program is designed for the processing of amplicon based resequencing data
generated on second generation sequencing platforms (eg Illumina HiSeq/MiSeq).

The approach used is to batch reads derived from each amplicon together in a
clustering step prior to aligning the sequence to a reference genome.

Steps:
    1)  Cluster similar reads
    2)  Align for overlap
    3)  Align to reference target sequences
    4)  Consolidate calls
    5)  Output VCF format calls

It is specifically designed not to share any of its codebase with other variant
calling pipelines or software, run quickly and minimize false positives.

Our intended purpose for this software is as lightweight backup and quality assurance
complementing a more traditional variant calling pipeline, and to provide amplicon level
phasing of multiple mutations that exist on the same amplicon.

All lines of code are covered by unit tests unless marked with #pragma no cover

Created by Matthew Wakefield and Graham Taylor.
Copyright (c) 2013-2015  Matthew Wakefield and The University of Melbourne. All rights reserved.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

"""
# ambivert noun PSYCHOLOGY
# a person who has a balance of extrovert and introvert features in their personality

import sys, os, io
import warnings, logging
import itertools, difflib, argparse
import hashlib, pickle
from collections import defaultdict
from multiprocessing.dummy import Pool #Threaded rather than process version of multiprocessing
import ambivert.align
from ambivert.sequence_utilities import *
from ambivert.truseq_manifest import make_sequences, parse_truseq_manifest
from ambivert.call_variants import call_variants, call_variants_to_vcf, make_vcf_header

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2013-2015,  Matthew Wakefield and The University of Melbourne"
__credits__ = ["Matthew Wakefield","Graham Taylor","Toby Sargeant"]
__license__ = "GPLv3"
__version__ = "0.5.2"
__maintainer__ = "Matthew Wakefield"
__email__ = "matthew.wakefield@unimelb.edu.au"
__status__ = "Development"

logfile = sys.stderr
logging.basicConfig(format='%(levelname)s:AmBiVErT:%(message)s',level=logging.DEBUG)
logging.disable(logging.DEBUG)
logging.captureWarnings(True)


def smith_waterman(seq1,seq2):
    """Wrapper method for Smith-Waterman alignment using
    the plumb.bob C implementation
    
    Arguments:
        - seq1 : an ascii DNA sequence string.  This is the query
                 sequence and must be all upper case
        - seq2 : an ascii DNA sequence string.  This is the reference
                 sequence and may contain lower case soft masking """
    alignment =  ambivert.align.local_align(bytes(seq1,'ascii'), len(seq1),
                                bytes(seq2.upper(),'ascii'), len(seq2),
                                ambivert.align.DNA_MAP[0],
                                ambivert.align.DNA_MAP[1], 
                                ambivert.align.DNA_SCORE,
                                -7, -1 #gap open, gap extend
                                )
    if '-' in seq1 or '-' in seq2:
        raise RuntimeError('Aligning Sequences with gaps is not supported',seq1,seq2)
    start_seq1 = alignment.contents.align_frag.contents.sa_start
    start_seq2 = alignment.contents.align_frag.contents.sb_start
    frag = alignment[0].align_frag
    align_seq1 = ''
    align_seq2 = ''
    while frag:
        frag = frag[0]
        if frag.type == ambivert.align.MATCH:
            f1 = seq1[frag.sa_start:frag.sa_start + frag.hsp_len]
            f2 = seq2[frag.sb_start:frag.sb_start + frag.hsp_len]
            align_seq1 += f1
            align_seq2 += f2
        elif frag.type == ambivert.align.A_GAP:
            align_seq1 += '-' * frag.hsp_len
            align_seq2 += seq2[frag.sb_start:frag.sb_start + frag.hsp_len]
        elif frag.type == ambivert.align.B_GAP:
            align_seq1 += seq1[frag.sa_start:frag.sa_start + frag.hsp_len]
            align_seq2 += '-' * frag.hsp_len
        frag = frag.next
    assert len(align_seq1) == len(align_seq2)
    ambivert.align.alignment_free(alignment)
    return align_seq1,align_seq2,start_seq1,start_seq2

def needleman_wunsch(seq1,seq2):
    """Wrapper method for Needleman-Wunsch alignment using
    the plumb.bob C implementation
    
    Arguments:
        - seq1 : an ascii DNA sequence string.  This is the query
                 sequence and must be all upper case
        - seq2 : an ascii DNA sequence string.  This is the reference
                 sequence and may contain lower case soft masking """
    alignment =  ambivert.align.global_align(bytes(seq1,'ascii'), len(seq1),
                                bytes(seq2.upper(),'ascii'), len(seq2),
                                ambivert.align.DNA_MAP[0],
                                ambivert.align.DNA_MAP[1], 
                                ambivert.align.DNA_SCORE,
                                -7, -1 #gap open, gap extend
                                )
    if '-' in seq1 or '-' in seq2:
        raise RuntimeError('Aligning Sequences with gaps is not supported',seq1,seq2)
    start_seq1 = 0
    start_seq2 = 0
    frag = alignment[0].align_frag
    align_seq1 = ''
    align_seq2 = ''
    while frag:
        frag = frag[0]
        if frag.type == ambivert.align.MATCH:
            f1 = seq1[frag.sa_start:frag.sa_start + frag.hsp_len]
            f2 = seq2[frag.sb_start:frag.sb_start + frag.hsp_len]
            align_seq1 += f1
            align_seq2 += f2
        elif frag.type == ambivert.align.A_GAP:
            align_seq1 += '-' * frag.hsp_len
            align_seq2 += seq2[frag.sb_start:frag.sb_start + frag.hsp_len]
        elif frag.type == ambivert.align.B_GAP:
            align_seq1 += seq1[frag.sa_start:frag.sa_start + frag.hsp_len]
            align_seq2 += '-' * frag.hsp_len
        frag = frag.next
    assert len(align_seq1) == len(align_seq2)
    ambivert.align.alignment_free(alignment)
    return align_seq1,align_seq2,start_seq1,start_seq2

def calculate_global_score(args):
    """Calculate alignment score for amplicon to reference
    comparison using global Needleman-Wunsch alignment.
    Argument:
        a tuple of (ref_key, query_sequence, reference_sequence)
        The arguments are constructed as tuple to support use 
        with the map function allowing pool.map for concurrency
    Returns:
        a tuple of (ref_key,score)
    """
    ref_key, seq1, reference_sequence = args
    alignment =  ambivert.align.global_align(bytes(seq1, 'ascii'), len(seq1),
                    bytes(reference_sequence.upper(),'ascii'), len(reference_sequence),
                    ambivert.align.DNA_MAP[0],
                    ambivert.align.DNA_MAP[1], 
                    ambivert.align.DNA_SCORE,
                    -7, -1 #gap open, gap extend
                    )
    score = int(alignment[0].score)
    ambivert.align.alignment_free(alignment)
    return (ref_key,score)

def calculate_local_score(args):
    """Calculate alignment score for amplicon to reference
    comparison using local Smith-Waterman alignment.
    Argument:
        a tuple of (ref_key, query_sequence, reference_sequence)
        The arguments are constructed as tuple to support use 
        with the map function allowing pool.map for concurrency
    Returns:
        a tuple of (ref_key,score)
    """
    ref_key, seq1, reference_sequence = args
    alignment =  ambivert.align.global_align(bytes(seq1, 'ascii'), len(seq1),
                    bytes(reference_sequence.upper(),'ascii'), len(reference_sequence),
                    ambivert.align.DNA_MAP[0],
                    ambivert.align.DNA_MAP[1], 
                    ambivert.align.DNA_SCORE,
                    -7, -1 #gap open, gap extend
                    )
    score = int(alignment[0].score)
    ambivert.align.alignment_free(alignment)
    return (ref_key,score)


class AmpliconData(object):
    """A Class for holding read data from amplicon experiments
    Reads are indexed by hashes of both the forward and reverse reads
    Arguments:
        These arguments are for use if technical variation occurs at ends
        - trim5 : number of bases to trim from 5' for key default = None
        - trim3   : number of bases to trim from 3' for key default = None
    """
    def __init__(self, trim5=None, trim3=None):
        self.data = defaultdict(list) #{md5hexdigest:[((f_name, f_seq, f_qual),(r_name, r_seq, r_qual)),]}
        self.merged = {} # a dictionary of overlapped merged sequences with same key as self.data
        self.unmergable = [] # a list of keys (to self.data) for which merging failed
        self.aligned = {} # a dictionary of aligned (ref,sample) sequence tuples with same key as self.data
        self.location = {} # dictionary of (chromosome, start, end, strand) with same key as self.data
        self.reference = {} # dictionary of reference sequence keys with same key as self.data
        self.reference_sequences = {} # dictionary of reference sequences. Key is (name, chromosome, start, end, strand)
        self.potential_variants = [] # list of keys (to self.data) for which putative variants were identified
        self.called_variants = {} # dictionary of lists of called variants with same key as self.data
        self.consolidated_variants = []
        self.trim5 = trim5
        self.trim3 = None if trim3==None else -1*abs(trim3)
        self.threshold = 0
        self.readpairs = 0
        self.input_filenames = []
        self.reference_filenames = []
        pass
    
    def __str__(self):
        """Returns a description of object content and settings"""
        result = ""
        result += "Data file names: \n"
        for filename in self.input_filenames:
            result += filename + '\n'
        result += "Reference file names: \n"
        for filename in self.reference_filenames:
            result += filename + '\n'
        result += "Readpairs: {0}\n".format(self.readpairs)
        result += "Threshold: {0}\n".format(self.threshold)
        result += "Unique read groups: {0}\n".format(len(self.merged))
        result += "Unmerged pairs: {0}\n".format(len(self.unmergable))
        return result
    
    def add_reads(self, f_name, f_seq, f_qual,r_name, r_seq, r_qual):
        """Add read sequence pairs to amplicon object representation
        Reads are stored internally as tuples of tuples in a dictionary of lists.
        The dictionary is keyed by the 16 character MD5 hexdigest of the concatinated
        forward and reverse sequences.
        ie  {md5hexdigest:[((f_name, f_seq, f_qual),(r_name, r_seq, r_qual)),]}
        
        Arguments:
            f_name  :   forward read name (ascii string)
            f_seq   :   forward sequence (ascii string)
            f_qual  :   forward quality scores (ascii string)
            r_name  :   reverse read name (ascii string)
            r_seq   :   reverse sequence (ascii string)
            r_qual  :   reverse quality scores (ascii string)
        """
        amplicon_key = hashlib.md5(bytes(f_seq[self.trim5:self.trim3]+r_seq[self.trim5:self.trim3], 'ascii')).hexdigest()
        self.data[amplicon_key].append(((f_name, f_seq, f_qual),(r_name, r_seq, r_qual)))
        self.readpairs += 1
        pass
    
    def process_twofile_readpairs(self, forward_file, reverse_file, parser=parse_fastq):
        """Parse two file fastq format read pairs and add readpairs to amplicon object
        by calling AmpliconData.add_reads
        
        Arguments:
            forwardfile : a file object containing fastq formatted forward reads
            reversefile : a file object containing matched reverse reads
            parser      : an parser function that yields (name,sequence,quality) from
                          a fastqfile.  Default ambivert.sequence_utilities.parse_fastq
        """
        logging.info('Reading data...')
        #sys.stderr.flush()
        if hasattr(forward_file,'name'):
            self.input_filenames.append(forward_file.name)
            logging.info('    Forward read file: {0}'.format(forward_file.name))
            #sys.stderr.flush()
        if hasattr(reverse_file,'name'):
            self.input_filenames.append(reverse_file.name)
            logging.info('    Reverse read file: {0}'.format(reverse_file.name))
            #sys.stderr.flush()
        for (f_name, f_seq, f_qual),(r_name, r_seq, r_qual) in zip(parser(forward_file), parser(reverse_file)):
            assert f_name.split()[0] == r_name.split()[0]
            self.add_reads(f_name, f_seq, f_qual,r_name, r_seq, r_qual)
        logging.info('Read {0} read pairs'.format(self.readpairs))
        #sys.stderr.flush()
        pass
    
    def get_above_threshold(self, threshold=1):
        """Returns groups of reads with identical forward and reverse
        sequences that occur more times than a threshold
        Arguments:
            threshold   : The number of occurences that must be exceeded
        Returns:
            a list of md5 hexdigest keys to read groups in AmpliconData.data
        """
        return [x for x in self.data if len(self.data[x]) > threshold]
    
    def merge_overlaps(self, threshold=0, minimum_overlap=10):
        """Uses Smith-Waterman alignment to overlap the ends of paired reads
        from AmpliconData.data creating AmpliconData.merged which maintains
        the same md5 hexdigest key as AmpliconData.data and a string of
        the merged sequence.
        Reads that do not overlap between forward and reverse are excluded and
        their md5 hexdigest key added to AmpliconData.unmergable
        Arguments:
            threshold   : The number of occurences that must be exceeded before
                          attempting to merge.  Default = 0
            minimum_overlap : The minimum number of bases that must overlap
                          for a merge to occur
        """
        #uses reverse_complement & flatten_paired_alignment from ambivert.sequence_utilities
        self.threshold = threshold #record in object data
        total = len(self.get_above_threshold(threshold))
        completed = 0
        logging.info('Merging overlaps for {0} unique pairs'.format(total))
        #sys.stderr.flush()
        for key in self.get_above_threshold(threshold):
            fwd = self.data[key][0][0][1]
            rev = reverse_complement(self.data[key][0][1][1])
            fwd_aligned, rev_aligned, fwd_start, rev_start = smith_waterman(fwd,rev)
            if len(fwd_aligned) < minimum_overlap:
                self.unmergable.append(key)
            else:
                self.merged[key] = fwd[:fwd_start]+flatten_paired_alignment(fwd_aligned,rev_aligned)+rev[rev_start+len(rev_aligned.replace('-','')):]
            completed += 1
        logging.info("\nSuccessfully merged {0} reads. {1} reads could not be merged".format(len(self.merged),len(self.unmergable)))
        #sys.stderr.flush()
        pass
    
    def add_references_from_fasta(self, fastafile):
        """Add reference sequences from a fasta file
        Reference sequences should include the entire amplicon
        sequence.  Regions that should not be called (eg primers)
        should be soft masked in lower case.  All other bases
        should be upper case.
        Additional annotation details are required on the name line
        and should be tab separated in the format:
        >amplicon_name    chromosome    start    end    strand(+/-)
        
        Arguments:
            fastafile : a file object containg fasta entries
            
        Data stored in AmpliconData.reference_sequences with key
        (name, chromosome, start, end, strand)
        """
        for name,sequence in parse_fasta(fastafile):
            self.reference_sequences[tuple(name.split())] = sequence
        pass
    
    def add_references_from_manifest(self, manifestfile):
        """Add reference sequences from an Illumina TruSeq amplicon
        manifest file.  Lightly tested for impact on variation in the
        implementation.
        Currently strictly asserts for presence of expected column header
        even when these columns are not required
        See ambivert.truseq_manifest for details
        """
        for name,sequence in make_sequences(*parse_truseq_manifest(manifestfile), with_probes=True, softmask_probes=True, all_plus=False):
            self.reference_sequences[tuple(name.split())] = sequence
        pass
    
    def match_to_reference(self, min_score=0.1, trim_primers=0, global_align=True, show_progress=False, threads=1):
        """Assign read groups to reference sequences by competitive alignment
        Does not do reverse complementing or alignment formatting and relies
        only on score to match
        Checks AmpliconData.reference for cached results to avoid
        redundant calculation
        Stores results in AmpliconData.reference keyed by md5 hexdigest in
        as keys to AmpliconData.reference_sequences in format 
        (name, chromosome, start, end, strand)
        
        Arguments:
            min_score :         The minimum alignment score to be considered a match
                                Default = 0.1
            trim_primers:       A fixed size to be trimmed prior to matching
                                Default = 0
            global_alignment:   Use global Needleman-Wunsch alignment instead of
                                Smith-Waterman.  Default = True
        """
        def match_by_edit(merged_key,ref_keys): #pragma no cover
            """original matching method using edit distances
            inaccurate handeling of indels results in frequent
            misassignment.
            DEPRECIATED - DO NOT USE"""
            #depreciated
            best_score = 0
            best_hit = ''
            for ref_key in ref_keys:
                if trim_primers:
                    score = difflib.SequenceMatcher(None, self.merged[merged_key][trim_primers:-trim_primers], self.reference_sequences[ref_key]).ratio()
                else:
                    score = difflib.SequenceMatcher(None, self.merged[merged_key], self.reference_sequences[ref_key]).ratio()
                if score > best_score:
                    best_score = score
                    best_hit = ref_key
            return best_hit,best_score
        
        def match_by_smith_waterman(merged_key,ref_keys):#pragma no cover
            """Match merged amplicon to reference by local alignment
            Arguments:
                merged_key : a md5 hexdigest key to AmpliconData.merged
                ref_keys   : a list of keys to AmpliconData.reference_sequences
                             in format [(name, chromosome, start, end, strand),]
            DEPRECIATED - USE PARALLEL VERSION
            """
            #DEPRECIATED
            best_score = 0
            best_hit = ''
            for ref_key in ref_keys:
                if trim_primers:
                    seq1 = self.merged[merged_key][trim_primers:-trim_primers]
                else:
                    seq1 = self.merged[merged_key]
                alignment =  ambivert.align.local_align(bytes(seq1, 'ascii'), len(seq1),
                                bytes(self.reference_sequences[ref_key].upper(),'ascii'), len(self.reference_sequences[ref_key]),
                                ambivert.align.DNA_MAP[0],
                                ambivert.align.DNA_MAP[1], 
                                ambivert.align.DNA_SCORE,
                                -7, -1 #gap open, gap extend
                                )
                score = alignment[0].score
                ambivert.align.alignment_free(alignment)
                if score > best_score:
                    best_score = score
                    best_hit = ref_key
            return best_hit,best_score

        def match_by_needleman_wunsch(merged_key,ref_keys):#pragma no cover
            """Match merged amplicon to reference by global alignment
            Arguments:
                merged_key : a md5 hexdigest key to AmpliconData.merged
                ref_keys   : a list of keys to AmpliconData.reference_sequences
                             in format [(name, chromosome, start, end, strand),]
            DEPRECIATED - USE PARALLEL VERSION
            """
            #DEPRECIATED
            best_score = int(-1000000)
            best_hit = ''
            for ref_key in ref_keys:
                if trim_primers:
                    seq1 = self.merged[merged_key][trim_primers:-trim_primers]
                else:
                    seq1 = self.merged[merged_key]
                alignment =  ambivert.align.global_align(bytes(seq1, 'ascii'), len(seq1),
                                bytes(self.reference_sequences[ref_key].upper(),'ascii'), len(self.reference_sequences[ref_key]),
                                ambivert.align.DNA_MAP[0],
                                ambivert.align.DNA_MAP[1], 
                                ambivert.align.DNA_SCORE,
                                -7, -1 #gap open, gap extend
                                )
                score = int(alignment[0].score)
                ambivert.align.alignment_free(alignment)
                if score > best_score:
                    best_score = score
                    best_hit = ref_key
            return best_hit,best_score

        def parallel_match_by_needleman_wunsch(merged_key,ref_keys,threads=threads):
            """Match merged amplicon to reference by global alignment
            using pythons multiprocessing to match in parallel
            Arguments:
                merged_key : a md5 hexdigest key to AmpliconData.merged
                ref_keys   : a list of keys to AmpliconData.reference_sequences
                             in format [(name, chromosome, start, end, strand),]
            """
            if trim_primers:
                seq1 = self.merged[merged_key][trim_primers:-trim_primers]
            else:
                seq1 = self.merged[merged_key]
                        
            pool = Pool(threads)
            scores = pool.map(calculate_global_score, [(x, seq1, self.reference_sequences[x]) for x in ref_keys])
            pool.close()
            
            best_score = int(-1000000)
            best_hit = ''
            
            for ref_key, score in scores:
                if score > best_score:
                    best_score = score
                    best_hit = ref_key
            return best_hit,best_score
        
        def parallel_match_by_smith_waterman(merged_key,ref_keys,threads=threads):
            """Match merged amplicon to reference by global alignment
            using pythons multiprocessing to match in parallel
            Arguments:
                merged_key : a md5 hexdigest key to AmpliconData.merged
                ref_keys   : a list of keys to AmpliconData.reference_sequences
                             in format [(name, chromosome, start, end, strand),]
            """
            if trim_primers:
                seq1 = self.merged[merged_key][trim_primers:-trim_primers]
            else:
                seq1 = self.merged[merged_key]
                        
            pool = Pool(threads)
            scores = pool.map(calculate_local_score, [(x, seq1, self.reference_sequences[x]) for x in ref_keys])
            pool.close()
            
            best_score = 0
            best_hit = ''
            
            for ref_key, score in scores:
                if score > best_score:
                    best_score = score
                    best_hit = ref_key
            return best_hit,best_score
        
        logging.info('Matching read bins to references - commas represent cached matches, periods matching by alignment')
        #sys.stderr.flush()
        for merged_key in self.merged:
            if merged_key in self.reference:
                if show_progress: #pragma no cover
                    print(',',end='',file=logfile)
                    logfile.flush()
                logging.debug("Matched in cached hashtable {0} to \n\t{1}\n\t{2}".format(
                            self.merged[merged_key],self.reference[merged_key],
                            self.reference_sequences[self.reference[merged_key]]))
                continue
            if global_align:
                best_hit,best_score = parallel_match_by_needleman_wunsch(merged_key,self.reference_sequences, threads=threads)
            else:
                best_hit,best_score = parallel_match_by_smith_waterman(merged_key,self.reference_sequences, threads=threads)
            if best_hit and best_score > min_score:
                self.reference[merged_key]= best_hit
                if show_progress: #pragma no cover
                    print('.',end='',file=logfile)
                    logfile.flush()
                logging.debug("Matched {0} to \n\t{1}\n\t{2}".format(
                            self.merged[merged_key],best_hit,self.reference_sequences[best_hit]))
            else:
                logging.warning("\nWARNING NO MATCH FOR {0}\n\
                                Use ambivert --fastqamplicon {1} to retrieve reads from this amplicon\
                                ".format(self.merged[merged_key],merged_key))
        if show_progress:
            print(file=logfile)
        pass
    
    def align_to_reference(self, global_align=True):
        """Perform reverse complementation if necessary and 
        full alignment to produce left aligned plus strand
        reference and query sequences.
        Arguments:
            global_alignment:   Use global Needleman-Wunsch alignment instead of
                                Smith-Waterman.  Default = True
        
        Stores (aligned_sample_seq, aligned_ref_seq) in AmpliconData.aligned
        and (chromosome, int(start), int(end), ref_start) in AmpliconData.location
        ref_start will be 0 for global and alignment offset in reference for local
        keyed by the same md5 hexdigest key as AmpliconData.data and AmpliconData.merged
        """
        for merged_key in self.reference: #use reference to restrict to matched merged pairs
            if not merged_key in self.merged:
                continue #occurs when matched in cached file but not present in data set
            if self.reference[merged_key][-1] == '-': #last element of self.reference_sequence key indicates minus strand probe
                #both the reference sequence and the merged sequence are on the minus strand
                query_seq = reverse_complement(self.merged[merged_key])
                ref_seq = reverse_complement(self.reference_sequences[self.reference[merged_key]])
            else:
                query_seq = self.merged[merged_key]
                ref_seq = self.reference_sequences[self.reference[merged_key]]
            if global_align:
                aligned_sample_seq,aligned_ref_seq,sample_start,ref_start = needleman_wunsch(query_seq,ref_seq)
            else:
                aligned_sample_seq,aligned_ref_seq,sample_start,ref_start = smith_waterman(query_seq,ref_seq)
            if aligned_ref_seq.upper() != aligned_sample_seq:
                self.potential_variants.append(merged_key)
            self.aligned[merged_key] = (aligned_sample_seq, aligned_ref_seq) #plus strand
            self.location[merged_key] = (self.reference[merged_key][1], int(self.reference[merged_key][2]), int(self.reference[merged_key][3]), ref_start)
        pass
    
    def get_amplicon_count(self, key):
        """Return the number of times a specific amplicon sequence
        occurs in the data.
        Arguments:
            key:    An md5 hexdigest amplicon identifier
                    Note that these are unique for full sequence
                    of the amplicon including trimmed sequence.
        Returns:
            integer count of amplicons if greater than threshold.
            zero if less than threshold.
        """
        return len(self.data[key])
    
    def get_amplicon_counts(self):
        """Get a dictionary of amplicon sequence occurance counts
        for all amplicons that have been assigned to the reference.
        Returns:
            a dictionary of {key:occurance} where
            key = An md5 hexdigest amplicon identifier from self.reference
            occurance = integer count of amplicons if greater than threshold.
            zero if less than threshold.
        """
        amplicon_counts = {key:0 for key in self.reference_sequences}
        if not self.reference:
            logging.info('Amplicons were not aligned')
            self.align_to_reference()
        for merged_key in self.reference:
            counts = len(self.data[merged_key])
            if counts > self.threshold:
                amplicon_counts[self.reference[merged_key]] += len(self.data[merged_key])
        return amplicon_counts
    
    def save_hash_table(self,newhashfile):
        """Save a lookup table with precomputed assignment
        of amplicon sequences to reference locations.
        Only saves assignments that have no potential variants.
        Arguments:
            newhashfile:    a writable python3 binary file
        """
        reference_sha224 = hashlib.sha224(repr(sorted(self.reference_sequences)).encode('latin-1')).hexdigest()
        with newhashfile as outfile:
            hash_dictionary = {merged_key:self.reference[merged_key] for merged_key in self.reference if merged_key not in self.potential_variants} 
            pickle.dump((reference_sha224,hash_dictionary),outfile)
        pass
    
    def load_hash_table(self,hashfile):
        """Load a lookup table with precomputed assignment
        of amplicon sequences to reference locations.
        Arguments:
            newhashfile:    a readable python3 binary file
                            with a pickled data structure
                            from AmpliconData.save_hash_table
        """
        with hashfile as infile:
            reference_sha224,refdict = pickle.load(infile)
            if reference_sha224 == hashlib.sha224(repr(sorted(self.reference_sequences)).encode('latin-1')).hexdigest():
                self.reference = refdict
            else:
                warnings.warn('WARNING: loaded read to reference hash library does not match reference sequences\n'+ \
                    'I really hope you know what you are doing... Check and if in doubt use --savehashtable\n'+ \
                    'without specifying an existing hash file.',UserWarning)
        pass
    
    def print_variants_as_alignments(self, outfile=sys.stdout):
        """Print a human readable alignment of potentially
        variant containing amplicons.  These amplicons are
        pre variant calling and filtering and contain all
        amplicons that differ from the reference - including
        in softmasked regions.
        
        Arguments:
            outfile:    a writable python3 text file object
                        Default = sys.stdout
        """
        for key in sorted(self.potential_variants):
            aligned_sample_seq, aligned_ref_seq = self.aligned[key]
            print(key,file=outfile)
            print(self.reference[key],file=outfile)
            print(aligned_ref_seq,file=outfile)
            matches = ''
            for a,b in zip(aligned_ref_seq, aligned_sample_seq):
                if a == b or a in 'abcdghkmnrstuvwy': #softmasked or match
                    matches += '.'
                else:
                    matches += b
            print(matches,file=outfile)
            print(aligned_sample_seq,file=outfile)
            print(file=outfile)
        pass
    
    def call_amplicon_variants(self):
        """Call variants in amplicons
        Takes input from AmpliconData.potential_variants
        Outputs to AmpliconData.called_variants
        """
        # uses call_variants from ambivert.call_variants
        for key in self.potential_variants:
            aligned_sample_seq, aligned_ref_seq = self.aligned[key]
            name, chromosome, amplicon_position, end, strand = self.reference[key]
            self.called_variants[key] = call_variants(aligned_sample_seq, aligned_ref_seq, chromosome, int(amplicon_position))
        pass
    
    def get_variant_positions(self):
        """Get variants position in reference coordinates
        Returns:
            a sorted list of tuples of the form
            (chromosome, start_positions, end_position)
        """
        variant_positions = []
        for amplicon_id in self.called_variants:
            for variant in self.called_variants[amplicon_id]:
                variant_positions.append((variant.chromosome, variant.start, variant.end))
        return sorted(list(set(variant_positions)))
    
    def consolidate_variants(self, exclude_softmasked_coverage=True):
        """Consolidate and normalize variant calls
        Takes input from AmpliconData.called_variants
        Outputs to AmpliconData.consolidated_variants
        
        Arguments:
            exclude_softmasked_coverage: boolian flag to exclude
                                         variants that overlap softmasked
                                         (lower case in reference) sequence
                                         Default = True
        """
        positions = self.get_variant_positions()
        for (chrom, start, end) in positions:
            if exclude_softmasked_coverage:
                amplicon_ids = self.get_amplicons_overlapping_without_softmasking(chrom,start,end-start)
            else:
                amplicon_ids = self.get_amplicons_overlapping(chrom,start,end-start)
            ref = [] #list of amplicon keys
            alt = {} #dictionary of lists of amplicons keyed by alternative alleles
            for amplicon_id in amplicon_ids:
                overlapping_variant = False
                if amplicon_id in self.called_variants:
                    for variant in self.called_variants[amplicon_id]:
                        if not ((variant.end < start) or (variant.start > end)):
                            overlapping_variant = True
                            if variant in alt:
                                alt[variant].append(amplicon_id)
                            else:
                                alt[variant] = [amplicon_id,]
                if not overlapping_variant:
                    ref.append(amplicon_id)
            
            total_depth = sum([self.get_amplicon_count(amplicon_id) for amplicon_id in amplicon_ids])
            
            #indels are not described at the same coordinate as snps so we deal with these first
            #we dont do 'correct' VCF with all alt alleles on one line, we repeat positions instead.
            #in rare cases of two variants overlapping a deletion we just duplicate an identical deletion entry
            #exact duplicate records are then cleaned up at end of this function
            #this is structured to allow compound calling in a subsequent version
            indels = [variant for variant in alt.keys() if len(variant.alt_allele) > 1 or len(variant.ref_allele) > 1 ]
            for indel in indels:
                variant_depth = sum([self.get_amplicon_count(amplicon_id) for amplicon_id in alt[indel]])
                logging.debug(str((indel,variant_depth, total_depth, variant_depth/total_depth, ref, alt[indel])))
                self.consolidated_variants.append((indel,variant_depth, total_depth, variant_depth/total_depth, tuple(ref), tuple(alt[indel])))
            
            #for now we do one variant per line for snps
            #this should change to at least an optional correct VCF formatting of all alleles on one line
            snvs = [key for key in alt.keys() if len(key.alt_allele) == 1 and len(key.ref_allele) == 1 ]
            for snv in snvs:
                variant_depth = sum([self.get_amplicon_count(key) for key in alt[snv]])
                logging.debug(str((snv,variant_depth, total_depth, variant_depth/total_depth, ref, alt[snv])))
                self.consolidated_variants.append((snv,variant_depth, total_depth, variant_depth/total_depth, tuple(sorted(ref)), tuple(sorted(alt[snv]))))
        #logic above does not preclude calling the same variant more than once
        #so we remove identical records
        self.consolidated_variants = sorted(list(set(self.consolidated_variants)))
        pass
    
    def get_filtered_variants(self, min_cover=0, min_reads=0, min_freq=0.1):
        """Get variants that match filtering conditions
        Arguments:
            min_cover:  minimum total coverage of all amplicons
                        at the variant postion (int) Default = 0
            min_reads:  minimum coverage of variant amplicons (int)
                        Default = 0
            min_freq:   minimum frequence of variant (float 0.0 - 1.0)
                        Default = 0.1
        
        Note that the minimum effective value of min_cover and min_reads
        is bound by AmpliconData.threshold.  Setting min_cover or 
        min_reads to less than AmpliconData.threshold will have no effect.
        """
        def is_softmasked(allele):
            return bool([base for base in allele if base.islower()])
        
        def is_ambiguous(allele):
            return bool([base for base in allele if not base in '-GATC']) #compound events may have gaps
        
        for variant in self.consolidated_variants:
            if variant[1] >= min_reads and \
               variant[2] >= min_cover and \
               variant[3] >= min_freq and \
               not is_softmasked(variant[0].alt_allele) and \
               not is_softmasked(variant[0].ref_allele) and \
               not is_ambiguous(variant[0].alt_allele):
               yield variant
    
    def print_consolidated_vcf(self, min_cover=0, min_reads=0, min_freq=0.1, outfile=sys.stdout):
        """Print variants that match filtering conditions in VCF Format
        
        Arguments:
            min_cover:  minimum total coverage of all amplicons
                        at the variant postion (int) Default = 0
            min_reads:  minimum coverage of variant amplicons (int)
                        Default = 0
            min_freq:   minimum frequence of variant (float 0.0 - 1.0)
                        Default = 0.1
            outfile:    a writeable text format file. Default = sys.stdout
        
        Note that the minimum effective value of min_cover and min_reads
        is bound by AmpliconData.threshold.  Setting min_cover or 
        min_reads to less than AmpliconData.threshold will have no effect
        and the values for these parameters will be reported as threshold
        in the VCF header.
        """
        
        if not self.called_variants:
            self.call_amplicon_variants()
        if not self.consolidated_variants:
            self.consolidate_variants()
        
        vcf_header = [
        "##fileformat=VCFv4.2",
        "##source=AmBiVeRT{version}".format(version=__version__),
        '##FILTER=<ID=depth,Description="more than {threshold} variant supporting reads">'.format(threshold=max(self.threshold,min_reads)),
        '##FILTER=<ID=cover,Description="more than {cover} reads at variant position">'.format(cover=max(self.threshold,min_cover)),
        '##FILTER=<ID=freq,Description="more than {min_freq}% of reads support variant">'.format(min_freq=min_freq*100),
        '##FILTER=<ID=primer,Description="involves primer sequence">',
        '##FILTER=<ID=homopoly5,Description="homopolymer of length > 5 at mutation site">',
        '##FILTER=<ID=homopoly10,Description="homopolymer of length > 10 at mutation site">',
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth excluding soft masked primers at variant site">',
        '##INFO=<ID=AC,Number=1,Type=Integer,Description="Alt allele supporting read count">',
        '##INFO=<ID=AF,Number=1,Type=Float,Description="Alt allele frequency">',
        '##INFO=<ID=ALTAMPS,Number=.,Type=String,Description="Unique identifiers for amplicons supporting alt allele">',
        '##INFO=<ID=REFAMPS,Number=.,Type=String,Description="Unique identifiers for amplicons supporting other alleles">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        ]
        print("\n".join(vcf_header),file=outfile)
        for variant in self.get_filtered_variants(min_cover=min_cover,min_reads=min_reads,min_freq=min_freq):
            soft_filter_name = 'PASS'
            if self.is_homopolymer_at_position(variant[0].chromosome,variant[0].start, minimum=5):
                soft_filter_name = 'homopoly5'
            if self.is_homopolymer_at_position(variant[0].chromosome,variant[0].start, minimum=10):
                soft_filter_name = 'homopoly10'
            print(variant[0].chromosome,
                variant[0].vcf_start,
                ".",
                variant[0].ref_allele,
                variant[0].alt_allele,
                '.',soft_filter_name,
                'DP={DP};AC={AC};AF={AF:.3};ALTAMPS={ALTAMPS};REFAMPS={REFAMPS}'.format(DP=variant[2],
                                                                                    AC=variant[1],
                                                                                    AF=variant[3],
                                                                                    ALTAMPS=",".join(variant[5]),
                                                                                    REFAMPS=",".join(variant[4]),
                                                                                    ),
                sep='\t',file=outfile)
        pass
    
    def get_amplicons_overlapping(self,chrom, pos, length=1):
        """Returns a list of amplicon ids for amplicons that overlap a reference position
        Arguments:
            chrom:  The reference chromosome
            pos:    The reference chromosome base position
            length: Length of the region for overlap
                    Default = 1
        Returns:
            a list of md5 hexdigest keys to AmpliconData.location
            These keys are shared with AmpliconData.data, .aligned and .called_variants
        """
        result = []
        pos = int(pos)
        length = int(length)
        for key in self.location:
            if chrom != self.location[key][0]:
                continue
            start = int(self.location[key][1])
            end = int(self.location[key][2])
            if not ( (pos+length-1 < start) or (pos > end) ):
                result.append(key)
        return result
    
    def get_reference_overlapping(self,chrom, pos, length=1):
        """Returns a list of reference identifiers that overlap a chromosomal position
        Arguments:
            chrom:  The reference chromosome
            pos:    The reference chromosome base position
            length: Length of the region for overlap
                    Default = 1
        Returns:
            a list of keys to AmpliconData.reference_sequences of the
            format (name, chromosome, start, end, strand)
        """
        #self.reference_sequences = {(name, chromosome, start, end, strand): ref_seq,}
        result = []
        pos = int(pos)
        length = int(length)
        for key in self.reference_sequences:
            if chrom != key[1]:
                continue
            start = int(key[2])
            end = int(key[3])
            if not ( (pos+length-1 < start) or (pos > end) ):
                result.append(key)
        return result
    
    def get_amplicons_overlapping_without_softmasking(self,chrom, pos, length=1):
        """Returns a list of amplicon ids for amplicons that overlap a reference position
        where the reference bases are not softmasked (represented by lower case)
        Arguments:
            chrom:  The reference chromosome
            pos:    The reference chromosome base position
            length: Length of the region for overlap
                    Default = 1
        Returns:
            a list of md5 hexdigest keys to AmpliconData.location
            These keys are shared with AmpliconData.data, .aligned and .called_variants
        """
        result = []
        pos = int(pos)
        length = int(length)
        for amplicon_id in self.get_amplicons_overlapping(chrom, pos, length):
            start = self.location[amplicon_id][1]
            end = self.location[amplicon_id][2]
            ungapped_reference = self.aligned[amplicon_id][1].replace('-','')
            if not [base for base in ungapped_reference[pos-start:pos-start+length+1] if base.islower()]:
                result.append(amplicon_id)
        return result
    
    def print_to_fastq(self, key, forwardfile=sys.stdout, reversefile=sys.stdout):
        """Print fastq of amplicons reads - intended for debugging purposes
        Arguments:
            key:            an md5 hexdigest key to AmpliconData.data
            forwardfile:    a writable text format file
                            Default: sys.stdout
            reversefile:    a writable text format file
                            Default: sys.stdout
        """
        for ((f_name, f_seq, f_qual),(r_name, r_seq, r_qual)) in self.data[key]:
            print('@{0}\n{1}\n+\n{2}'.format(f_name,f_seq,f_qual),file=forwardfile)
            print('@{0}\n{1}\n+\n{2}'.format(r_name,r_seq,r_qual),file=reversefile)
        pass
    
    def print_to_sam(self, key, samfile=sys.stdout):
        """Output aligned sequences as SAM with constant quality
        Arguments:
            key:        an md5 hexdigest key to AmpliconData.data
            samfile:    a writable text format file
                        Default: sys.stdout
        """
        # This may be expanded to deal with quality scores in future versions
        # Implemented to allow comparisons on IGV with other programs
        # Output should be run through samtools calmd/fillmd and converted to bam
        if not self.aligned:
            logging.info('Amplicons were not aligned')
            self.align_to_reference()
        aligned_sample_seq, aligned_ref_seq = self.aligned[key]
        cigar, start, length = gapped_alignment_to_cigar(aligned_ref_seq,aligned_sample_seq)
        rname, amplicon_start, amplicon_end, strand = self.location[key]
        pos = amplicon_start + start
        flag = '0'
        seq = aligned_sample_seq.replace('-','')
        qual = 'I'*len(seq)
        for x in range(self.get_amplicon_count(key)):
            qname = '{0}_{1}'.format(key,x+1)
            print(qname,flag,rname,pos,'40',cigar,'*','0',length,seq,qual,sep='\t',file=samfile)
        pass
    
    def printall_to_sam(self, samfile=sys.stdout,
                        sq_header = [
                            '@SQ\tSN:chr1\tLN:249250621',
                            '@SQ\tSN:chr2\tLN:243199373',
                            '@SQ\tSN:chr3\tLN:198022430',
                            '@SQ\tSN:chr4\tLN:191154276',
                            '@SQ\tSN:chr5\tLN:180915260',
                            '@SQ\tSN:chr6\tLN:171115067',
                            '@SQ\tSN:chr7\tLN:159138663',
                            '@SQ\tSN:chr8\tLN:146364022',
                            '@SQ\tSN:chr9\tLN:141213431',
                            '@SQ\tSN:chr10\tLN:135534747',
                            '@SQ\tSN:chr11\tLN:135006516',
                            '@SQ\tSN:chr12\tLN:133851895',
                            '@SQ\tSN:chr13\tLN:115169878',
                            '@SQ\tSN:chr14\tLN:107349540',
                            '@SQ\tSN:chr15\tLN:102531392',
                            '@SQ\tSN:chr16\tLN:90354753',
                            '@SQ\tSN:chr17\tLN:81195210',
                            '@SQ\tSN:chr18\tLN:78077248',
                            '@SQ\tSN:chr19\tLN:59128983',
                            '@SQ\tSN:chr20\tLN:63025520',
                            '@SQ\tSN:chr21\tLN:48129895',
                            '@SQ\tSN:chr22\tLN:51304566',
                            '@SQ\tSN:chrX\tLN:155270560',
                            '@SQ\tSN:chrY\tLN:59373566',
                            '@SQ\tSN:chrM\tLN:16571',
                            '@SQ\tSN:chr3_PATCH_NW_003871060.1\tLN:173151',
                        ]):
        """Output all aligned sequences as SAM with constant quality
        Arguments:
            samfile:    a writable text format file
                        Default: sys.stdout
            sq_header:  a list of SAM @SQ header lines in the format
                        ['@SQ\tSN:chr1\tLN:249250621',]
                        Default: human GRCh37
        """
        print('@HQ\tVN:1.5\tSO:unsorted',file=samfile)
        for sq_line in sq_header:
            print(sq_line,file=samfile)
        print('@PG\tID:AmBiVErT\tVN:{version}'.format(version=__version__),file=samfile)
        for key in sorted(list(self.aligned)):
            self.print_to_sam(key, samfile=samfile)
        pass
    
    def is_homopolymer_at_position(self,chrom, pos, minimum = 5):
        """Test if a homopolymer extends to the right of a sequence position
        Arguments:
            chrom:      Reference chromosome
            pos:        Position on reference chromosome
            minimum:    the minimum number of repeats to be 
                        considered a homopolymer
                        Default = 5
        """
        #get identifiers for reference amplicons that overlap
        reference_ids = self.get_reference_overlapping(chrom, pos, length=1)
        try:
            ref = sorted(reference_ids)[-1] # get right most overlapping reference
        except IndexError: #pragma no cover
            raise IndexError('No reference amplicons overlapping {0} {1}'.format(chrom, pos))
        ref_start = int(ref[2])
        if ref[4] == '-':
            ref_seq = reverse_complement(self.reference_sequences[ref])
        else:
            ref_seq = self.reference_sequences[ref]
        i = pos - ref_start
        firstbase = ref_seq[i].upper()
        while i < len(ref_seq) and ref_seq[i].upper() == firstbase:
            i +=1
        if i - (pos - ref_start) >= minimum:
            return True
        else:
            return False
    
def process_commandline_args(): #pragma no cover
    """Processes command line arguments and returns args object
    Use ambivert --help for details of arguments
    """
    parser = argparse.ArgumentParser(description="""AmBiVErT: A program for binned analysis of amplicon data
        AmBiVErT clusters identical amplicon sequences and thresholds based on read frequency to remove technical errors.
        Due to sequencing errors occuring with a more random distribution than low frequency variants this approach
        reduces the number of amplicon products that must be assigned to target regions & assessed for variant calls.
        AmBiVErT overlaps forward and reverse reads from the same amplicon and preserves local phasing information.
        Typical running time for first use is several hours, which reduces to less than 10 minutes when the
        hash table calculated on a previous run is supplied for analysis of subsequent samples with the same amplicons.
        
        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
        
        This software is provided for research use only.
        """)
    parser.add_argument('-f','--forward',
                        type=str,
                        help='a fastq format file of forward direction amplicon reads. \
                             May be compressed with gzip with the appropriate suffix (.gz)')
    parser.add_argument('-r','--reverse',
                        type=str,
                        help='a fastq format file of reverse direction amplicon reads. \
                             May be compressed with gzip with the appropriate suffix (.gz)')
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
    parser.add_argument('--threshold',
                        type=int,
                        default=20,
                        help='the minimum occurance threshold.  Unique amplicon sequence variants that occur fewer than threshold times are ignored.\
                             For complex variable mixtures of sequences (eg FFPE samples) this parameter should be set to zero.  Default 20')    
    parser.add_argument('--min_cover',
                        type=int,
                        default=0,
                        help='the minimum coverage at a site required to call a variant. \
                            This parameter only has effect if it is > threshold. Default 0')    
    parser.add_argument('--min_reads',
                        type=int,
                        default=0,
                        help='the minimum number of variant containing reads required to call a variant. \
                            This parameter only has effect if it is > threshold. Default 0')    
    parser.add_argument('--min_freq',
                        type=float,
                        default=0.1,
                        help='the minimum proportion of mutated reads. Default 0.1 (ten percent)')    
    parser.add_argument('--overlap',
                        type=int,
                        default=20,
                        help='The minimum overlap required between forward and reverse sequences to merge. Default 20bp')    
    parser.add_argument('--hashtable',
                        type=argparse.FileType('rb'),
                        help='Filename for a precomputed hash table of exact matches of amplicons to references.  Generate with --savehashtable')    
    parser.add_argument('--savehashtable',
                        type=argparse.FileType('wb'),
                        help='Output a precomputed hash table that matches amplicons exactly to references.  Used to speed up matching with --hashtable')    
    parser.add_argument('--threads',
                        type=int,
                        default=1,
                        help='Number of concurrent threads to use for mapping amplicons to references.  Default 1')    
    parser.add_argument('--alignments',
                        type=str,
                        default='',
                        help='Print a formatted text version of variant containing alignments to a file. "-" will print to stderr')    
    parser.add_argument('--sam',
                        type=argparse.FileType('wt'),
                        help='Output a sam format alignment of the aligned amplicons.  Uses constant quality scores,\
                             ignores amplicon strand and does not produce md tag calls')    
    parser.add_argument('--fastqamplicon',
                        type=str,
                        default='',
                        help='Retrieve fastq output for specified amplicon (specified by MD5 identifier).\
                          Analysis is not run, only importing of data.  Must specify --fastqfilename')    
    parser.add_argument('--fastqfilename',
                        type=str,
                        default='',
                        help='Base file name for forward and reverse fastqfiles for use with --fastqamplicon\
                           Defaults to amplicon MD5 identifier')    
    parser.add_argument('--prefix',
                        type=str,
                        default='',
                        help='Shorthand specification of --forward, --reverse, --countfile and --outfile. \
                                --forward = <prefix>_R1.fastq.gz \
                                --reverse = <prefix>_R2.fastq.gz \
                                --outfile = <prefix>.vcf \
                                --countfile = <prefix>.counts')
    parser.add_argument('--version',
                        action='store_true',
                        help='print version information and exit')
    parser.add_argument('--test',
                        action='store_true',
                        help='run tests before starting, and quit if they don\'t pass')
    args = parser.parse_args()
    if args.test:
        import unittest
        import ambivert.tests.test_ambivert
        logging.info('Running test suite to confirm correct functionality')
        logging.disable(logging.CRITICAL)
        loader = unittest.TestLoader()
        tests = loader.loadTestsFromModule(ambivert.tests.test_ambivert)
        ambivert_test_output=io.StringIO()
        runner = unittest.TextTestRunner(ambivert_test_output)
        result = runner.run(tests)
        logging.disable(logging.NOTSET)
        if not result.wasSuccessful():
            raise RuntimeError('Error in executing tests - this is most likely an installation issue\n{0}'.format(ambivert_test_output.getvalue()))
        else:
            logging.info('Tests completed successfully')
    if args.version:
        print(__version__)
        sys.exit()
    if args.prefix:
        args.countfile = open(args.prefix + '.counts','w')
        args.output = open(args.prefix + '.vcf','w')
        args.forward = open(args.prefix + '_R1.fastq.gz','rb')
        args.reverse = open(args.prefix + '_R2.fastq.gz','rb')
    elif not (args.forward and args.reverse):
        parser.print_help()
        raise RuntimeError("ERROR: You must specify --prefix or both --forward and --reverse")
    if args.alignments:
        if args.alignments == '-':
            args.alignments = sys.stderr
        else:
            args.alignments = open(args.alignments,'wt')
    if not args.manifest and not args.fasta and not args.fastqamplicon:
        parser.print_help()
        raise RuntimeError("ERROR: You must specify at least one of --manifest and --fasta")
    if args.fastqamplicon and not args.fastqfilename:
        args.fastqfilename = args.fastqamplicon
    return args

def process_amplicon_data(forward_file, reverse_file,
                          manifest=None, fasta=None,
                          threshold=50, overlap=20, 
                          savehashtable=None, hashtable=None,
                          threads=1):
    """Read amplicon data from fastq files and return
    an AmpliconData object with overlapped and aligned
    to reference amplicons
    
    Arguments:
        forward_file : a filename string or file object containing fastq
                       format forward amplicon reads.  Can be gzipped.
        reverse_file : a filename string or file object containing fastq
                       format reverse amplicon reads.  Can be gzipped.
        manifest     : an Illumina TruSeq amplicon manifest file
        fasta        : a file object containg fasta entries of reference
                       sequences.  Reference sequences should include the 
                       entire amplicon sequence.  Regions that should not
                       be called (eg primers) should be soft masked in
                       lower case.  All other bases should be upper case.
                       Additional annotation details are required on the 
                       name line and should be tab separated in the format:
                       >amplicon_name chromosome start end strand(+/-)
                       
                       One of manifest or fasta must be specified.
                       If both manifest and fasta are specified all
                       reference sequences from both files will be used.
        
        threshold    : the minimum number of occurances of an amplicon
                       with a specific sequence required to process.
                       This value limits the calling parameters min_reads
                       and min_depth.  Smaller values will increase time
                       required to process data by limiting complexity
                       reduction.
                       Default: 50.
        overlap      : the minimum number of bases required for forward and
                       reverse reads to be considered overlapping.
        savehashtable: a writable binary format file for saving the mapping
                       to reference lookup table
        hashtable    : a readable binary format file containing a previously
                       saved mapping to reference lookup table.
        threads      : Number of concurrent threads to use for mapping 
                       amplicons to references. Default: 1
    
    Returns:
        An ambivert.ambivert.AmpliconData object.
    """
    amplicons = AmpliconData()
    amplicons.process_twofile_readpairs(open_potentially_gzipped(forward_file),
                                        open_potentially_gzipped(reverse_file))
    amplicons.merge_overlaps(minimum_overlap=overlap, threshold=threshold)
    if manifest:
        amplicons.add_references_from_manifest(manifest)
    if fasta: #pragma no cover
        amplicons.add_references_from_fasta(fasta) 
    if hashtable: #pragma no cover
        amplicons.load_hash_table(hashtable)
    amplicons.match_to_reference(show_progress = True, threads = threads)
    amplicons.align_to_reference()
    if savehashtable: #pragma no cover
        amplicons.save_hash_table(savehashtable)
    return amplicons

def call_variants_per_amplicon(amplicons, args): #pragma no cover
    """Depreciated.  Use AmpliconData.print_consolidated_vcf or 
    AmpliconData.call_amplicon_variants() and .consolidate_variants()
    """
    #deprecated
    warnings.warn('call_variants_per_amplicon is deprecated',DeprecationWarning)
    print(make_vcf_header(args.threshold),file=args.output)
    for key in amplicons.potential_variants:
            aligned_ref_seq, aligned_sample_seq = amplicons.aligned[key]
            name, chromosome, amplicon_position, end, strand = amplicons.reference[key]
            try:
                call_variants_to_vcf(aligned_sample_seq, aligned_ref_seq, chromosome, int(amplicon_position), outfile=args.output)
            except:
                print('WARNING: ',name,' FAILED to CALL',file=sys.stderr)
                print(aligned_ref_seq,file=sys.stderr)
                print(aligned_sample_seq,file=sys.stderr)
                print(file=sys.stderr)
                raise
    pass

def main(): #pragma no cover
    args = process_commandline_args()
    
    if args.fastqamplicon:
        amplicons = AmpliconData()
        amplicons.process_twofile_readpairs(open_potentially_gzipped(args.forward),
                                            open_potentially_gzipped(args.reverse))
        amplicons.merge_overlaps(minimum_overlap=args.overlap, threshold=args.threshold)
        amplicons.print_to_fastq(args.fastqamplicon,
                                forwardfile=open(args.fastqfilename+'_R1.fastq','wt'),
                                reversefile=open(args.fastqfilename+'_R2.fastq','wt'))
        sys.exit()
        
        
    amplicons = process_amplicon_data(args.forward,args.reverse,
                                      args.manifest,args.fasta,
                                      args.threshold,args.overlap,
                                      args.savehashtable,args.hashtable,
                                      args.threads)
    
    if args.countfile:
        with args.countfile as outfile:
            amplicon_counts = amplicons.get_amplicon_counts()
            for key in sorted(amplicon_counts.keys()):
                outfile.write('{0}\t{1}\n'.format(key,amplicon_counts[key]))
    
    if args.alignments:
        amplicons.print_variants_as_alignments(outfile=args.alignments)
    
    #These are called in amplicons.print_consolidated_vcf
    #amplicons.call_amplicon_variants()
    #amplicons.consolidate_variants()
    
    amplicons.print_consolidated_vcf(min_cover=args.min_cover, min_reads=args.min_reads, min_freq=args.min_freq, outfile=args.output)
    
    if args.sam:
        amplicons.printall_to_sam(samfile=args.sam)
    pass

if __name__ == '__main__': #pragma no cover
    main()
