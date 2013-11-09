#!/usr/bin/env python2.7
# encoding: utf-8
"""
ambivert.py

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
complementing a more traditional variant calling pipeline.

All lines of code [will be/are] covered by unit tests unless marked with #pragma no cover

Created by Matthew Wakefield and Graham Taylor.
Copyright (c) 2013  Matthew Wakefield and The University of Melbourne. All rights reserved.
"""
from __future__ import print_function, division
import sys, os
import itertools, difflib, argparse
import hashlib, cPickle
from collections import defaultdict
#from cogent.align.algorithm import nw_align, sw_align
from sequence_utilities import parse_fastq, parse_fasta, reverse_complement, flatten_paired_alignment, format_alignment, open_potentially_gzipped
from truseq_manifest import parse_truseq_manifest, make_sequences
from call_mutations import call_mutations, call_mutations_to_vcf, make_vcf_header
import plumb.bob
    

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
    alignment =  plumb.bob.local_align(bytes(seq1), len(seq1),
                                bytes(seq2.upper()), len(seq2),
                                plumb.bob.DNA_MAP[0],
                                plumb.bob.DNA_MAP[1], 
                                plumb.bob.DNA_SCORE,
                                -7, -1 #gap open, gap extend
                                )
    start_seq1 = alignment.contents.align_frag.contents.sa_start
    start_seq2 = alignment.contents.align_frag.contents.sb_start
    frag = alignment[0].align_frag
    align_seq1 = ''
    align_seq2 = ''
    while frag:
        frag = frag[0]
        if frag.type == plumb.bob.MATCH:
            f1 = seq1[frag.sa_start:frag.sa_start + frag.hsp_len]
            f2 = seq2[frag.sb_start:frag.sb_start + frag.hsp_len]
            align_seq1 += f1
            align_seq2 += f2
        elif frag.type == plumb.bob.A_GAP:
            align_seq1 += '-' * frag.hsp_len
            align_seq2 += seq2[frag.sb_start:frag.sb_start + frag.hsp_len]
        elif frag.type == plumb.bob.B_GAP:
            align_seq1 += seq1[frag.sa_start:frag.sa_start + frag.hsp_len]
            align_seq2 += '-' * frag.hsp_len
        frag = frag.next
    plumb.bob.alignment_free(alignment)
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
        self.location = {} # dictionary of (chromosome, start, end, strand) with same key as self.data
        self.reference = {} # dictionary of reference sequence keys with same key as self.data
        self.reference_sequences = {} # dictionary of reference sequences. Key is (name, chromosome, start, end, strand)
        self.potential_variants = [] # list of keys (to self.data) for which putative variants were identified
        self.called_mutations = {} # dictionary of lists of called mutations with same key as self.data
        self.consolidated_mutations = []
        self.trim5 = trim5
        self.trim3 = None if trim3==None else -1*abs(trim3)
        self.threshold = 0
        self.readpairs = 0
        pass
    
    def __str__(self):
        return self.data.__str__()
    
    def add_reads(self, f_name, f_seq, f_qual,r_name, r_seq, r_qual):
        amplicon_key = hashlib.md5(f_seq[self.trim5:self.trim3]+r_seq[self.trim5:self.trim3]).hexdigest()
        self.data[amplicon_key].append(((f_name, f_seq, f_qual),(r_name, r_seq, r_qual)))
        self.readpairs += 1
        pass
    
    def process_twofile_readpairs(self, forward_file, reverse_file, parser=parse_fastq):
        print('Reading data...', file=logfile)
        if hasattr(forward_file,'name'):
            print('    Forward read file: ', forward_file.name, file=logfile)
        if hasattr(reverse_file,'name'):
            print('    Reverse read file: ', reverse_file.name, file=logfile)
        for (f_name, f_seq, f_qual),(r_name, r_seq, r_qual) in itertools.izip(parser(forward_file), parser(reverse_file)):
            assert f_name.split()[0] == r_name.split()[0]
            self.add_reads(f_name, f_seq, f_qual,r_name, r_seq, r_qual)
        print('Read',self.readpairs,'read pairs', file=logfile)
        pass
    
    #def process_interleaved_readpairs(parser=parse_fastq, filename):
    #    pass
    
    def get_above_threshold(self, threshold=1):
        return [x for x in self.data if len(self.data[x]) > threshold]
    
    def merge_overlaps(self, threshold=0, minimum_overlap=10):
        self.threshold = threshold #record in object data
        total = len(self.get_above_threshold(threshold))
        completed = 0
        print('Merging overlaps for {0} unique pairs'.format(total),file=logfile)
        for key in self.get_above_threshold(threshold):
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
        self.reference_sequences = {tuple(x[0].split()):x[1] for x in parse_fasta(fastafile)}
        pass
    
    def add_references_from_manifest(self, manifestfile):
        for name,sequence in make_sequences(*parse_truseq_manifest(manifestfile), with_probes=True, softmask_probes=True, all_plus=False):
            self.reference_sequences[tuple(name.split())] = sequence
        pass
    
    def match_to_reference(self, min_score = 0.1, trim_primers=0):
        def match_by_edit(merged_key,ref_keys):
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
            return best_hit,score
        
        def match_by_smith_waterman(merged_key,ref_keys):
            best_score = 0
            best_hit = ''
            for ref_key in ref_keys:
                if trim_primers:
                    seq1 = self.merged[merged_key][trim_primers:-trim_primers]
                else:
                    seq1 = self.merged[merged_key]
                alignment =  plumb.bob.local_align(seq1, len(seq1),
                                self.reference_sequences[ref_key].upper(), len(self.reference_sequences[ref_key]),
                                plumb.bob.DNA_MAP[0],
                                plumb.bob.DNA_MAP[1], 
                                plumb.bob.DNA_SCORE,
                                -7, -1 #gap open, gap extend
                                )
                score = alignment[0].score
                plumb.bob.alignment_free(alignment)
                if score > best_score:
                    best_score = score
                    best_hit = ref_key
            return best_hit,score

                
        print('Matching read bins to references',file=logfile)
        for merged_key in self.merged:
            if merged_key in self.reference:
                print(',',end='',file=logfile)
                continue
            best_hit,best_score = match_by_smith_waterman(merged_key,self.reference_sequences)
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
            if self.reference[merged_key][-1] == '-': #last element of self.reference_sequence key indicates minus strand probe
                #both the reference sequence and the merged sequence are on the minus strand
                query_seq = reverse_complement(self.merged[merged_key])
                ref_seq = reverse_complement(self.reference_sequences[self.reference[merged_key]])
            else:
                query_seq = self.merged[merged_key]
                ref_seq = self.reference_sequences[self.reference[merged_key]]
            aligned_sample_seq,aligned_ref_seq,sample_start,ref_start = smith_waterman(query_seq,ref_seq)
            if aligned_ref_seq.upper() != aligned_sample_seq:
                self.potential_variants.append(merged_key)
                #print(aligned_ref_seq, file=logfile)
                #print(aligned_sample_seq, file=logfile)
            self.aligned[merged_key] = (aligned_sample_seq, aligned_ref_seq) #plus strand
            self.location[merged_key] = (self.reference[merged_key][1], int(self.reference[merged_key][2]), int(self.reference[merged_key][3]), ref_start)
        pass
    
    def get_amplicon_count(self, key):
        return len(self.data[key])
    
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
            hash_dictionary = {merged_key:self.reference[merged_key] for merged_key in self.reference if merged_key not in self.potential_variants} 
            cPickle.dump((reference_sha224,hash_dictionary),outfile)
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
    
    def print_variants_as_alignments(self, outfile=sys.stdout):
        for key in self.potential_variants:
            aligned_sample_seq, aligned_ref_seq = self.aligned[key]
            print(self.reference[key],file=outfile)
            print(aligned_ref_seq,file=outfile)
            matches = ''
            for a,b in itertools.izip(aligned_ref_seq, aligned_sample_seq):
                if a == b or a in 'abcdghkmnrstuvwy':
                    matches += '.'
                else:
                    matches += b
            print(matches,file=outfile)
        pass
    
    def call_amplicon_mutations(self):
        for key in self.potential_variants:
            aligned_ref_seq, aligned_sample_seq = self.aligned[key]
            name, chromosome, amplicon_position, end, strand = self.reference[key]
            self.called_mutations[key] = call_mutations(aligned_sample_seq, aligned_ref_seq, chromosome, int(amplicon_position))
        pass
    
    def get_variant_positions(self):
        variant_positions = []
        for amplicon_id in self.called_mutations:
            for variant in self.called_mutations[amplicon_id]:
                variant_positions.append((variant.chromosome, variant.start, variant.end))
        return sorted(list(set(variant_positions)))
    
    def consolidate_mutations(self, exclude_softmasked_coverage=True):
        positions = self.get_variant_positions()
        for (chrom, start, end) in positions:
            if exclude_softmasked_coverage:
                amplicon_ids = self.get_amplicons_overlapping_without_softmasking(chrom,start,end-start)
            else:
                amplicon_ids = self.get_amplicons_overlapping(chrom,start,end-start)
            ref = [] #list of amplicon keys
            alt = {} #dictionary of lists of amplicons keyed by alternative alleles
            for amplicon_id in amplicon_ids:
                overlapping_mutation = False
                if amplicon_id in self.called_mutations:
                    for mutation in self.called_mutations[amplicon_id]:
                        if not ((mutation.end < start) or (mutation.start > end)):
                            overlapping_mutation = True
                            if mutation in alt:
                                alt[mutation].append(amplicon_id)
                            else:
                                alt[mutation] = [amplicon_id,]
                if not overlapping_mutation:
                    ref.append(amplicon_id)
            
            total_depth = sum([self.get_amplicon_count(amplicon_id) for amplicon_id in amplicon_ids])
            
            #indels are not described at the same coordinate as snps so we deal with these first
            #we dont do 'correct' VCF with all alt alleles on one line, we repeat positions instead.
            #in rare cases of two mutations overlapping a deletion we just duplicate an identical deletion entry
            #exact duplicate records are then cleaned up at end of this function
            #this is structured to allow compound calling in a subsequent version
            indels = [variant for variant in alt.keys() if len(variant.alt_allele) > 1 or len(variant.ref_allele) > 1 ]
            for indel in indels:
                variant_depth = sum([self.get_amplicon_count(amplicon_id) for amplicon_id in alt[indel]])
                #print(indel,variant_depth, total_depth, variant_depth/total_depth, ref, alt[indel])
                self.consolidated_mutations.append((indel,variant_depth, total_depth, variant_depth/total_depth, tuple(ref), tuple(alt[indel])))
            
            #for now we do one variant per line for snps
            #this should change to at least an optional correct VCF formatting of all alleles on one line
            snvs = [key for key in alt.keys() if len(key.alt_allele) == 1 and len(key.ref_allele) == 1 ]
            for snv in snvs:
                variant_depth = sum([self.get_amplicon_count(key) for key in alt[snv]])
                #print(snv,variant_depth, total_depth, variant_depth/total_depth, ref, alt[snv])
                self.consolidated_mutations.append((snv,variant_depth, total_depth, variant_depth/total_depth, tuple(ref), tuple(alt[snv])))
        #logic above does not preclude calling the same mutation more than once
        #so we remove identical records
        self.consolidated_mutations = sorted(list(set(self.consolidated_mutations)))
        pass
        
    def get_filtered_mutations(self, min_cover=0, min_reads=0, min_freq=0.1):
        for variant in self.consolidated_mutations:
            if variant[1] >= min_reads and \
               variant[2] >= min_cover and \
               variant[3] >= min_freq and \
               not [base for base in variant[0].alt_allele if base.islower()] and \
               not [base for base in variant[0].ref_allele if base.islower()]:
               yield variant
        
    def print_consolidated_vcf(self, min_cover=0, min_reads=0, min_freq=0.1, outfile=sys.stdout):
        if not self.called_mutations:
            self.call_amplicon_mutations()
        if not self.consolidated_mutations:
            self.consolidate_mutations()
        
        vcf_header = [
        "##fileformat=VCF4.1",
        "##source=AmBiVeRT0.1.0",
        '##FILTER=<ID=depth,Description="more than {threshold} variant supporting reads">'.format(threshold=max(self.threshold,min_reads)),
        '##FILTER=<ID=cover,Description="more than {cover} reads at variant position">'.format(cover=max(self.threshold,min_cover)),
        '##FILTER=<ID=freq,Description="more than {min_freq}% of reads support variant">'.format(min_freq=min_freq*100),
        '##FILTER=<ID=primer,Description="involves primer sequence">',
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth excluding soft masked primers at variant site">',
        '##INFO=<ID=AC,Number=1,Type=Integer,Description="Alt allele supporting read count">',
        '##INFO=<ID=AF,Number=1,Type=Float,Description="Alt allele frequency">',
        '##INFO=<ID=ALTAMPS,Number=.,Type=String,Description="Unique identifiers for amplicons supporting alt allele">',
        '##INFO=<ID=REFAMPS,Number=.,Type=String,Description="Unique identifiers for amplicons supporting other alleles">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        ]
        print("\n".join(vcf_header),file=outfile)
        for variant in self.get_filtered_mutations(min_cover=min_cover,min_reads=min_reads,min_freq=min_freq):
            print(variant[0].chromosome,
                variant[0].vcf_start,
                ".",
                variant[0].ref_allele,
                variant[0].alt_allele,
                '.','PASS',
                'DP={DP};AC={AC};AF={AF:.3};ALTAMPS={ALTAMPS};REFAMPS={REFAMPS}'.format(DP=variant[2],
                                                                                    AC=variant[1],
                                                                                    AF=variant[3],
                                                                                    ALTAMPS=",".join(variant[5]),
                                                                                    REFAMPS=",".join(variant[4]),
                                                                                    ),
                sep='\t',file=outfile)
        pass
            

    def get_amplicons_overlapping(self,chrom, pos, length):
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

    def get_amplicons_overlapping_without_softmasking(self,chrom, pos, length):
        result = []
        pos = int(pos)
        length = int(length)
        for amplicon_id in self.get_amplicons_overlapping(chrom, pos, length):
            start = self.location[amplicon_id][1]
            end = self.location[amplicon_id][2]
            query_start_in_amplicon = pos - start
            length_in_ref = 0
            position_in_ref = query_start_in_amplicon
            while length_in_ref < length and position_in_ref < len(self.aligned[amplicon_id][0]):
                if self.aligned[amplicon_id][0][position_in_ref] != '-':
                    length_in_ref += 1
                position_in_ref +=1
            query_end_in_amplicon = position_in_ref
            if not [base for base in self.aligned[amplicon_id][0][query_start_in_amplicon:query_end_in_amplicon] if base.islower()]:
                result.append(amplicon_id)
        return result
    
    def print_to_fastq(self, key, forwardfile=sys.stdout, reversefile=sys.stdout):
        """Print fastq of amplicons reads - intended for debugging purposes"""
        for ((f_name, f_seq, f_qual),(r_name, r_seq, r_qual)) in self.data[key]:
            print('@'+f_name,file=forwardfile)
            print(f_seq,file=forwardfile)
            print('+',file=forwardfile)
            print(f_qual,file=forwardfile)
            print('@'+r_name,file=reversefile)
            print(r_seq,file=reversefile)
            print('+',file=reversefile)
            print(r_qual,file=reversefile)
        pass
    
def process_commandline_args():
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
    #parser.add_argument('--mpileup',
    #                    type=argparse.FileType('w'),
    #                    help='output mpileup files with filenames of the format "<name_provided_as_argument>_<amplicon_name>.mpileup"')    
    parser.add_argument('--threshold',
                        type=int,
                        default=20,
                        help='the minimum occurance threshold.  Unique amplicon sequence variants that occur fewer than threshold times are ignored.')    
    parser.add_argument('--min_cover',
                        type=int,
                        default=0,
                        help='the minimum coverage at a site required to call a mutation. \
                            This parameter only has effect if it is > threshold. Default 0')    
    parser.add_argument('--min_reads',
                        type=int,
                        default=0,
                        help='the minimum number of mutation containing reads required to call a mutation. \
                            This parameter only has effect if it is > threshold. Default 0')    
    parser.add_argument('--min_freq',
                        type=float,
                        default=0.1,
                        help='the minimum proportion of mutated reads. Default 0.1 (ten percent)')    
    parser.add_argument('--overlap',
                        type=int,
                        default=20,
                        help='The minimum overlap required between forward and reverse sequences to merge')    
    parser.add_argument('--primer',
                        type=int,
                        default=0,
                        help='The size of the smallest primer.  This number of bases is trimmed from the end of the merged sequences \
                             to reduce the possibility that small amplicons will fail to match due to primer mismatch')    
    parser.add_argument('--hashtable',
                        type=argparse.FileType('U'),
                        help='Filename for a precomputed hash table of exact matches of amplicons to references.  Generate with --savehashtable')    
    parser.add_argument('--savehashtable',
                        type=argparse.FileType('w'),
                        help='Output a precomputed hash table that matches amplicons exactly to references.  Used to speed up matching with --hashtable')    
    parser.add_argument('--prefix',
                        type=str,
                        default='',
                        help='Shorthand specification of --forward, --reverse, --countfile and --outfile \
                                --forward = <prefix>_R1.fastq.gz \
                                --reverse = <prefix>_R2.fastq.gz \
                                --outfile = <prefix>.vcf \
                                --countfile = <prefix>.counts')
    args = parser.parse_args()
    if args.prefix:
        args.countfile = open(args.prefix + '.counts','w')
        args.output = open(args.prefix + '.vcf','w')
        args.forward = open(args.prefix + '_R1.fastq.gz','r')
        args.reverse = open(args.prefix + '_R2.fastq.gz','r')
    return args

def process_amplicon_data(forward_file, reverse_file,
                          manifest=None, fasta=None,
                          threshold=50, overlap=20, primer=15, 
                          savehashtable=None, hashtable=None,
                          ):
    amplicons = AmpliconData()
    amplicons.process_twofile_readpairs(open_potentially_gzipped(forward_file),
                                        open_potentially_gzipped(reverse_file))
    amplicons.merge_overlaps(minimum_overlap=20, threshold=threshold)
    if manifest:
        amplicons.add_references_from_manifest(manifest)
    if fasta:
        amplicons.add_references_from_fasta(fasta)
    if hashtable:
        amplicons.load_hash_table(hashtable)
    amplicons.match_to_reference()
    amplicons.align_to_reference()
    if savehashtable:
        amplicons.save_hash_table(savehashtable)
    return amplicons

def call_mutations_per_amplicon(amplicons, args):
    print(make_vcf_header(args.threshold),file=args.output)
    for key in amplicons.potential_variants:
            aligned_ref_seq, aligned_sample_seq = amplicons.aligned[key]
            name, chromosome, amplicon_position, end, strand = amplicons.reference[key]
            try:
                call_mutations_to_vcf(aligned_sample_seq, aligned_ref_seq, chromosome, int(amplicon_position), outfile=args.output)
            except:
                print('WARNING: ',name,' FAILED to CALL',file=sys.stderr)
                print(aligned_ref_seq,file=sys.stderr)
                print(aligned_sample_seq,file=sys.stderr)
                print(file=sys.stderr)
                raise
    pass

    
    
def main():
    args = process_commandline_args()
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
    
    #amplicons.print_variants_as_alignments(outfile=sys.stderr)
    
    #amplicons.call_amplicon_mutations()
    
    #print(amplicons.get_variant_positions())
    
    #amplicons.consolidate_mutations()
    
    #call_mutations_per_amplicon(amplicons,args)
    
    #TODO parse args to this function
    amplicons.print_consolidated_vcf(min_cover=args.min_cover, min_reads=args.min_reads, min_freq=args.min_freq, outfile=sys.stdout)
    
    pass

if __name__ == '__main__':
    main()
