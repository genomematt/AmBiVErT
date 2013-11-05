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
from __future__ import print_function
import sys, os
import itertools, difflib, argparse
import hashlib, cPickle
from collections import defaultdict
#from cogent.align.algorithm import nw_align, sw_align
from sequence_utilities import parse_fastq, parse_fasta, reverse_complement, flatten_paired_alignment, format_alignment
from truseq_manifest import parse_truseq_manifest, make_sequences
from call_mutations import call_mutations, make_vcf_header
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
    alignment =  plumb.bob.local_align(seq1, len(seq1),
                                seq2.upper(), len(seq2),
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
        print('Reading data...', file=logfile)
        if hasattr(forward_file,'name'):
            print('    Forward read file: ', forward_file.name, file=logfile)
        if hasattr(reverse_file,'name'):
            print('    Reverse read file: ', reverse_file.name, file=logfile)
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

    
    def get_amplicons_overlapping(self,chrom, pos, length):
        result = []
        for key in self.location:
            if chrom != self.location[key][0]:
                continue
            start = self.location[key][1]
            end = self.location[key][2]
            if not ( (int(pos)+int(length)-1 < int(start)) or (int(pos) > int(end)) ):
                result.append(key)
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
                        default=0,
                        help='The size of the smallest primer.  This number of bases is trimmed from the end of the merged sequences \
                             to reduce the possibility that small amplicons will fail to match due to primer mismatch')    
    parser.add_argument('--hashtable',
                        type=argparse.FileType('U'),
                        help='Filename for a precomputed hash table of exact matches of amplicons to references.  Generate with --savehashtable')    
    parser.add_argument('--savehashtable',
                        type=argparse.FileType('w'),
                        help='Output a precomputed hash table that matches amplicons exactly to references.  Used to speed up matching with --hashtable')    
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
    amplicons.align_to_reference()
    if savehashtable:
        amplicons.save_hash_table(savehashtable)
    return amplicons

def main():
    #testargs = ['--forward','/Users/wakefield/Software/ABBA/130503_M00267_0038_L001_MS0302_12015B_E0023_TruSeq_R1.fastq',
    #            '--reverse','/Users/wakefield/Software/ABBA/130503_M00267_0038_L001_MS0302_12015B_E0023_TruSeq_R2.fastq',
    #            '--manifest','/Users/wakefield/Software/ABBA/TruSeq_CAT_Manifest_Allocate-CAT.txt',
    #            #'--hashtable','/Users/wakefield/Software/ABBA/cached_read_to_reference_table',
    #            '--savehashtable','/Users/wakefield/Software/ABBA/new_cached_read_to_reference_table',
    #            '--threshold','20',
    #            '--overlap','20',
    #            '--primer','15',
    #            '--countfile','/Users/wakefield/Software/ABBA/130503_M00267_0038_L001_MS0302_12015B_E0023_TruSeq_AmpliconCounts',
    #            ]
    #args = process_commandline_args(testargs)
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
    
    amplicons.print_variants_as_alignments(outfile=sys.stderr)
    
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
    
    ##debug code
    #keys = amplicons.get_amplicons_overlapping('chr17',41243125,41247176)
    #print(keys)
    #forwardfile = open('testdata_R1.fastq','w')
    #reversefile = open('testdata_R2.fastq','w')
    #for key in keys:
    #    if amplicons.reference[key] in [('BRCA1_Exon9_UserDefined_(9825051)_7473609_chr17_41243125_41243349', 'chr17', '41243125', '41243349', '+'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473610_chr17_41243267_41243491', 'chr17', '41243267', '41243491', '-'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473611_chr17_41243405_41243630', 'chr17', '41243405', '41243630', '+'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473612_chr17_41243559_41243783', 'chr17', '41243559', '41243783', '-'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473613_chr17_41243701_41243927', 'chr17', '41243701', '41243927', '+'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473614_chr17_41243841_41244065', 'chr17', '41243841', '41244065', '-'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473615_chr17_41243981_41244206', 'chr17', '41243981', '41244206', '+'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473616_chr17_41244123_41244347', 'chr17', '41244123', '41244347', '-'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473617_chr17_41244261_41244486', 'chr17', '41244261', '41244486', '+'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473618_chr17_41244399_41244625', 'chr17', '41244399', '41244625', '-'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473619_chr17_41244539_41244764', 'chr17', '41244539', '41244764', '+'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473620_chr17_41244679_41244909', 'chr17', '41244679', '41244909', '-'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473621_chr17_41244855_41245081', 'chr17', '41244855', '41245081', '+'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473622_chr17_41245027_41245253', 'chr17', '41245027', '41245253', '-'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473623_chr17_41245195_41245420', 'chr17', '41245195', '41245420', '+'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473624_chr17_41245363_41245588', 'chr17', '41245363', '41245588', '-'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473625_chr17_41245533_41245757', 'chr17', '41245533', '41245757', '+'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473626_chr17_41245703_41245928', 'chr17', '41245703', '41245928', '-'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473627_chr17_41245865_41246105', 'chr17', '41245865', '41246105', '+'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473628_chr17_41246051_41246277', 'chr17', '41246051', '41246277', '-'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473629_chr17_41246219_41246443', 'chr17', '41246219', '41246443', '+'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473630_chr17_41246387_41246623', 'chr17', '41246387', '41246623', '-'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473631_chr17_41246571_41246801', 'chr17', '41246571', '41246801', '+'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473632_chr17_41246745_41246989', 'chr17', '41246745', '41246989', '-'),
    #    ('BRCA1_Exon9_UserDefined_(9825051)_7473633_chr17_41246933_41247176', 'chr17', '41246933', '41247176', '+'),
    #    ]:
    #        amplicons.print_to_fastq(key, forwardfile, reversefile)
    
    

    pass

if __name__ == '__main__':
    main()
