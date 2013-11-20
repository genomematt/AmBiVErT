#!/usr/bin/env python2.7
# encoding: utf-8
"""
ambivert.py

Created by Matthew Wakefield and Graham Taylor.
Copyright (c) 2013  Matthew Wakefield and The University of Melbourne. All rights reserved.
"""
from __future__ import print_function
import sys, os
from collections import namedtuple
#import itertools, difflib, argparse
#import hashlib, cPickle
#from collections import defaultdict
##from cogent.align.algorithm import nw_align, sw_align
#from sequence_utilities import parse_fastq, parse_fasta, reverse_complement, flatten_paired_alignment, format_alignment
#from truseq_manifest import parse_truseq_manifest, make_sequences
#import plumb.bob
#    

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2013,  Matthew Wakefield and The University of Melbourne"
__credits__ = ["Matthew Wakefield","Graham Taylor"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Matthew Wakefield"
__email__ = "matthew.wakefield@unimelb.edu.au"
__status__ = "Development"

logfile = sys.stderr

def caller(mutant, reference, gap='-',softmask = True):
    mismatch = ''
    deletion = ''
    insertion = ''
    onebased_pos_in_ungapped_ref = 0
    softmasked = True
    
    for site in range(len(mutant)):
        #ignore softmasked sites or gaps in softmasked sequence
        if softmask and softmasked and reference[site] == reference[site].lower():
            #in upstream softmasking
            onebased_pos_in_ungapped_ref += 1
            continue
        elif softmask and softmasked and reference[site] == gap:
            gaplength = 0
            while reference[site+gaplength] == gap:
                gaplength +=1
            if reference[site+gaplength] == reference[site+gaplength].lower():
                #gap is within softmasking
                continue
            #gap abutts softmasking but is in unmasked sequence
        elif softmask and (not softmasked) and not reference[site] == gap and reference[site] == reference[site].lower():
            #we have reached downstream soft masking
            break
        else:
            #we have reached end of softmasking
            softmasked = False

        if mutant[site].upper() == reference[site].upper() and not reference[site] == gap and not mutant[site] == gap:
            onebased_pos_in_ungapped_ref += 1
            if not mismatch and not deletion and not insertion:
                continue
            elif mismatch:
                yield 'X', site-len(mismatch), onebased_pos_in_ungapped_ref-len(mismatch), mismatch
                mismatch = ''
            elif deletion:
                yield 'D', site-len(deletion), onebased_pos_in_ungapped_ref-len(deletion), deletion
                deletion = ''
            elif insertion:
                yield 'I', site-len(insertion), onebased_pos_in_ungapped_ref, insertion
                insertion = ''
        elif reference[site] == gap:
            insertion += mutant[site]
            if mismatch:
                yield 'X', site-len(mismatch), onebased_pos_in_ungapped_ref-len(mismatch), mismatch
                mismatch = ''
            elif deletion:
                yield 'D', site-len(deletion), onebased_pos_in_ungapped_ref-len(deletion), deletion
                deletion = ''
        elif mutant[site] == gap:
            deletion += reference[site]
            onebased_pos_in_ungapped_ref += 1
            if mismatch:
                yield 'X', site-len(mismatch), onebased_pos_in_ungapped_ref-len(mismatch), mismatch
                mismatch = ''
            elif insertion:
                yield 'I', site-len(insertion), onebased_pos_in_ungapped_ref, insertion
                insertion = ''
        else:
            onebased_pos_in_ungapped_ref += 1
            if mismatch:
                yield 'X', site-len(mismatch), onebased_pos_in_ungapped_ref-len(mismatch), mismatch
                mismatch = ''
            mismatch += mutant[site]
            if insertion:
                yield 'I', site-len(insertion), onebased_pos_in_ungapped_ref, insertion
                insertion = ''
            elif deletion:
                yield 'D', site-len(deletion), onebased_pos_in_ungapped_ref-len(deletion), deletion
                deletion = ''
                
    #clean up non match states at end of sequence
    if mismatch:
        yield 'X', site-len(mismatch), onebased_pos_in_ungapped_ref-len(mismatch), mismatch
        mismatch = ''
    if deletion:
        yield 'D', site-len(deletion), onebased_pos_in_ungapped_ref-len(deletion), deletion
        deletion = ''
    if insertion:
        yield 'I', site-len(insertion), onebased_pos_in_ungapped_ref, insertion
        insertion = ''
        
def call_mutations(mutant, reference, chromosome, ref_start=1, **kw):
    """Call mutations and print VCF formatted results
    """
    AmpliconMutation = namedtuple('AmpliconMutation', 'chromosome, start, end, vcf_start, ref_allele, alt_allele')
    result = []
    for mutation in caller(mutant, reference, **kw):
        if mutation[0] == 'X':
            #snv
            assert mutant[mutation[1]] == mutation[3]
            start = ref_start+mutation[2]-1 #subtract 1 as adding two one based coordinates
            vcf_start = start
            end = start
            ref_allele = reference[mutation[1]]
            alt_allele = mutation[3]
            
        if mutation[0] == 'D':
            #deletion
            start = ref_start+mutation[2]-1 #subtract 1 as adding two one based coordinates
            vcf_start = start - 1 #subtract 1 for context base
            end = start + len(mutation[3])-1
            ref_allele = reference[mutation[1]-1:mutation[1]+len(mutation[3])]
            alt_allele = reference[mutation[1]-1]
            
        if mutation[0] == 'I':
            #insertion
            start = ref_start+mutation[2]-1 #subtract 1 as adding two one based coordinates
            vcf_start = start - 1 #subtract 1 for context base
            end = start
            ref_allele = reference[mutation[1]-1]
            alt_allele = reference[mutation[1]-1]+mutation[3]
        result.append(AmpliconMutation(chromosome, start, end, vcf_start, ref_allele, alt_allele))
    return result

def format_as_vcf(mutations, outfile=sys.stdout):
    for mutation in mutations:
        print(mutation.chromosome,
                mutation.vcf_start,
                ".",
                mutation.ref_allele,
                mutation.alt_allele,
                '.','PASS','.',
                sep='\t',file=outfile)
    pass

def call_mutations_to_vcf(mutant, reference, chromosome, ref_start=1, outfile=sys.stdout, **kw):
    format_as_vcf(call_mutations(mutant, reference, chromosome, ref_start=ref_start, **kw), outfile=outfile)
    pass

def make_vcf_header(threshold):
    return "##fileformat=VCF4.1\n\
##source=AmBiVeRT0.1.0\n\
##FILTER=<ID=depth,Description=more than {threshold} variant supporting reads>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".format(threshold=threshold)

if __name__ == '__main__':
    pass