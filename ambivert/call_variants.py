#!/usr/bin/env python2.7
# encoding: utf-8
"""
call_variants.py

Created by Matthew Wakefield and Graham Taylor.
Copyright (c) 2013-2015  Matthew Wakefield and The University of Melbourne. All rights reserved.

   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


"""
from __future__ import print_function
import sys, os
import warnings,logging
from collections import namedtuple

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2013-2015,  Matthew Wakefield and The University of Melbourne"
__credits__ = ["Matthew Wakefield","Graham Taylor"]
__license__ = "GPLv3"
__version__ = "0.5.2"
__maintainer__ = "Matthew Wakefield"
__email__ = "matthew.wakefield@unimelb.edu.au"
__status__ = "Development"

logfile = sys.stderr

def caller(variant_amplicon, reference, gap='-',softmask = True):
    """Mutation calling code
    Arguments:
        variant_amplicon: aligned query sequence as an interable of bases
                            ie a string or list of single characters
        reference:        aligned reference sequence as an iterable of bases
                          softmasked sequence is represented as lower case.
                          
    Returns:
        a list of tuples of:
            Muation type: 'I': Insertion, 'D': Deletion, 'X': Substitition
            Start:        Zero based position in gapped sequence
            Pos:          One based position in ungapped reference sequence
            Bases:        Substituted, inserted or deleted sequence
    """
    
    mismatch = ''
    deletion = ''
    insertion = ''
    onebased_pos_in_ungapped_ref = 0
    softmasked = True
    
    if len(variant_amplicon) != len(reference): #pragma no cover #works correctly but assertRaises test fails
        raise RuntimeError('Sequence lengths are not equal:\n{0}\n{1}'.format(variant_amplicon, reference))
    
    for site in range(len(variant_amplicon)):
        #ignore softmasked sites or gaps in softmasked sequence
        if softmask and softmasked and reference[site] == reference[site].lower() and not reference[site] == gap:
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
        
        if variant_amplicon[site].upper() == reference[site].upper() and not reference[site] == gap and not variant_amplicon[site] == gap:
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
            insertion += variant_amplicon[site]
            if mismatch:
                yield 'X', site-len(mismatch), onebased_pos_in_ungapped_ref-len(mismatch), mismatch
                mismatch = ''
            elif deletion:
                yield 'D', site-len(deletion), onebased_pos_in_ungapped_ref-len(deletion)+1, deletion
                deletion = ''
        elif variant_amplicon[site] == gap:
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
            mismatch += variant_amplicon[site]
            if insertion:
                yield 'I', site-len(insertion), onebased_pos_in_ungapped_ref, insertion
                insertion = ''
            elif deletion:
                yield 'D', site-len(deletion), onebased_pos_in_ungapped_ref-len(deletion), deletion
                deletion = ''
                
    #clean up non match states at end of sequence
    onebased_pos_in_ungapped_ref += 1
    if mismatch:
        yield 'X', site-len(mismatch)+1, onebased_pos_in_ungapped_ref-len(mismatch), mismatch
        mismatch = ''
    if deletion:
        yield 'D', site-len(deletion)+1, onebased_pos_in_ungapped_ref-len(deletion), deletion
        deletion = ''
    if insertion:
        yield 'I', site-len(insertion)+1, onebased_pos_in_ungapped_ref, insertion
        insertion = ''
        
def call_variants(variant_amplicon, reference, chromosome, ref_start=1, **kw):
    """Call variants and print VCF formatted results
    """
    AmpliconVariant = namedtuple('AmpliconVariant', 'chromosome, start, end, vcf_start, ref_allele, alt_allele')
    result = []
    for variant in caller(variant_amplicon, reference, **kw):
        if variant[0] == 'X':
            #snv
            assert variant_amplicon[variant[1]] == variant[3]
            start = ref_start+variant[2]-1 #subtract 1 as adding two one based coordinates
            vcf_start = start
            end = start
            ref_allele = reference[variant[1]]
            alt_allele = variant[3]
            
        if variant[0] == 'D':
            #deletion
            start = ref_start+variant[2]-1 #subtract 1 as adding two one based coordinates
            vcf_start = start - 1 #subtract 1 for context base
            end = start + len(variant[3])-1
            ref_allele = reference[variant[1]-1:variant[1]+len(variant[3])]
            alt_allele = reference[variant[1]-1]
            
        if variant[0] == 'I':
            #insertion
            start = ref_start+variant[2]-1 #subtract 1 as adding two one based coordinates
            vcf_start = start - 1 #subtract 1 for context base
            end = start
            ref_allele = reference[variant[1]-1]
            alt_allele = reference[variant[1]-1]+variant[3]
        result.append(AmpliconVariant(chromosome, start, end, vcf_start, ref_allele, alt_allele))
    return result

def format_as_vcf(variants, outfile=sys.stdout):
    for variant in variants:
        print(variant.chromosome,
                variant.vcf_start,
                ".",
                variant.ref_allele,
                variant.alt_allele,
                '.','PASS','.',
                sep='\t',file=outfile)
    pass

def call_variants_to_vcf(variant_amplicon, reference, chromosome, ref_start=1, outfile=sys.stdout, **kw):
    format_as_vcf(call_variants(variant_amplicon, reference, chromosome, ref_start=ref_start, **kw), outfile=outfile)
    pass

def make_vcf_header(threshold): #pragma no cover
    """depricated"""
    return "##fileformat=VCF4.1\n\
##source=AmBiVeRT{version}\n\
##FILTER=<ID=depth,Description=more than {threshold} variant supporting reads>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".format(version=__version__,threshold=threshold)

if __name__ == '__main__': #pragma no cover
    pass