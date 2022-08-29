#!/usr/bin/env python
# encoding: utf-8
"""
mpileup.py

Created by Matthew Wakefield.
Copyright (c) 2013-2015  Matthew Wakefield and The Walter and The University of Melbourne. All rights reserved.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

"""
from __future__ import print_function
import sys
import os
import argparse
from cogent.align.algorithm import nw_align

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2013-2015,  Matthew Wakefield and The University of Melbourne"
__credits__ = ["Matthew Wakefield","Graham Taylor"]
__license__ = "GPL"
__version__ = "0.5.2"
__maintainer__ = "Matthew Wakefield"
__email__ = "matthew.wakefield@unimelb.edu.au"
__status__ = "Development"

def fasta(f, token='>'):
    """fasta file generator
    Usage: for name,seq in fasta(open(filename)):"""
    if type(f) == str:
        f = open(f)
    seq = None
    name = None   
    for line in f:
        line = line.strip()
        if line.startswith(token):
            if name:
                yield (name, seq)
            seq = ''
            name = line[1:]
        elif seq != None:
            seq += line
    if name:
        yield (name, seq)


def to_mpileup_column(aligned_ref_seq, aligned_sample_seq):
    result = []
    insertion = ''
    deletion = ''
    assert len(aligned_ref_seq) == len(aligned_sample_seq)
    assert aligned_sample_seq[0] != '-' #leading gaps should be handled upstream
    for i in range(len(aligned_ref_seq)):
        if insertion != '':
            if aligned_ref_seq[i] == '-':
                #stay in insertion state
                insertion += aligned_sample_seq[i]
            else:
                #finalize insertion
                result[-1] += '+{0}{1}'.format(len(insertion),insertion)
                insertion = ''
        elif deletion != '':
            if aligned_sample_seq[i] == '-':
                #stay in deletion state
                deletion += aligned_ref_seq[i]
            else:
                #finalize deletion
                result[-1] += '-{0}{1}'.format(len(deletion),deletion)
                result.extend(['*',]*len(deletion))
                deletion = ''
    
        if not insertion and not deletion :
            if aligned_sample_seq[i] == aligned_ref_seq[i]:
                #in pairwise alignment this cant be a gap
                result.append('.')
            elif aligned_sample_seq[i] == '-':
                deletion += aligned_ref_seq[i]
            elif aligned_ref_seq[i] == '-':
                insertion += aligned_sample_seq[i]
            else:
                result.append(aligned_sample_seq[i])
    #finalize if in insertion or deletion state
    if insertion:
        result[-1] += '+{0}{1}'.format(len(insertion),insertion)
        insertion = ''
    if deletion:
        print('WARNING: Terminal deletion found! This may not produce valid mpileup format')
        #most valid response to this state is to return truncated result i.e. do nothing here
    result[0] = '^I'+result[0]
    result[-1] = result[-1]+'$'
    return result

def mpileup_from_redudndant_alignment(chromosome, start, aligned_ref_seq, aligned_sample_seq, number=1):
    while aligned_sample_seq[0] == '-':
        start +=1
        aligned_sample_seq = aligned_sample_seq[1:]
        aligned_ref_seq = aligned_ref_seq[1:]
    ungapped_ref_seq = [x for x in aligned_ref_seq if x != '-']
    mpileup_column = to_mpileup_column(aligned_ref_seq, aligned_sample_seq)
    result = []
    for i in range(len(mpileup_column)):
        result.append((chromosome,start+i,ungapped_ref_seq[i],number,''.join(mpileup_column[i]*number),'I'*number))
    return '\n'.join(['\t'.join([str(x) for x in row]) for row in result])

def test():
    print(to_mpileup_column('CACTCCATCGAGATTTCACT','CACTCCATCGAGATTTCTCT'))
    print(to_mpileup_column('CACTCCATC---ATTTCACT','CACTCCATCGAGATTTCTCT'))
    print(to_mpileup_column('CACTCCATCGAGATTTCACT','CACTCCATC---ATTTCTCT'))
    print()
    print(mpileup_from_redudndant_alignment(7,  140453119, 'CACTCCATCGAGATTTCACT','CACTCCATCGAGATTTCTCT',8))
    print()
    print(mpileup_from_redudndant_alignment(7,  140453119, 'CACTCCATC---ATTTCACT','CACTCCATCGAGATTTCTCT',8))
    print()
    print(mpileup_from_redudndant_alignment(7,  140453119, 'CACTCCATCGAGATTTCACT','CACTCCATC---ATTTCTCT',8))
    print()
    seq1,seq2=nw_align('CACTCCATCGAGATTTCACT','CACTCCATCATTTCTCT')
    print(mpileup_from_redudndant_alignment(7,  140453119, seq1,seq2,8))
    print()
    pass

def command_line_interface(*args,**kw):
    parser = argparse.ArgumentParser(description='A script for converting redundant/binned amplicon counts to mpileup format')
    parser.add_argument('--reference',
                        type=argparse.FileType('r'),
                        default=None,
                        help='a fasta format reference sequence')
    parser.add_argument('--query',
                        type=argparse.FileType('r'),
                        default=sys.stdin,
                        help='A one sequence per line tab separated file with the fields: readname Chromosome start-end\\tSequence\\tCount')    
    return parser.parse_args(*args,**kw)

if __name__ == '__main__':
    args = command_line_interface()
    #python mpileup.py --reference='BRAF_amplicon.fasta' --query='TGx120044_S6_L001_R1_001_fastq_fna_1_num_grp_flanked_positive_r_92_BRAF1_locus_1_num_grp_srt.txt',
    if args.reference:
        refname,ref_seq = fasta(args.reference).next()
    else:
        try:
            refname,ref_seq,ref_count = args.query.readline().strip('\n').split('\t')
        except:
            print("Could not read first sequence of query as reference",file=sys.stderr)
            args.print_help()
    maxcount = 0
    threshold = 0
    for name,seq,count in (x.strip('\n').split('\t') for x in args.query):
        count = int(count)
        if count > maxcount:
            maxcount = count
            threshold = 0.01*maxcount
        if count < threshold:
            print('Skipping ',name,' as is present at less than 1% of reads')
            continue
        outfile = open(name[1:].replace(' ','_')+'_mpileup','w')
        aligned_ref_seq, aligned_sample_seq = nw_align(ref_seq,seq)
        chromosome = name.split()[2].split(':')[0]
        start = int(name.split()[2].split(':')[1].split('-')[0])
        print(mpileup_from_redudndant_alignment(chromosome,  start, aligned_ref_seq, aligned_sample_seq, count),file=outfile)
        outfile.close()        
        
        

