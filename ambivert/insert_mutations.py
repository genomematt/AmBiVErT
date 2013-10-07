#!/usr/bin/env python
# encoding: utf-8
"""
insert_mutations.py

Created by Matthew Wakefield on 2013-09-26.
Copyright (c) 2013  Matthew Wakefield and The University of Melbourne. All rights reserved.
"""
from __future__ import print_function
import sys
import os
from sequence_utilities import *
from simulate_mutations import expand_cigar, engap

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2013,  Matthew Wakefield and The University of Melbourne"
__credits__ = ["Matthew Wakefield",]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Matthew Wakefield"
__email__ = "matthew.wakefield@unimelb.edu.au"
__status__ = "Development"

FORWARD_FLAGS = ['99',]
REVERSE_FLAGS = ['147',]

def split_location_string(location):
    chromosome = location.split(':')[0]
    start = location.split(':')[-1].split('-')[0]
    end = location.split(':')[-1].split('-')[-1] #may be same as start
    return chromosome, start, end

def overlaps(pos,length,start,end):
    return not ( (int(pos)+int(length)-1 < int(start)) or (int(pos) > int(end)) )
    
def get_readnames_overlapping_position(samfile, chromosome, start, end):
    first_read = samfile.tell()
    
    read_names = []
    
    for line in samfile:
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.split('\t')[:11]
        if rname == chromosome:
            xcigar = str(expand_cigar(cigar))
            length_in_reference = len(seq) + xcigar.count('D') - xcigar.count('I')
            if overlaps(pos,len(seq),start,end):
                read_names.append(qname)
    
    samfile.seek(first_read)
    
    return read_names

def get_readpairs_overlapping_position(samfile, chromosome, start, end):
    """returns overlapping pairs with two passes of the file.
    First pass identifies all matching reads.
    Second pass collects reads from matching pairs.
    """
        
    samfile.seek(0)
    header = get_sam_header(samfile)
    
    read_names = get_readnames_overlapping_position(samfile,chromosome,start,end)
    forward_reads = {}
    reverse_reads = {}
    
    for line in samfile:
        if line.split('\t')[0] in read_names:
            qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.split('\t')[:11]
            if flag in FORWARD_FLAGS:
                forward_reads[qname] = (pos,cigar,seq,qual)
            if flag in REVERSE_FLAGS:
                reverse_reads[qname] = (pos,cigar,seq,qual)
    
    result = []
    
    for qname in forward_reads:
        if qname in reverse_reads:
            result.append(((qname,)+forward_reads[qname],(qname,)+reverse_reads[qname]))
        else:
            result.append(((qname,)+forward_reads[qname],None))
    
    for qname in reverse_reads:
        if not qname in forward_reads:
            result.append((None,(qname,)+reverse_reads[qname]))
    
    return result

def get_reads_not_overlapping_position(samfile, chromosome, start, end):
    # TODO check memory usage and speed of this function on our data sizes
    # Is a two or one dictionary based version faster?
    samfile.seek(0)
    header = get_sam_header(samfile)
    
    read_names = get_readnames_overlapping_position(samfile,chromosome,start,end)
    
    for line in samfile:
        if not line.split('\t')[0] in read_names:
            qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.split('\t')[:11]
            if flag in FORWARD_FLAGS:
                yield (qname,seq,qual)
            if flag in REVERSE_FLAGS:
                yield (qname, reverse_complement(seq),qual[::-1])

def point_mutate_read_pair(read_pair,):
    if read_pair[0]:
        yield mutate_read(*read_pair[0])
    if read_pair[1]:
        name,seq,qual = mutate_read(*read_pair[1])
        yield name, reverse_complement(seq), qual[::-1]

def get_mutation_position(pos,cigar,one_based_mutation_site):
    pos = int(pos)
    xcigar = str(expand_cigar(cigar))
    if xcigar.count('D') or xcigar.count('I') or xcigar.count('S'):
        pos_in_ref = pos
        pos_in_read = 0
        pos_in_cigar = 0
        while pos_in_ref != one_based_mutation_site:
            if xcigar[pos_in_cigar] == 'M':
                pos_in_read += 1
                pos_in_ref += 1
            elif xcigar[pos_in_cigar] in ['I', 'S']:
                pos_in_read += 1
                #extra read base dont so incriment pos_in_ref
            elif xcigar[pos_in_cigar] == 'D':
                pos_in_ref += 1
                # missing read base dont so incriment pos_in_read
                # mutated bases will occur on first non deletion
                # base before deletion if they overlap
                # TODO change this behaviour so deletion reduced.
            pos_in_cigar += 1
            assert pos_in_cigar <= len(xcigar)
        assert pos_in_ref == one_based_mutation_site
        return pos_in_read
    else:
        return one_based_mutation_site-pos

def point_mutate_read(name,pos,cigar,seq,qual,one_based_mutation_site,mutation_base):
    mutation_site_in_read = get_mutation_position(pos,cigar,one_based_mutation_site)
    
    seq = list(seq)
    
    seq[mutation_site_in_read] = mutation_base
    
    return name,"".join(seq),qual

def deletion_mutate_read(name,pos,cigar,seq,qual,one_based_mutation_site,length):
    mutation_start_site_in_read = get_mutation_position(pos,cigar,one_based_mutation_site)
    mutation_end_site_in_read = get_mutation_position(pos,cigar,one_based_mutation_site+length)
    
    seq = str(seq)
    qual = str(qual)
    
    mut_seq = seq[:mutation_start_site_in_read]+seq[mutation_end_site_in_read:]
    mut_qual = qual[:mutation_start_site_in_read]+qual[mutation_end_site_in_read:]
    
    return name, mut_seq, mut_qual

    
    

if __name__ == '__main__':
	main()

