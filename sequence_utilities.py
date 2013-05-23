#!/usr/bin/env python
# encoding: utf-8
"""
sequence_utilities.py

Created by Matthew Wakefield on 2013-05-02.
Copyright (c) 2013  Matthew Wakefield and The University of Melbourne. All rights reserved.
"""
from __future__ import print_function
import sys, os
import itertools

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2013,  Matthew Wakefield and The University of Melbourne"
__credits__ = ["Matthew Wakefield",]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Matthew Wakefield"
__email__ = "matthew.wakefield@unimelb.edu.au"
__status__ = "Development"


def open_compressed(filename, mode='r'):
    extention = os.path.splitext(filename)
    if extention == 'bz2':
        return BZ2File(filename, mode)
    elif extention == 'gz':
        return GzipFile(filename, mode)
    else:
        return open(filename, mode)

def parse_fastq(filename):
    with open_compressed(filename) as fastqfile:
        name = True
        while name:
            name = fastqfile.readline().strip('\n')
            seq = fastqfile.readline().strip('\n')
            fastqfile.readline()
            qual = fastqfile.readline().strip('\n')
            if name:
                yield name[1:], seq, qual

def parse_fasta(filename, token='>'):
    """fasta and multi-fasta file parser
    Usage: for name,seq in fasta(open(filename)):
                do something
           parse_fasta(open(filename)).next()
    """
    with open_compressed(filename) as f:
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

def reverse_complement(seqstring):
    """Case sensitive IUPAC complement
        Arguments:
        o seqstring - an iterable object of characters (eg a string)
        o reverse   - optional boolean flag for reversing order of returned result default=False
        
        Returns a string
    """
    result=[]
    complement={'a':'t', 'c':'g', 'g':'c', 't':'a',
                'A':'T', 'C':'G', 'G':'C', 'T':'A',
                '-':'-', '.':'.', 'n':'n', 'N':'N',
                'U':'A', 'M':'K', 'R':'Y', 'W':'W',
                'S':'S', 'Y':'R', 'K':'M', 'V':'B',
                'H':'D', 'D':'H', 'B':'V',
                'u':'a', 'm':'k', 'r':'y', 'w':'w',
                's':'s', 'y':'r', 'k':'m', 'v':'b',
                'h':'d', 'd':'h', 'b':'v',
                }
    for base in seqstring:
        try:
            result.append(complement[base])
        except KeyError:
            result.append('n')
    return ''.join(result[::-1])

def encode_ambiguous(bases):
    bases = tuple(set([x.upper() for x in bases]))
    iupac = {('A','G'):'R',
            ('C','T'):'Y',
            ('C','G'):'S',
            ('A','T'):'W',
            ('G','T'):'K',
            ('A','C'):'M',
            ('C','G','T'):'B',
            ('A','G','T'):'D',
            ('A','C','T'):'H',
            ('A','C','G'):'V',
            ('A','C','G','T'):'N',}
    return iupac[bases]

def flatten_paired_alignment(seq1,seq2,gap='-'):
    result = []
    for base1,base2 in itertools.izip(seq1,seq2):
        if base1 == gap:
            result.append(base2)
        elif base2 == gap:
            result.append(base1)
        elif base1 == base2:
            result.append(base1)
        else:
            result.append(encode_ambiguous(base1,base2))
    return ''.join(result)

def make_blocklist(seqstring, block_size=80):
    """format sequence into a list of blocks"""
    blocklist = []
    seqlength = len(seqstring)
    for block in range(0, seqlength, block_size):
        if block + block_size < seqlength:
            blocklist.append(seqstring[block: block + block_size])
        else:
            blocklist.append(seqstring[block:])
    return blocklist

def slice_string_in_blocks(seqstring, block_size=80):
    """slice string into block_size[=80] lines"""
    blocklist = make_blocklist(seqstring, block_size)
    return '\n'.join(blocklist) + '\n'

def format_fasta(name,seq, block_size=80):
    """returns a string in fasta format"""
    return '>'+name+'\n'+slice_string_in_blocks(seq,block_size)

def format_alignment(sequences):
    align={}
    result = ''
    for (name, seq) in sequences:
        align[name]=make_blocklist(seq, block_size=50)
    for i in range(len(align.values()[0])):
        for key in align:
            result +='%30s\t' % key
            if align[key][i]:
                result += align[key][i]+'\n'
        result += '\n'
    return result

#if __name__ == '__main__':
#	main()

