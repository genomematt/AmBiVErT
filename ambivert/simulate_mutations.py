#!/usr/bin/env python
# encoding: utf-8
"""
test_mutation_mapping.py

Created by Matthew Wakefield on 2013-05-03.
Copyright (c) 2013  Matthew Wakefield and The Walter and Eliza Hall Institute. All rights reserved.
"""
from __future__ import print_function
import os, sys
import re
import argparse
from sequence_utilities import *
from truseq import parse_truseq_manifest


__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2013,  Matthew Wakefield and The Walter and Eliza Hall Institute"
__credits__ = ["Matthew Wakefield",]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Development"


def point_mutate_sequence(sequence,chromosome=None,one_based_start=None,one_based_site=1):
    site = one_based_site - 1
    reference = sequence[site]
    for base in ['A','C','G','T']:
        if base == reference:
            continue
        mdz_tag = '{pre}{mut}{post}'.format(pre = site if site else '',
                                            mut = base,
                                            post = (len(sequence)-1)-site if (len(sequence)-1)-site else '')
        name = '{chromosome}_{start}_{cigar_len}M_{mdz_tag}'.format(chromosome=chromosome,
                                                                    start=one_based_start,
                                                                    cigar_len=len(sequence),
                                                                    mdz_tag=mdz_tag)
        yield name,sequence[:site]+base+sequence[site+1:]

def make_all_point_mutations(sequence, chromosome=None, one_based_start=None, primer=0, skip_softmasked=True):
    for i in range(primer,len(sequence)):
        if sequence[i] in ['A','C','G','T']:
            for x in point_mutate_sequence(sequence, chromosome, one_based_start, i+1):
                yield x

def get_sam_header(samfile):
    line = "@"
    header = []
    pointer = 0
    while line[0] == '@':
        pointer = samfile.tell()
        line = samfile.readline().strip('\n')
        if line[0] == '@':
            header.append(line)
    samfile.seek(pointer)
    return header

def check_point_mutate_sequence(samline):
    line = samline.readline().strip('\n').split()
    true_chromosome,true_start,true_cigar,true_mdz_tag = line[0].split('_')
    tags = {x.split(":")[0]:x.split(":")[2] for x in line[11:]}
    result = []
    comparisons = [(line[2],true_chromosome),
                  (line[3],true_start),
                  (line[5],true_cigar),
                  (tags['MD'],true_mdz_tag),]
    for comparison in comparisons:
        if comparison[0] != comparison[1]:
            result.append(comparison)
    if result:
        return (result,comparisons)
    else:
        return None

def constant_quality(sequence, value="I"):
    return [value,]*len(sequence)

def format_fastq_entry(name,sequence,q_generator=constant_quality):
    return "@{name}\n{seq}\n+\n{qual}\n".format(name=name,seq=sequence,qual="".join(q_generator(sequence)))

def expand_cigar(cigar):
    return "".join([x[1]*int(x[0]) for x in re.findall(r'([0-9]+)([MIDNSHPX=])',cigar)])

def engap(seq, cigar, delete='D', insert='I', match='M', gap='-'):
    """Convert a match/delete/insert string and sequence into gapped sequence
    To convert the target sequence swap delete and insert symbols.
        Arguments:
            o   seq : a sequence string
            o   cigar : a cigar string eg 80M5D10M10I
            o   delete : deletion state character Default = 'D'
            o   insert : insertion state character Default = 'I'
            o   match : match state character Default = 'M'
            o   gap : output gap character Default = '-'
        Returns:
            o   string : gapped seqeunce
    """
    gapped = []
    xcigar = expand_cigar(cigar)
    if len(seq) != xcigar.count(match) + xcigar.count(insert):
        raise AssertionError, "Sequence length mismatch Seq:{0} Cigar:{1}".format(
                    len(seq),xcigar.count(match) + xcigar.count(insert))
    seq = list(seq)
    for symbol in xcigar:
        if symbol == delete:
            gapped.append('-')
        else:
            gapped.append(seq.pop(0))
    return "".join(gapped)

def cigar_trimmer(cigar,trim_from_start=0,trim_from_end=0):
    xcigar = expand_cigar(cigar)
    result = []
    sequence_length = len(xcigar) - xcigar.count('D')
    position_in_sequence = 0
    for state in xcigar:
        if not (state == 'D' or state == 'S'):
            position_in_sequence +=1
        if position_in_sequence > trim_from_start and position_in_sequence <= sequence_length - trim_from_end:
            result.append(state)
    return compact_cigar(result)

def compact_cigar(expanded_cigar):
    result = []
    last_state = expanded_cigar[0]
    count = 0
    for state in expanded_cigar:
        if state == last_state:
            count +=1
        else:
            result.append(str(count)+last_state)
            last_state = state
            count = 1
    result.append(str(count)+state)
    return "".join(result)

def expand_mdtag_tokens(mdtag_tokens):
    result = []
    for x in mdtag_tokens:
        if x[0] in '1234567890':
            result.extend(['',]*int(x))
        elif x[0] == '^':
            result.append(x)
        else:
            result.extend(list(x))
    return result

def compact_expanded_mdtag_tokens(expanded_mdtag_tokens):
    result = []
    count = 0
    for token in expanded_mdtag_tokens:
        if token == '':
            count += 1
            continue
        elif count:
            result.append(str(count))
            result.append(token)
            count = 0
        else:
            result.append(token)
    if count: 
            result.append(str(count))
    return "".join(result)

def mutation_detection_tag_trimmer(mdtag,trim_from_start=0,trim_from_end=0):
    #when trimming reads only non-standard state is deletions
    #mutations and insertions maintain 1 md token per string postition.
    mdtag_tokens = re.findall(r'(\d+|\D+)',mdtag)
    if len(mdtag_tokens) == 1:
        #no mutation states - just need to resize
        return str(int(mdtag_tokens[0])-(trim_from_start+trim_from_end))
    #if (mdtag_tokens[-1][0] in '1234567890') and (int(mdtag_tokens[0]) > trim_from_start) and (int(mdtag_tokens[-1]) < trim_from_end):
    #    #check that last entry is a number.  If so len(tokens) >= 3.
    #    #trim does not include mutations - resize ends
    #    new_first_token = str(int(mdtag_tokens[0])-trim_from_start)
    #    new_last_token = str(int(mdtag_tokens[-1])-trim_from_end)
    #    return "".join([new_first_token,]+mdtag_tokens[1:-1]+[new_last_token,])
    xmdtag_tokens = expand_mdtag_tokens(mdtag_tokens)
    if trim_from_end:
        return compact_expanded_mdtag_tokens(xmdtag_tokens[trim_from_start:-trim_from_end])
    else:
        return compact_expanded_mdtag_tokens(xmdtag_tokens[trim_from_start:])

def mutated_amplicon_to_paired_reads(sequence,chromosome,start,cigar,mdtag,quality='',readlength=150):
    #needs to work even when len(sequence) < readlength and return trimmed reads
    forward_sequence = sequence[:readlength]
    reverse_sequence = reverse_complement(sequence[-readlength:])
    forward_quality = quality[:readlength]
    reverse_quality = "".join(reversed(quality[-readlength:]))
    trimsize = max(0,len(sequence)-readlength)
    reverse_start = str(int(start) + trimsize)
    reverse_cigar = cigar_trimmer(cigar,trim_from_start=trimsize)
    reverse_mdtag = mutation_detection_tag_trimmer(mdtag,trim_from_start=trimsize)
    forward_start = start
    forward_cigar = cigar_trimmer(cigar,trim_from_end=trimsize)
    forward_mdtag = mutation_detection_tag_trimmer(mdtag,trim_from_end=trimsize)
    return ((forward_sequence, forward_quality, chromosome, forward_start, forward_cigar, forward_mdtag),
            (reverse_sequence, reverse_quality, chromosome, reverse_start, reverse_cigar, reverse_mdtag))

def sequence_from_fasta_file(infile):
    for name,sequence in parse_fasta(infile):
        name, chromosome, start, end = name.split()
        yield chromosome, start, sequence

def sequence_from_manifest_file(infile):
    targets = parse_truseq_manifest(infile)[2]
    for target in targets:
        yield target.Chromosome, target.Start_Position, target.Sequence

def amplicons_to_mutated_reads(forward_outfile=sys.stdout,reverse_outfile=sys.stderr,sequences=sequence_from_fasta_file(sys.stdin)):
    for chromosome, start, sequence in sequences:
        for mutant_name, mutant_sequence in list(make_all_point_mutations(sequence,chromosome=chromosome,one_based_start=start)):
                mutant_chromosome,mutant_start,cigar,mdtag = mutant_name.split('_')
                reads = mutated_amplicon_to_paired_reads(mutant_sequence,mutant_chromosome,mutant_start,cigar,mdtag)
                readname = "_".join(reads[0][2:]+reads[1][2:])
                print(format_fastq_entry(readname,reads[0][0]),end='',file=forward_outfile)
                print(format_fastq_entry(readname,reads[1][0]),end='',file=reverse_outfile)
    pass

def command_line_interface(*args,**kw):
    parser = argparse.ArgumentParser(description='A script for simulating paired end reads with mutations\
                                                 from an amplicon target file')
    parser.add_argument('--manifest',
                        type=argparse.FileType('U'),
                        default=None,
                        help='an Illumina TrueSeq Amplicon manifest file. Default: None')
    parser.add_argument('--fasta',
                        type=argparse.FileType('U'),
                        default=None,
                        help='a fasta file of amplicons with description lines ">name chromosome start end" Default: stdin')
    parser.add_argument('--read1',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='a fastq output file of forward reads. Default: stdout')
    parser.add_argument('--read2',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='a fastq output file of reverse reads. Default: stdout')
    return parser.parse_args(*args,**kw)


if __name__ == '__main__':
    args = command_line_interface()
    if args.manifest:
        amplicons_to_mutated_reads(args.read1,args.read2,sequence_from_manifest_file(args.manifest))
    elif args.fasta:
        amplicons_to_mutated_reads(args.read1,args.read2,sequence_from_fasta_file(args.fasta))
    else:
        print("You must supply a manifest or fasta file")
    
    #print(mutation_detection_tag_trimmer('140G9',trim_from_start=50,trim_from_end=0))
        
        
    #seq = "ATCAGAGATGTAGTACAACGTCGTTTCAGTCTGAGATAATCTTCTGAACTGGTGGGAGCAGTCCTAGTGGATTCACTGACAGATATAAATTGTTTTTCTCCTGTTGAACCAGACAAAA"
    #chromosome = "13"
    #one_based_start = "32972677"
    #for x in make_all_point_mutations(seq,"13","32972677"):
    #    print(x[0].split('_')[0],x[0].split('_')[1],x[0].split('_')[2],x[0].split('_')[3])
    #    print(mutated_amplicon_to_paired_reads(x[1],x[0].split('_')[0],x[0].split('_')[1],x[0].split('_')[2],x[0].split('_')[3], readlength=20))
    #    


    #print(cigar_trimmer('20M'))
    #print(cigar_trimmer('20M',trim_from_start=1 ))
    #print(cigar_trimmer('25M5D5M4I6M', trim_from_start=7, trim_from_end=6))
    ##[print(x) for x in point_mutate_sequence('GATC')]
    #print()
    #[print(x) for x in point_mutate_sequence('GATC',chromosome='test',one_based_start='1230789',one_based_site=2)]
    #print()
    #[print(x) for x in point_mutate_sequence('GATC',chromosome='test',one_based_start='1230789',one_based_site=4)]
    #print()
    #[[print(y) for y in x] for x in make_all_point_mutations('GATC',chromosome='test',one_based_start='1230789')]
    #print()
    #[[print(y) for y in x] for x in make_all_point_mutations('GATC',chromosome='test',one_based_start='1230789', primer=2)]
    #print()
    #[[print(format_fastq_entry(*y),end='') for y in x] for x in make_all_point_mutations('GATC',chromosome='test',one_based_start='1230789', primer=2)]

