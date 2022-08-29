#!/usr/bin/env python3
# encoding: utf-8
"""
sequence_utilities.py

Created by Matthew Wakefield.
Copyright (c) 2013-2015  Matthew Wakefield and The University of Melbourne. All rights reserved.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

"""

import sys, os, io, re
from string import ascii_uppercase, ascii_lowercase
import gzip

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2013-2015,  Matthew Wakefield and The University of Melbourne"
__credits__ = ["Matthew Wakefield",]
__license__ = "GPLv3"
__version__ = "0.5.2"
__maintainer__ = "Matthew Wakefield"
__email__ = "matthew.wakefield@unimelb.edu.au"
__status__ = "Development"


def open_potentially_gzipped(thefile):#pragma no cover
    """Use the name of a provided file to detect gzipped files and return
    an appropriately handled file object.
    Arguments:
        thefile: a file object or filename string
    Returns:
        file object
    Note that if this opens a file it will be in 'rb' mode
    For concistency it is recommended you pass 'rb' mode files
    as input arguments.
    """
    if type(thefile) == str:
        thefile = io.open(thefile, mode='rb')
    if not hasattr(thefile,'name'):
        #likely a stringio object or other file like stream
        return thefile
    extention = os.path.splitext(thefile.name)[-1]
    if extention == '.gz':
        return gzip.open(thefile,  mode='rb')
    else:
        return thefile

def parse_fastq(thefile):
    """A basic non-validating fastq file parser
    Does not support multi-line sequence and/or multi-line quality entries.
    Arguments:
        thefile:    a binary mode file object
            (although uncompressed fastq is text, gzipped files will be 
             binary data. All files are delt with as binary for concistency.)
    Returns:
        (readname, sequence, quality)
    """
    with thefile as fastqfile:
        name = True
        while name:
            name = str(fastqfile.readline(), encoding='latin-1').strip('\n')
            seq = str(fastqfile.readline(), encoding='latin-1').strip('\n')
            fastqfile.readline() #could check starts with '+' as soft validation
            qual = str(fastqfile.readline(), encoding='latin-1').strip('\n')
            if name:
                yield name[1:], seq, qual

def parse_fasta(filename, token='>'):
    """fasta and multi-fasta file parser
    Usage: for name,seq in fasta(open(filename)):
                do something
           parse_fasta(open(filename)).next()
    Arguments:
        filename:   a read mode text file
        token:      a character that indicates the header line
                    Default: '>' (FASTA)
    """
    with filename as f:
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

def complement(seqstring):
    """Case sensitive IUPAC complement
        Arguments:
        o seqstring - an iterable object of characters (eg a string)
        
        Returns a string of complemented sequence
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
    return ''.join(result)

def reverse_complement(seqstring):
    """Case sensitive IUPAC reverse complement
        Arguments:
        o seqstring - an iterable object of characters (eg a string)
        
        Returns a string of reverse complemented sequence
    """
    return complement(seqstring)[::-1]

def encode_ambiguous(bases):
    """Encode a group of ambiguous bases with IUPAC symbol
    Arguments:
        bases:  a list, tuple or iterable of single characters
    Returns:
        IUPAC_code: a single upper case character from the set
                RYSWKMBDHVN.  Will return N if any base is N
                Will return input base if len(set(bases)) == 1
    """
    bases = tuple(sorted(set([x.upper() for x in bases])))
    if len(bases) == 1:
        return bases[0]
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
    if 'N' in bases:
        return 'N'
    else:
        return iupac[bases]

def flatten_paired_alignment(seq1,seq2,gap='-'):
    result = []
    for base1,base2 in zip(seq1,seq2):
        if base1 == gap:
            result.append(base2)
        elif base2 == gap:
            result.append(base1)
        elif base1 == base2:
            result.append(base1)
        else:
            result.append(encode_ambiguous((base1,base2)))
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

def get_tag(read,tag):
    for x in read[11:]:
        if x.startswith(tag):
            return x.split(':')[-1]
    return None

def expand_cigar(cigar):
    return "".join([x[1]*int(x[0]) for x in re.findall(r'([0-9]+)([MIDNSHPX=])',cigar)])

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
    assert len(seq) == xcigar.count(match) + xcigar.count(insert)
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

def fix_softmasked_expanded_cigar(expanded_cigar, start = 0):
    #fix SDS and SIS states
    if isinstance(expanded_cigar, list):
        expanded_cigar = ''.join(expanded_cigar)
    if 'S' in expanded_cigar:
        ##These versions trim insert delete states adjacent to softmask
        ##This produces more valid cigar strings, but may softmask variants
        #initial_SDI_block = re.findall(r'^([SDI]*)',expanded_cigar)[0]
        #final_SDI_block = re.findall(r'([SDI]*$)',expanded_cigar)[0]
        leading_matches = re.findall(r'^([SDI]*S)',expanded_cigar)
        initial_SDI_block = leading_matches[0] if leading_matches else ''
        trailing_matches = re.findall(r'(S[SDI]*$)',expanded_cigar)
        final_SDI_block = trailing_matches[0] if trailing_matches else ''
        initial_softmask_string = initial_SDI_block.replace('D','').replace('I','S')
        start += len(initial_softmask_string)
        final_softmask_string = final_SDI_block.replace('D','').replace('I','S')
        length_in_ref = len(expanded_cigar) - (len(initial_SDI_block) + len(final_SDI_block))
        fixed_expanded_cigar = initial_softmask_string + expanded_cigar[len(initial_SDI_block):len(expanded_cigar)-len(final_SDI_block)] + final_softmask_string
        ### TODO find any MSM blocks
        ### TODO replace internal S with M
    else:
        fixed_expanded_cigar = expanded_cigar
        length_in_ref = len([x for x in expanded_cigar if x in 'MD'])
    return fixed_expanded_cigar, start, length_in_ref

def gapped_alignment_to_cigar(aligned_reference,aligned_sample, gap='-', snv='M'):
    if len(aligned_reference) != len(aligned_sample):
        raise RuntimeError('Unequal sequences lengths - not correctly aligned \n    {0}\n    {1}'.format(aligned_reference,aligned_sample))
    xcigar = []
    for i in range(len(aligned_reference)):
        if aligned_reference[i] == aligned_sample[i]:
            if aligned_reference[i] == gap:
                raise RuntimeError('Invalid alignment state \n{0}\n{1}'.format(aligned_reference[i],aligned_sample[i]))
            xcigar.append('M')
        elif aligned_reference[i] in ascii_uppercase and \
            aligned_sample[i] in ascii_uppercase:
            xcigar.append(snv)
        elif aligned_reference[i] in ascii_lowercase:
            if not aligned_sample[i] == gap:
                xcigar.append('S')
            else:
                xcigar.append('D')
        elif aligned_reference[i] == gap:
            xcigar.append('I')
        elif aligned_sample[i] == gap:
            xcigar.append('D')
        else: #pragma no cover # I cant think of a way to get here but complex so dont want to default to a valid state
            raise RuntimeError('Invalid alignment state \n{0}\n{1}'.format(aligned_reference[i],aligned_sample[i]))
    fixed_xcigar, start, length = fix_softmasked_expanded_cigar(xcigar)
    return compact_cigar(fixed_xcigar), start, length
        
def expand_mdtag(mdtag):
    mdtag_tokens = re.findall(r'(\d+|\D+)',mdtag)
    result = []
    for x in mdtag_tokens:
        if x[0] in '1234567890':
            result.extend(['',]*int(x))
        elif x[0] == '^':
            result.extend(list(x.lower()))
        else:
            result.extend(list(x))
    return result

def compact_expanded_mdtag_tokens(expanded_mdtag_tokens):
    result = []
    count = 0
    in_deletion = False
    for token in expanded_mdtag_tokens:
        if token == '':
            count += 1
            in_deletion = False
        elif count:
            #exiting match
            result.append(str(count))
            count = 0
            if token == '^':
                in_deletion = True
            result.append(token)
        else:
            #in mismatch or deletion states
            if in_deletion and token == '^':
                #adjacent deletions to be merged
                continue
            if in_deletion and (token in 'CAGTN'):
                #have a mismatch adjacent to a deletion
                result.append('0')
                in_deletion = False
            result.append(token)
    if count: 
            result.append(str(count))
    return "".join(result).upper()
        


#if __name__ == '__main__':
#	main()

