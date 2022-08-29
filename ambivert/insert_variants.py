#!/usr/bin/env python3
# encoding: utf-8
"""
insert_variants.py

Created by Matthew Wakefield.
Copyright (c) 2013-2015  Matthew Wakefield and The University of Melbourne. All rights reserved.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

"""
from __future__ import absolute_import

import sys, os, argparse, logging, random
from .sequence_utilities import *
from .simulate_variants import expand_cigar, engap

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2013-2015,  Matthew Wakefield and The University of Melbourne"
__credits__ = ["Matthew Wakefield",]
__license__ = "GPL"
__version__ = "0.5.2"
__maintainer__ = "Matthew Wakefield"
__email__ = "matthew.wakefield@unimelb.edu.au"
__status__ = "Development"

logfile = sys.stderr
logging.basicConfig(format='%(levelname)s:AmBiVErT:%(message)s',level=logging.DEBUG)
logging.disable(logging.DEBUG)
logging.captureWarnings(True)

FORWARD_FLAGS = ['99', #    read paired, read mapped in proper pair, mate reverse strand, first in pair
                 '83', #    read paired, read mapped in proper pair, read reverse strand, first in pair
                 '73', #    read paired, mate unmapped, first in pair
                 '89', #    read paired, mate unmapped, first in pair, read reverse strand
                 ]
REVERSE_FLAGS = ['147',#    read paired, read mapped in proper pair, read reverse strand, second in pair
                 '163',#    read paired, read mapped in proper pair, mate reverse strand, second in pair
                 '137',#    read paired, mate unmapped, second in pair
                 '153', #    read paired, mate unmapped, second in pair, read reverse strand
                 ]
FORWARD_PAIR_FLAGS = ['99','83']
REVERSE_PAIR_FLAGS = ['147','163']
FORWARD_SINGLETON_FLAGS = ['73','89']
REVERSE_SINGLETON_FLAGS = ['137','153']
MINUS_STRAND_FLAGS = ['81','83','89','163','153','181','145','147']

def split_location_string(location):
    chromosome = location.split(':')[0]
    start = location.split(':')[-1].split('-')[0]
    end = location.split(':')[-1].split('-')[-1] #may be same as start
    return chromosome, start, end

def format_fastq(name,seq,qual,minus_strand = False,end='\n'):
    if minus_strand:
        return '@{0}\n{1}\n+\n{2}{3}'.format(name,reverse_complement(seq),qual[::-1],end)
    else:
        return '@{0}\n{1}\n+\n{2}{3}'.format(name,seq,qual,end)

def overlaps(pos,length,start,end):
    return not ( (int(pos)+int(length)-1 < int(start)) or (int(pos) > int(end)) )
    
def get_readnames_overlapping_position(samfile, chromosome, start, end, readlength=300):
    first_read = samfile.tell()
    
    read_names = []
    
    for line in samfile:
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.split('\t')[:11]
        if rname == chromosome:
            xcigar = str(expand_cigar(cigar))
            length_in_reference = len(seq) + xcigar.count('D') - xcigar.count('I')
            if overlaps(pos,length_in_reference,start,end):
                read_names.append(qname)
    
    samfile.seek(first_read)
    
    return read_names

def get_readpairs_overlapping_position(samfile, chromosome, start, end, read_names=None):
    """returns overlapping pairs with two passes of the file.
    First pass identifies all matching reads.
    Second pass collects reads from matching pairs.
    """
        
    samfile.seek(0)
    header = get_sam_header(samfile)
    
    if not read_names:
        read_names = get_readnames_overlapping_position(samfile,chromosome,start,end)
    forward_reads = {}
    reverse_reads = {}
    
    for line in samfile:
        if line.split('\t')[0] in read_names:
            qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.split('\t')[:11]
            minus = flag in MINUS_STRAND_FLAGS
            if flag in FORWARD_FLAGS:
                forward_reads[qname] = (pos,cigar,seq,qual,minus)
            if flag in REVERSE_FLAGS:
                reverse_reads[qname] = (pos,cigar,seq,qual,minus)
    
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

def get_reads_not_overlapping_position(samfile, chromosome, start, end, read_names=None):
    # TODO check memory usage and speed of this function on our data sizes
    # Is a two or one dictionary based version faster?
    samfile.seek(0)
    header = get_sam_header(samfile)
    
    if not read_names:
        read_names = get_readnames_overlapping_position(samfile,chromosome,start,end)
    
    forward_reads = []
    reverse_reads = []
    singleton_reads = []
    
    for line in samfile:
        if not line.split('\t')[0] in read_names:
            qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.split('\t')[:11]
            # these flags need to be selected so we only get each read once.  May need to use TAGs as well?
            minus = flag in MINUS_STRAND_FLAGS
            if flag in FORWARD_PAIR_FLAGS: #maps to plus strand and has correctly mapped mate
                forward_reads.append((qname,seq,qual,minus))
            elif flag in REVERSE_PAIR_FLAGS: #maps to minus strand and has correctly mapped mate
                reverse_reads.append((qname,seq,qual,minus))
            elif flag in FORWARD_SINGLETON_FLAGS:   #is a singleton
                singleton_reads.append((qname,seq,qual,minus))
            elif flag in REVERSE_SINGLETON_FLAGS:   #is a singleton
                singleton_reads.append((qname,seq,qual,minus))
    
                
    assert len(forward_reads) == len(reverse_reads) #should be same names in each therefore same length
    
    forward_read_names = [x[0] for x in forward_reads] #these two lines may be slow. Worth checking?
    assert len(set(forward_read_names)) == len(forward_read_names) #test each only present once
    
    forward_reads.sort()
    reverse_reads.sort()
    
    return forward_reads, reverse_reads, singleton_reads


def get_variant_position(pos,cigar,one_based_variant_site):
    pos = int(pos)
    if one_based_variant_site < pos:
        return one_based_variant_site-pos
    xcigar = str(expand_cigar(cigar))
    if xcigar.count('D') or xcigar.count('I') or xcigar.count('S'):
        pos_in_ref = pos
        pos_in_read = 0
        pos_in_cigar = 0
        while pos_in_ref != one_based_variant_site and pos_in_cigar < len(xcigar):
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
            if pos_in_cigar >= len(xcigar):
                return max(len(xcigar)+1,one_based_variant_site-pos) #ensure > length of read
        assert pos_in_ref == one_based_variant_site
        return pos_in_read
    else:
        return one_based_variant_site-pos

def point_mutate_read(name,pos,cigar,seq,qual,minus,one_based_variant_site,variant_base):
    try:
        variant_site_in_read = get_variant_position(pos,cigar,one_based_variant_site)
    except:
        print((pos,cigar,one_based_variant_site))
        raise
    
    if variant_site_in_read >= len(seq) or variant_site_in_read < 0:
        # postion does not exist in this read so return unmodified
        # occurs when this function is applied to a pair and only one read overlaps mutation
        return name,seq,qual,minus
        
    seq = list(seq)
    try:
        seq[variant_site_in_read] = variant_base
    except:
        print("".join(seq))
        print(' '*variant_site_in_read+'^',variant_site_in_read)
        raise
    
    return name,"".join(seq),qual,minus

def deletion_mutate_read(name,pos,cigar,seq,qual,minus,one_based_variant_site,length):
    variant_start_site_in_read = get_variant_position(pos,cigar,one_based_variant_site)
    variant_end_site_in_read = get_variant_position(pos,cigar,one_based_variant_site+length)
    
    seq = str(seq)
    qual = str(qual)
    
    # String slicing out of range returns correct strings.
    # No specific code needed to deal with mutation outside of sequence
    mut_seq = seq[:variant_start_site_in_read]+seq[variant_end_site_in_read:]
    mut_qual = qual[:variant_start_site_in_read]+qual[variant_end_site_in_read:]
    
    return name, mut_seq, mut_qual,minus

def insertion_mutate_read(name,pos,cigar,seq,qual,minus,one_based_variant_site,insertion_bases,insertion_quality=None):
    variant_site_in_read = get_variant_position(pos,cigar,one_based_variant_site)
    
    seq = str(seq)
    qual = str(qual)
    
    if variant_site_in_read > len(seq):
        # postion does not exist in this read so return unmodified
        # This can occur when this function is applied to pairs
        return name,seq,qual,minus
    
    mut_seq = seq[:variant_site_in_read]+insertion_bases+seq[variant_site_in_read:]
    
    if not insertion_quality:
        insertion_quality = 'I'*len(insertion_bases)
    mut_qual = qual[:variant_site_in_read]+insertion_quality+qual[variant_site_in_read:]
    
    return name, mut_seq[:len(seq)], mut_qual[:len(seq)],minus

    
def salt_variant_into_mapped_reads(samfile, chromosome, one_based_variant_site,
                                    variant_base='', variant_length=0, insertion_bases='',
                                    forward_fastq_name='mutated_R1.fastq', reverse_fastq_name='mutated_R2.fastq',
                                    unpaired_fastq_name='mutated_unpaired.fastq',
                                    variant_frequency=0.5,
                                    variant_name='',
                                    overwrite=False):
    if (os.path.exists(forward_fastq_name) or os.path.exists(reverse_fastq_name)) and not overwrite:
        logging.info('Skipping {0} as file {1} or {2} already exists'.format(variant_name,forward_fastq_name,reverse_fastq_name))
        return None
        
    forward_fastq = open(forward_fastq_name,'w')
    reverse_fastq = open(reverse_fastq_name,'w')
    unpaired_fastq = open(unpaired_fastq_name,'w')
    
    length = max(len(variant_base), int(variant_length), len(insertion_bases))
    
    read_names = get_readnames_overlapping_position(samfile, chromosome, one_based_variant_site, int(one_based_variant_site)+length)
    forward_reads, reverse_reads, singleton_reads = get_reads_not_overlapping_position(samfile,chromosome, one_based_variant_site, int(one_based_variant_site)+length, read_names=read_names)
    logging.info('{0} read pairs do not overlap the mutation site'.format(len(forward_reads)))
    for read in forward_reads:
        forward_fastq.write(format_fastq(*read))
    for read in reverse_reads:
        reverse_fastq.write(format_fastq(*read))
    for read in singleton_reads:
        unpaired_fastq.write(format_fastq(*read))
    
    templates = get_readpairs_overlapping_position(samfile,chromosome, one_based_variant_site, one_based_variant_site+length,read_names=read_names)
    logging.info('{0} read pairs overlap the mutation site'.format(len(templates)))
    
    total_templates = len(templates)
    number_of_templates_to_mutate = min(int(variant_frequency * total_templates) + 1, total_templates)
    logging.info('inserting mutations into {0} read pairs'.format(number_of_templates_to_mutate))
    
    
    random.shuffle(templates)
    
    templates_to_mutate = templates[:number_of_templates_to_mutate]
    
    if len(variant_base) == 1: #SNV
        for i,template in enumerate(templates_to_mutate):
            name = "{variant_name}_{index}_{variants}_{total}".format(variant_name=variant_name,index=i,variants=number_of_templates_to_mutate,total=total_templates)
            if template[0] and template[1]:
                pos,cigar,seq,qual,minus = template[0][1:6]
                forward_fastq.write(format_fastq(*point_mutate_read(name,pos,cigar,seq,qual,minus,one_based_variant_site,variant_base)))
                pos,cigar,seq,qual,minus = template[1][1:6]
                reverse_fastq.write(format_fastq(*point_mutate_read(name,pos,cigar,seq,qual,minus,one_based_variant_site,variant_base)))
            elif template[0]:
                pos,cigar,seq,qual,minus = template[0][1:6]
                unpaired_fastq.write(format_fastq(*point_mutate_read(name,pos,cigar,seq,qual,minus,one_based_variant_site,variant_base)))
            elif template[1]:
                pos,cigar,seq,qual,minus = template[1][1:6]
                unpaired_fastq.write(format_fastq(*point_mutate_read(name,pos,cigar,seq,qual,minus,one_based_variant_site,variant_base)))
        
    elif variant_length:
        for i,template in enumerate(templates_to_mutate):
            name = "{variant_name}_{index}_{variants}_{total}".format(variant_name=variant_name,index=i,variants=number_of_templates_to_mutate,total=total_templates)
            if template[0] and template[1]:
                pos,cigar,seq,qual,minus = template[0][1:6]
                forward_fastq.write(format_fastq(*deletion_mutate_read(name,pos,cigar,seq,qual,minus,one_based_variant_site,variant_length)))
                pos,cigar,seq,qual,minus = template[1][1:6]
                reverse_fastq.write(format_fastq(*deletion_mutate_read(name,pos,cigar,seq,qual,minus,one_based_variant_site,variant_length)))
            elif template[0]:
                pos,cigar,seq,qual,minus = template[0][1:6]
                unpaired_fastq.write(format_fastq(*deletion_mutate_read(name,pos,cigar,seq,qual,minus,one_based_variant_site,variant_length)))
            elif template[1]:
                pos,cigar,seq,qual,minus = template[1][1:6]
                unpaired_fastq.write(format_fastq(*deletion_mutate_read(name,pos,cigar,seq,qual,minus,one_based_variant_site,variant_length)))
        
    elif insertion_bases and not variant_length:
        for i,template in enumerate(templates_to_mutate):
            name = "{variant_name}_{index}_{variants}_{total}".format(variant_name=variant_name,index=i,variants=number_of_templates_to_mutate,total=total_templates)
            if template[0] and template[1]:
                pos,cigar,seq,qual,minus = template[0][1:6]
                forward_fastq.write(format_fastq(*insertion_mutate_read(name,pos,cigar,seq,qual,minus,one_based_variant_site,insertion_bases)))
                pos,cigar,seq,qual,minus = template[1][1:6]
                reverse_fastq.write(format_fastq(*insertion_mutate_read(name,pos,cigar,seq,qual,minus,one_based_variant_site,insertion_bases)))
            elif template[0]:
                pos,cigar,seq,qual,minus = template[0][1:6]
                unpaired_fastq.write(format_fastq(*insertion_mutate_read(name,pos,cigar,seq,qual,minus,one_based_variant_site,insertion_bases)))
            elif template[1]:
                pos,cigar,seq,qual,minus = template[1][1:6]
                unpaired_fastq.write(format_fastq(*insertion_mutate_read(name,pos,cigar,seq,qual,minus,one_based_variant_site,insertion_bases)))

    elif insertion_bases and variant_length:
        #this twostep implementation will trim len(deletion) - len(insertion) bases from end
        for i,template in enumerate(templates_to_mutate):
            name = "{variant_name}_{index}_{variants}_{total}".format(variant_name=variant_name,index=i,variants=number_of_templates_to_mutate,total=total_templates)
            if template[0] and template[1]:
                pos,cigar,seq,qual,minus = template[0][1:6]
                name,delseq,delqual = deletion_mutate_read(name,pos,cigar,seq,qual,minus,one_based_variant_site,variant_length)
                forward_fastq.write(format_fastq(*insertion_mutate_read(name,pos,cigar,delseq,delqual,minus,one_based_variant_site,insertion_bases)))
                pos,cigar,seq,qual,minus = template[1][1:6]
                name,delseq,delqual = deletion_mutate_read(name,pos,cigar,seq,qual,minus,one_based_variant_site,variant_length)
                reverse_fastq.write(format_fastq(*insertion_mutate_read(name,pos,cigar,delseq,delqual,minus,one_based_variant_site,insertion_bases)))
            elif template[0]:
                pos,cigar,seq,qual,minus = template[0][1:6]
                name,delseq,delqual = deletion_mutate_read(name,pos,cigar,seq,qual,minus,one_based_variant_site,variant_length)
                unpaired_fastq.write(format_fastq(*insertion_mutate_read(name,pos,cigar,delseq,delqual,minus,one_based_variant_site,insertion_bases)))
            elif template[1]:
                pos,cigar,seq,qual,minus = template[1][1:6]
                name,delseq,delqual = deletion_mutate_read(name,pos,cigar,seq,qual,minus,one_based_variant_site,variant_length)
                unpaired_fastq.write(format_fastq(*insertion_mutate_read(name,pos,cigar,delseq,delqual,minus,one_based_variant_site,insertion_bases)))
    
    else: #pragma no cover
        raise RuntimeError('Unexpected state in insert variants')
    
    for template in templates[number_of_templates_to_mutate:]:
            if template[0] and template[1]:
                forward_fastq.write(format_fastq(template[0][0],template[0][3],template[0][4],template[0][5]))
                reverse_fastq.write(format_fastq(template[1][0],template[1][3],template[1][4],template[1][5]))
            elif template[0]:
                unpaired_fastq.write(format_fastq(template[0][0],template[0][3],template[0][4],template[0][5]))
            elif template[1]:
                unpaired_fastq.write(format_fastq(template[1][0],template[1][3],template[1][4],template[1][5]))
    pass

def get_mutations_from_vcf(vcffile):
    """Lightweight VCF parser that yields CHROM POS ID REF ALT"""
    for line in vcffile:
        if line and line[0] == '#':
            continue
        line=line.split('\t')
        yield line[0],int(line[1]),line[2],line[3],line[4]

def confirm_mutations_from_vcf(vcffile,variant_frequency=0.5):
    """Extracts mutations from a vcf file and looks for matching mutations in files named <identifier>_<frequency>.vcf"""
    if hasattr(vcffile,'name'):
        logging.info('Reading mutations to confirm from {0}'.format(vcffile.name))
    i=0
    correct = 0
    incorrect = 0
    untested = 0
    for chrom,one_based_variant_site,ident,ref,alt in get_mutations_from_vcf(vcffile):
        #logging.info('Mutation to confirm {4} {0}:{1} {2}>{3}'.format(chrom,one_based_variant_site,ref,alt,ident))
        i+=1
        if ident == '.':
            ident = str(i)
        testvcfname = '{0}_{1}.vcf'.format(ident,variant_frequency)
        if not os.path.exists(testvcfname):
            untested += 1
            #logging.info('Could not find a test vcf file for {4} {0}:{1} {2}>{3}'.format(chrom,one_based_variant_site,ref,alt,ident))
            continue
        found = False
        with open(testvcfname,'rt') as testvcffile:
            close_calls = []
            for test_entry in get_mutations_from_vcf(testvcffile):
                if test_entry == (chrom,one_based_variant_site,'.',ref,alt):
                    found = True
                    correct += 1
                    logging.info('Correctly identified {4} {0}:{1} {2}>{3}'.format(chrom,one_based_variant_site,ref,alt,ident))
                elif test_entry[0] == chrom and int(one_based_variant_site)-int(test_entry[1]) > 0 and int(one_based_variant_site)-int(test_entry[1]) < 20:
                    close_calls.append(test_entry)
            if found == False:
                incorrect += 1
                logging.info('Failed to identify {4} {0}:{1} {2}>{3}'.format(chrom,one_based_variant_site,ref,alt,ident))
                if close_calls: [logging.info('    {0}'.format(x)) for x in close_calls]
                
                
    print('Identifed correctly: {0}\nIncorrect: {1}\nUntested: {2}'.format(correct, incorrect, untested))
    return (correct, incorrect, untested)
            
def salt_mutations_from_vcf(vcffile,samfile,variant_frequency=0.5):
    if hasattr(vcffile,'name'):
        logging.info('Reading mutations to insert from {0}'.format(vcffile.name))
    i=0
    for chrom,one_based_variant_site,ident,ref,alt in get_mutations_from_vcf(vcffile):
        logging.info('Mutation to insert {4} {0}:{1} {2}>{3}'.format(chrom,one_based_variant_site,ref,alt,ident))
        i+=1
        if ident == '.':
            ident = str(i)
        if len(ref) == 1 and len(alt) == 1:
            #SNV
            salt_variant_into_mapped_reads(samfile, chrom, one_based_variant_site,
                                                variant_base=alt, variant_length=0, insertion_bases='',
                                                forward_fastq_name='{0}_{1}_R1.fastq'.format(ident,variant_frequency),
                                                reverse_fastq_name='{0}_{1}_R2.fastq'.format(ident,variant_frequency),
                                                unpaired_fastq_name='{0}_{1}_unpaired.fastq'.format(ident,variant_frequency),
                                                variant_frequency=variant_frequency,
                                                variant_name='')
        elif len(ref) > len(alt):
            #deletion
            if len(alt) == 1:
                #simple deletion
                salt_variant_into_mapped_reads(samfile, chrom, one_based_variant_site+1,
                                                variant_base='', variant_length=len(ref) - 1, insertion_bases='',
                                                forward_fastq_name='{0}_{1}_R1.fastq'.format(ident,variant_frequency),
                                                reverse_fastq_name='{0}_{1}_R2.fastq'.format(ident,variant_frequency),
                                                unpaired_fastq_name='{0}_{1}_unpaired.fastq'.format(ident,variant_frequency),
                                                variant_frequency=variant_frequency,
                                                variant_name='')
            else:
                #further checks here for context other than single base
                #compound insertion deletion
                salt_variant_into_mapped_reads(samfile, chrom, one_based_variant_site+1,
                                                variant_base='', variant_length=len(ref) - 1, insertion_bases=alt[1:],
                                                forward_fastq_name='{0}_{1}_R1.fastq'.format(ident,variant_frequency),
                                                reverse_fastq_name='{0}_{1}_R2.fastq'.format(ident,variant_frequency),
                                                unpaired_fastq_name='{0}_{1}_unpaired.fastq'.format(ident,variant_frequency),
                                                variant_frequency=variant_frequency,
                                                variant_name='')
        elif len(alt) >= len(ref):
            if len(ref) == 1:
                #simple insertion
                salt_variant_into_mapped_reads(samfile, chrom, one_based_variant_site+1,
                                                variant_base='', variant_length=0, insertion_bases=alt[1:],
                                                forward_fastq_name='{0}_{1}_R1.fastq'.format(ident,variant_frequency),
                                                reverse_fastq_name='{0}_{1}_R2.fastq'.format(ident,variant_frequency),
                                                unpaired_fastq_name='{0}_{1}_unpaired.fastq'.format(ident,variant_frequency),
                                                variant_frequency=variant_frequency,
                                                variant_name='')
            else:
                #further checks here for context other than single base
                #compound insertion deletion
                salt_variant_into_mapped_reads(samfile, chrom, one_based_variant_site+1,
                                                variant_base='', variant_length=len(ref) - 1, insertion_bases=alt[1:],
                                                forward_fastq_name='{0}_{1}_R1.fastq'.format(ident,variant_frequency),
                                                reverse_fastq_name='{0}_{1}_R2.fastq'.format(ident,variant_frequency),
                                                unpaired_fastq_name='{0}_{1}_unpaired.fastq'.format(ident,variant_frequency),
                                                variant_frequency=variant_frequency,
                                                variant_name='')
        else:
            raise RuntimeError('Unexpected state when parsing VCF {0}'.format((chrom,one_based_variant_site,ident,ref,alt)))

def process_commandline_args(): #pragma no cover
    """Processes command line arguments and returns args object
    Use insert_variants --help for details of arguments
    """
    parser = argparse.ArgumentParser(description="""A script for inserting
        variants from a vcffile into a mapped reads.
        Intended for making test data sets with known mutations in a
        background of real technical variability.
        
        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
        
        This software is provided for research use only.
        """)
    parser.add_argument('--sam',
                        type=argparse.FileType('rt'),
                        default=sys.stdin,
                        help='a sam file or sam stream.  Default: stdin\
                            Use "samtools view <yourfilename> | insert_variants" for bam files')    
    parser.add_argument('--vcf',
                        type=argparse.FileType('rt'),
                        help='a vcf file')
    parser.add_argument('--freq',
                        type=float,
                        default=0.5,
                        help='the frequency of the inserted mutation')    
    parser.add_argument('--confirm',
                        action='store_true',
                        help='check the output of a caller run on salted data.\
                            You must specify the the same vcf file and frequency as the original salting run.\
                            Looks for files named <identifier>_<frequency>.vcf')    
    args = parser.parse_args()
    if not args.vcf:
        raise RuntimeError('You must provide a VCF file')
    return args

def main():
    args = process_commandline_args()
    if args.confirm:
        confirm_mutations_from_vcf(args.vcf,variant_frequency=args.freq)
    else:
        salt_mutations_from_vcf(vcffile=args.vcf,samfile=args.sam,variant_frequency=args.freq)
    pass

if __name__ == '__main__':
    main()

