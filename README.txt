Experimental Alpha Release of AmBiVErT
--------------------------------------

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


Created by Matthew Wakefield and Graham Taylor.
Copyright (c) 2013  Matthew Wakefield and The University of Melbourne. All rights reserved.

$ ambivert --help
usage: ambivert [-h] [-f FORWARD] [-r REVERSE] [-m MANIFEST] [--fasta FASTA]
                [--output OUTPUT] [--countfile COUNTFILE]
                [--threshold THRESHOLD] [--min_cover MIN_COVER]
                [--min_reads MIN_READS] [--min_freq MIN_FREQ]
                [--overlap OVERLAP] [--hashtable HASHTABLE]
                [--savehashtable SAVEHASHTABLE] [--alignments ALIGNMENTS]
                [--prefix PREFIX]

AmBiVErT: A program for binned analysis of amplicon data AmBiVErT clusters
identical amplicon sequences and thresholds based on read frequency to remove
technical errors. Due to sequencing errors occuring with a more random
distribution than low frequency variants this approach reduces the number of
amplicon products that must be assigned to target regions & assessed for
variant calls. AmBiVErT overlaps forward and reverse reads from the same
amplicon and preserves local phasing information. Typical running time for
first use is several hours, which reduces to less than 10 minutes when the
hash table calculated on a previous run is supplied for analysis of subsequent
samples with the same amplicons.

optional arguments:
  -h, --help            show this help message and exit
  -f FORWARD, --forward FORWARD
                        a fastq format file of forward direction amplicon
                        reads. May be compressed with gzip with the
                        appropriate suffix (.gz)
  -r REVERSE, --reverse REVERSE
                        a fastq format file of reverse direction amplicon
                        reads. May be compressed with gzip with the
                        appropriate suffix (.gz)
  -m MANIFEST, --manifest MANIFEST
                        an Illumina TrueSeq Amplicon manifest file.
  --fasta FASTA         an fasta amplicon manifest file. Sequences should be
                        limited to the regions to be called and exclude
                        primers & adaptors. This file can be provided in
                        addition to an Illumina manifest to specify additional
                        off target regions.
  --output OUTPUT       output of alignments with variants. Default: stdout
  --countfile COUNTFILE
                        output of occurance counts per amplicon. Includes all
                        counts that map to the reference amplicon. This count
                        does not include reads that occured at frequencies
                        below --threshold <default=20>
  --threshold THRESHOLD
                        the minimum occurance threshold. Unique amplicon
                        sequence variants that occur fewer than threshold
                        times are ignored. Default 20
  --min_cover MIN_COVER
                        the minimum coverage at a site required to call a
                        variant. This parameter only has effect if it is >
                        threshold. Default 0
  --min_reads MIN_READS
                        the minimum number of variant containing reads
                        required to call a variant. This parameter only has
                        effect if it is > threshold. Default 0
  --min_freq MIN_FREQ   the minimum proportion of mutated reads. Default 0.1
                        (ten percent)
  --overlap OVERLAP     The minimum overlap required between forward and
                        reverse sequences to merge. Default 20bp
  --hashtable HASHTABLE
                        Filename for a precomputed hash table of exact matches
                        of amplicons to references. Generate with
                        --savehashtable
  --savehashtable SAVEHASHTABLE
                        Output a precomputed hash table that matches amplicons
                        exactly to references. Used to speed up matching with
                        --hashtable
  --alignments ALIGNMENTS
                        Print a formatted text version of variant containing
                        alignments to a file. "-" will print to stderr
  --prefix PREFIX       Shorthand specification of --forward, --reverse,
                        --countfile and --outfile. --forward =
                        <prefix>_R1.fastq.gz --reverse = <prefix>_R2.fastq.gz
                        --outfile = <prefix>.vcf --countfile = <prefix>.counts
						
						
Additional tools for simulating amplicon data
---------------------------------------------

simulate_variants
------------------

usage: simulate_variants [-h] [--manifest MANIFEST] [--fasta FASTA]
                         [--read1 READ1] [--read2 READ2] [--skip_softmasked]
                         [--position_excludes_softmasked]
                         [--deletions DELETIONS] [--check_sam CHECK_SAM]

A script for simulating paired end reads with variants from an amplicon target
file

optional arguments:
  -h, --help            show this help message and exit
  --manifest MANIFEST   an Illumina TrueSeq Amplicon manifest file. Default:
                        None
  --fasta FASTA         a fasta file of amplicons all in the plus strand
                        orientation with description lines ">name chromosome
                        start end" Default: None
  --read1 READ1         a fastq output file of forward reads. Default: stdout
  --read2 READ2         a fastq output file of reverse reads. Default: stdout
  --skip_softmasked     dont generate variants in softmasked sequence.
                        Default: True
  --position_excludes_softmasked
                        Exclude softmasked sequence when calculating start
                        site of read, cigar and variant detection strings.
                        Default: True
  --deletions DELETIONS
                        The size deletions to insert instead of point
                        variants. Default: None
  --check_sam CHECK_SAM
                        a sam file for parsing to identify entries where name
                        does not match the sam file mapping location and
                        variant strings

truseq_manifest
---------------

usage: truseq_manifest [-h] [--manifest MANIFEST] [--output OUTPUT] [--probes]
                       [--adaptors] [--with_probes] [--softmask_probes]
                       [--all_plus]

A script for converting Illumina TruSeq Amplicon manifest files to fasta files
Produces either a fasta file of target sequences without primers or a file of
primer sequences suitable for use by a trimming program (eg Nesoni clip)

optional arguments:
  -h, --help           show this help message and exit
  --manifest MANIFEST  an Illumina TruSeq Amplicon manifest file. Default:
                       stdin
  --output OUTPUT      a multi fasta output file of sequence targets. Default:
                       stdout
  --probes             output only the ULSO and DLSO primer sequences
  --adaptors           append Illumina adaptor sequences to the primer
                       sequences
  --with_probes        append the ULSO and DLSO sequences to the fasta target
                       sequences
  --softmask_probes    append the ULSO and DLSO sequences to the fasta target
                       sequences
  --all_plus           reorient target sequences so they are all presented on
                       the plus strand

