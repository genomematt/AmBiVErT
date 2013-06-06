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
    4)  TODO: Consolidate calls
    5)  TODO: Output VCF format calls


Created by Matthew Wakefield and Graham Taylor.
Copyright (c) 2013  Matthew Wakefield and The University of Melbourne. All rights reserved.

usage: ambivert.py [-h] -f FORWARD -r REVERSE ( -m MANIFEST | --fasta FASTA)
                   [--output OUTPUT] [--countfile COUNTFILE]
                   [--threshold THRESHOLD] [--overlap OVERLAP]
                   [--primer PRIMER] [--hashtable HASHTABLE]
                   [--savehashtable SAVEHASHTABLE]

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
                        reads. May be compressed with gzip or bzip2 with the
                        appropriate suffix (.gz/.bz2)
  -r REVERSE, --reverse REVERSE
                        a fastq format file of reverse direction amplicon
                        reads. May be compressed with gzip or bzip2 with the
                        appropriate suffix (.gz/.bz2)
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
                        the minimum occurance threshold. Amplicons that occur
                        fewer than threshold times are ignored.
  --overlap OVERLAP     The minimum overlap required between forward and
                        reverse sequences to merge
  --primer PRIMER       The size of the smallest primer. This number of bases
                        is trimmed from the end of the merged sequences to
                        reduce the possibility that small amplicons will fail
                        to match due to primer mismatch
  --hashtable HASHTABLE
                        Filename for a precomputed hash table that matches
                        amplicons to references. Generate with --savehashtable
  --savehashtable SAVEHASHTABLE
                        Output a precomputed hash table that matches amplicons
                        to references. Use to speed up matching with
                        --hashtable

