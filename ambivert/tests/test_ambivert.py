#!/usr/bin/env python2.7
# encoding: utf-8
"""
test_ambivert.py

Created by Matthew Wakefield on 2013-11-04.
Copyright (c) 2013 Matthew Wakefield and The University of Melbourne. All rights reserved.
"""

from __future__ import division, print_function
import unittest, io
from hashlib import md5
from ambivert import ambivert

from pkg_resources import resource_stream

class test_ambivert(unittest.TestCase):
    def setUp(self):
        self.amplicons = ambivert.AmpliconData()
        pass
    
    def test_AmpliconData_addreads(self):
        f_name = '@M00267:63:000000000-A5GL7:1:1104:11288:4465 1:N:0:15'
        f_seq = 'TTACCTTCCATGAGTTGTAGGTTTCTGCTGTGCCTGACTGGCATTTGGTTGTACTTTTTTTTCTTTATCTCTTCACTGCTAGAACAACTATCAATTTGCAATTCAGTACAATTAGGTGGGCTTAGATTTCTACTGACTACTAGTTCAAGCG'
        f_qual = 'BBBBBFFFFFFFGG5GGGGGGFHHHHFHHGHHHHHGFHHHHHGHHHGHGGFEFHHHHGHGGGHHHEFGHFDHHGFHHHHHHHHGHFHFHHEFFFHHHHGHHFHFHHHFFEGHHFHGBFFCGHHHFHFFHHHHHHHHHGHHHHHHHHEDHHD'
        r_name = '@M00267:63:000000000-A5GL7:1:1104:11288:4465 2:N:0:15'
        r_seq = 'GTTAAATATCCACAATTCAAAAGCACCTAAAAAGAATAGGCTGAGGAGGAAGTCTTCTACCAGGCATATTCATGCGCTTGAACTAGTAGTCAGTAGAAATCTAAGCCCACCTAATTGTACTGAATTGCAAATTGATAGTTGTTCTAGCAGT'
        r_qual = 'AABAAFFFBFFFGFGGGGGGGGHHHHHHHHGHHGGHGHGGHFHGGEHGGGHFHGHHHHHFFHHEHHFHHHHHHFHGGGGGHHHHHHGGHHHGFHHGHGGHHGGHHGHGGFGHHHGHGHHHHHHGHGHFHHHHHGFHHHGHHHHHGFFHHHB'
        self.amplicons.add_reads(f_name, f_seq, f_qual,r_name, r_seq, r_qual)
        
        self.assertEqual(md5(str(self.amplicons.data[f_seq+r_seq])).hexdigest(),'a9f89bfd3d14021dd43afaefbdeeb0f6')
        pass
        
    
    def test_AmpliconData_process_twofile_readpairs(self):
        forward_file = resource_stream(__name__, 'testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'testdata_R2.fastq')
        self.amplicons.process_twofile_readpairs(forward_file, reverse_file)
        
        #print(str(self.amplicons.data))
        self.assertEqual(md5(str(self.amplicons.data)).hexdigest(),'0f852572f90c142103b89eb4961720c4')
        pass

    def test_AmpliconData_get_above_threshold(self):
        forward_file = resource_stream(__name__, 'testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'testdata_R2.fastq')
        self.amplicons.process_twofile_readpairs(forward_file, reverse_file)
        
        #self.assertEqual(md5(str(self.amplicons.data)).hexdigest(),'0f852572f90c142103b89eb4961720c4')
        
        self.assertEqual(md5(str(self.amplicons.get_above_threshold(1))).hexdigest(),'e1fdc5020c89f23bcb3f498f9240487a')
        self.assertEqual(md5(str(self.amplicons.get_above_threshold(75))).hexdigest(),'c7b7ae699cbc0630ba336b8716b6617b')
        self.assertEqual(self.amplicons.get_above_threshold(175),['ACTTCTATAAATAGACTGGGGCAAACACAAAAACCTGGTTCCAATACCTAAGTTTGAATCCATGCTTTGCTCTTCTTGATTATTTTCTTCCAAGCCCGTTCCTCTTTCTTCATCATCTGAAACCAATTCCTTGTCACTCAGACCAACTCCCGACTGCAAATACAAACACCCAGGATCCTTTCTTGATTGGTTCTTCCAAACAAATGAGGCATCAGTCTGAAAGCCAGGGAGTTGGTCTGAGTGACAAGGAATTGGTTTCAGATGATGAAGAAAGAGGAACGGGCTTGGAAGAAAATAATCAA'])
        
        pass

    def test_AmpliconData_merge_overlaps(self):
        
        f_name = '@M00267:63:000000000-A5GL7:1:1104:11288:4465 1:N:0:15'
        f_seq = 'TTACCTTCCATGAGTTGTAGGTTTCTGCTGTGCCTGACTGGCATTTGGTTGTACTTTTTTTTCTTTATCTCTTCACTGCTAGAACAACTATCAATTTGCAATTCAGTACAATTAGGTGGGCTTAGATTTCTACTGACTACTAGTTCAAGCG'
        f_qual = 'BBBBBFFFFFFFGG5GGGGGGFHHHHFHHGHHHHHGFHHHHHGHHHGHGGFEFHHHHGHGGGHHHEFGHFDHHGFHHHHHHHHGHFHFHHEFFFHHHHGHHFHFHHHFFEGHHFHGBFFCGHHHFHFFHHHHHHHHHGHHHHHHHHEDHHD'
        r_name = '@M00267:63:000000000-A5GL7:1:1104:11288:4465 2:N:0:15'
        r_seq = 'GTTAAATATCCACAATTCAAAAGCACCTAAAAAGAATAGGCTGAGGAGGAAGTCTTCTACCAGGCATATTCATGCGCTTGAACTAGTAGTCAGTAGAAATCTAAGCCCACCTAATTGTACTGAATTGCAAATTGATAGTTGTTCTAGCAGT'
        r_qual = 'AABAAFFFBFFFGFGGGGGGGGHHHHHHHHGHHGGHGHGGHFHGGEHGGGHFHGHHHHHFFHHEHHFHHHHHHFHGGGGGHHHHHHGGHHHGFHHGHGGHHGGHHGHGGFGHHHGHGHHHHHHGHGHFHHHHHGFHHHGHHHHHGFFHHHB'
        self.amplicons.add_reads(f_name, f_seq, f_qual,r_name, r_seq, r_qual)
        
        logfile = io.StringIO()
        self.amplicons.merge_overlaps()
        self.assertEqual(self.amplicons.merged,{'TTACCTTCCATGAGTTGTAGGTTTCTGCTGTGCCTGACTGGCATTTGGTTGTACTTTTTTTTCTTTATCTCTTCACTGCTAGAACAACTATCAATTTGCAATTCAGTACAATTAGGTGGGCTTAGATTTCTACTGACTACTAGTTCAAGCGGTTAAATATCCACAATTCAAAAGCACCTAAAAAGAATAGGCTGAGGAGGAAGTCTTCTACCAGGCATATTCATGCGCTTGAACTAGTAGTCAGTAGAAATCTAAGCCCACCTAATTGTACTGAATTGCAAATTGATAGTTGTTCTAGCAGT': 'TTACCTTCCATGAGTTGTAGGTTTCTGCTGTGCCTGACTGGCATTTGGTTGTACTTTTTTTTCTTTATCTCTTCACTGCTAGAACAACTATCAATTTGCAATTCAGTACAATTAGGTGGGCTTAGATTTCTACTGACTACTAGTTCAAGCGCATGAATATGCCTGGTAGAAGACTTCCTCCTCAGCCTATTCTTTTTAGGTGCTTTTGAATTGTGGATATTTAAC'})
        pass
    
    def test_process_amplicon_data(self):
        forward_file = resource_stream(__name__, 'testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'testdata_R2.fastq')
        manifest = resource_stream(__name__, 'testdatamanifest.txt')
        amplicons = ambivert.process_amplicon_data(forward_file, reverse_file,
                                  manifest=manifest, fasta=None,
                                  threshold=50, overlap=20, primer=15, 
                                  savehashtable=None, hashtable=None,
                                  )
        self.assertEqual(md5(str(amplicons.potential_variants)).hexdigest(),'f0a10362afd9d576c6082af8030f9cdf')
        pass
    
    def test_process_amplicon_data(self):
        forward_file = resource_stream(__name__, 'testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'testdata_R2.fastq')
        manifest = resource_stream(__name__, 'testdatamanifest.txt')
        amplicons = ambivert.process_amplicon_data(forward_file, reverse_file,
                                  manifest=manifest, fasta=None,
                                  threshold=50, overlap=20, primer=15, 
                                  savehashtable=None, hashtable=None,
                                  )
        self.assertEqual(md5(str(amplicons.potential_variants)).hexdigest(),'f0a10362afd9d576c6082af8030f9cdf')
    
    def test_AmpliconData_test_get_amplicon_count(self):
        forward_file = resource_stream(__name__, 'testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'testdata_R2.fastq')
        manifest = resource_stream(__name__, 'testdatamanifest.txt')
        amplicons = ambivert.process_amplicon_data(forward_file, reverse_file,
                                  manifest=manifest, fasta=None,
                                  threshold=50, overlap=20, primer=15, 
                                  savehashtable=None, hashtable=None,
                                  )
        print(amplicons.get_amplicon_counts())
        self.assertEqual(md5(str(amplicons.get_amplicon_counts())).hexdigest(),'aa84de03ec1385010cca644d3ce744a0')
        


    
if __name__ == '__main__':
    unittest.main()