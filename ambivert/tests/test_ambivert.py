#!/usr/bin/env python3
# encoding: utf-8
"""
test_ambivert.py

Created by Matthew Wakefield on 2013-11-04.
Copyright (c) 2013 Matthew Wakefield and The University of Melbourne. All rights reserved.
"""

#from __future__ import division, print_function, unicode_literals
import unittest, io
import hashlib
from ambivert import ambivert

from pkg_resources import resource_stream

def md5(data):
    #python3 compatibility
    return hashlib.md5(bytes(data,'ascii'))

def make_testdata_fastq(amplicons): #pragma no cover
    """This function can be called in ambivert.main to make testdata files restricted to brca1 exon9"""
    keys = amplicons.get_amplicons_overlapping('chr17',41243125,41247176)
    forwardfile = open('testdata_R1.fastq','w')
    reversefile = open('testdata_R2.fastq','w')
    for key in keys:
        if amplicons.reference[key] in [('BRCA1_Exon9_UserDefined_(9825051)_7473609_chr17_41243125_41243349', 'chr17', '41243125', '41243349', '+'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473610_chr17_41243267_41243491', 'chr17', '41243267', '41243491', '-'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473611_chr17_41243405_41243630', 'chr17', '41243405', '41243630', '+'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473612_chr17_41243559_41243783', 'chr17', '41243559', '41243783', '-'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473613_chr17_41243701_41243927', 'chr17', '41243701', '41243927', '+'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473614_chr17_41243841_41244065', 'chr17', '41243841', '41244065', '-'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473615_chr17_41243981_41244206', 'chr17', '41243981', '41244206', '+'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473616_chr17_41244123_41244347', 'chr17', '41244123', '41244347', '-'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473617_chr17_41244261_41244486', 'chr17', '41244261', '41244486', '+'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473618_chr17_41244399_41244625', 'chr17', '41244399', '41244625', '-'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473619_chr17_41244539_41244764', 'chr17', '41244539', '41244764', '+'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473620_chr17_41244679_41244909', 'chr17', '41244679', '41244909', '-'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473621_chr17_41244855_41245081', 'chr17', '41244855', '41245081', '+'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473622_chr17_41245027_41245253', 'chr17', '41245027', '41245253', '-'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473623_chr17_41245195_41245420', 'chr17', '41245195', '41245420', '+'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473624_chr17_41245363_41245588', 'chr17', '41245363', '41245588', '-'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473625_chr17_41245533_41245757', 'chr17', '41245533', '41245757', '+'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473626_chr17_41245703_41245928', 'chr17', '41245703', '41245928', '-'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473627_chr17_41245865_41246105', 'chr17', '41245865', '41246105', '+'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473628_chr17_41246051_41246277', 'chr17', '41246051', '41246277', '-'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473629_chr17_41246219_41246443', 'chr17', '41246219', '41246443', '+'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473630_chr17_41246387_41246623', 'chr17', '41246387', '41246623', '-'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473631_chr17_41246571_41246801', 'chr17', '41246571', '41246801', '+'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473632_chr17_41246745_41246989', 'chr17', '41246745', '41246989', '-'),
        ('BRCA1_Exon9_UserDefined_(9825051)_7473633_chr17_41246933_41247176', 'chr17', '41246933', '41247176', '+'),
        ]:
            amplicons.print_to_fastq(key, forwardfile, reversefile)
    pass


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
        
        self.assertEqual(md5(str(self.amplicons.data[f_seq+r_seq])).hexdigest(),'d751713988987e9331980363e24189ce')
        pass
        
    
    def test_AmpliconData_process_twofile_readpairs(self):
        forward_file = resource_stream(__name__, 'data/testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'data/testdata_R2.fastq')
        self.amplicons.process_twofile_readpairs(forward_file, reverse_file)
        
        #print(str(self.amplicons.data))
        self.assertEqual(md5(str(sorted(self.amplicons.data))).hexdigest(),'0b8d563beaeb05e6b6ea8615fbf8906b')
        pass

    def test_AmpliconData_get_above_threshold(self):
        forward_file = resource_stream(__name__, 'data/testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'data/testdata_R2.fastq')
        self.amplicons.process_twofile_readpairs(forward_file, reverse_file)
        
        #self.assertEqual(md5(str(self.amplicons.data)).hexdigest(),'0f852572f90c142103b89eb4961720c4')
        self.assertEqual(md5(str(sorted(self.amplicons.get_above_threshold(1)))).hexdigest(),'0b8d563beaeb05e6b6ea8615fbf8906b')
        self.assertEqual(md5(str(sorted(self.amplicons.get_above_threshold(75)))).hexdigest(),'6d062386c9f2c16dc7a245a38fbcf60f')
        self.assertEqual(self.amplicons.get_above_threshold(175),['6bb3477f87edfbf67d6ab286926bff99'])
        
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
        self.assertEqual(self.amplicons.merged,{md5('TTACCTTCCATGAGTTGTAGGTTTCTGCTGTGCCTGACTGGCATTTGGTTGTACTTTTTTTTCTTTATCTCTTCACTGCTAGAACAACTATCAATTTGCAATTCAGTACAATTAGGTGGGCTTAGATTTCTACTGACTACTAGTTCAAGCGGTTAAATATCCACAATTCAAAAGCACCTAAAAAGAATAGGCTGAGGAGGAAGTCTTCTACCAGGCATATTCATGCGCTTGAACTAGTAGTCAGTAGAAATCTAAGCCCACCTAATTGTACTGAATTGCAAATTGATAGTTGTTCTAGCAGT').hexdigest(): 'TTACCTTCCATGAGTTGTAGGTTTCTGCTGTGCCTGACTGGCATTTGGTTGTACTTTTTTTTCTTTATCTCTTCACTGCTAGAACAACTATCAATTTGCAATTCAGTACAATTAGGTGGGCTTAGATTTCTACTGACTACTAGTTCAAGCGCATGAATATGCCTGGTAGAAGACTTCCTCCTCAGCCTATTCTTTTTAGGTGCTTTTGAATTGTGGATATTTAAC'})
        pass
    
    def test_process_amplicon_data(self):
        forward_file = resource_stream(__name__, 'data/testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'data/testdata_R2.fastq')
        manifest = resource_stream(__name__, 'testdatamanifest.txt')
        amplicons = ambivert.process_amplicon_data(forward_file, reverse_file,
                                  manifest=manifest, fasta=None,
                                  threshold=50, overlap=20, primer=15, 
                                  savehashtable=None, hashtable=None,
                                  )
        self.assertEqual(md5(str(amplicons.potential_variants)).hexdigest(),'f0a10362afd9d576c6082af8030f9cdf')
        pass
    
    def test_process_amplicon_data(self):
        forward_file = resource_stream(__name__, 'data/testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'data/testdata_R2.fastq')
        manifest = resource_stream(__name__, 'data/testdatamanifest.txt')
        amplicons = ambivert.process_amplicon_data(forward_file, reverse_file,
                                  manifest=manifest, fasta=None,
                                  threshold=50, overlap=20, primer=15, 
                                  savehashtable=None, hashtable=None,
                                  )
        self.assertEqual(md5(str(sorted(amplicons.potential_variants))).hexdigest(),'1a28f10a7a1e2ea430e79453e367a342')
    
    def test_AmpliconData_test_get_amplicon_count(self):
        forward_file = resource_stream(__name__, 'data/testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'data/testdata_R2.fastq')
        manifest = resource_stream(__name__, 'data/testdatamanifest.txt')
        amplicons = ambivert.process_amplicon_data(forward_file, reverse_file,
                                  manifest=manifest, fasta=None,
                                  threshold=50, overlap=20, primer=15, 
                                  savehashtable=None, hashtable=None,
                                  )
        #print(amplicons.get_amplicon_counts())
        self.assertEqual(md5(str(sorted(amplicons.get_amplicon_counts()))).hexdigest(),'564424f9a9e909cb9323d5ceac93a559')
        


    
if __name__ == '__main__':
    unittest.main()