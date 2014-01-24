#!/usr/bin/env python3
# encoding: utf-8
"""
test_ambivert.py

Created by Matthew Wakefield on 2013-11-04.
Copyright (c) 2013 Matthew Wakefield and The University of Melbourne. All rights reserved.

Incomplete test coverage of some conditional code.  Excludes command line interface code.  Current status:
coverage3 report --show-missing
Name                           Stmts   Miss  Cover   Missing
---------------------------------------------------------------------------------------------
ambivert/ambivert                305     21    93%   78-79, 120, 123, 202, 223-224, 232, 239, 262, 284, 323, 348-350, 417, 435-437, 552, 554, 558
"""

#from __future__ import division, print_function, unicode_literals
import unittest, io, os
import hashlib
from tempfile import NamedTemporaryFile
from pkg_resources import resource_stream
from ambivert import ambivert

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
        self.assertEqual(md5(str(sorted(self.amplicons.data))).hexdigest(),'37e10d52aa874a8fbdb40330e1294f02')
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
        
    def test_smithwaterman_align_DNA(self):
        s1 = 'GCGCCAGCCACTCAACTATATTTTGGTTAATATCTCCTTGGCTGGTTTCATCTACTGCATCTTCTCGGTCTTCACTGTCTTCATCTCTAGCTCCCAAGGC'
        s2 = 'GCTGGTTTCATCTACTGCATCTTCTCGGTCTTCACTGTCTTCATCTCTAGCTCCCAAGGCTACTTCATCTTTGGCCGCCATGTTTGTGCTATGGAGGGTT'
    
        aligned_s1, aligned_s2, start_s1, start_s2 = ambivert.smith_waterman(s1,s2)
        self.assertEqual(aligned_s1,aligned_s2)
        self.assertEqual(aligned_s1,s1[40:])
        self.assertEqual(start_s1, 40)
        self.assertEqual(start_s2, 0)

        s1 = 'GCGCCAGCCACTCAACTATATTTTGGTTAATATCTCCTTCACACACACACACTTCTCGGTCGTCGTCGTCTTCACTGTCTTCATCTCTAGCTCCCAAGGC'
        s2 = 'GCGCCAGCCACTCAACTATATTTGGTTAATATCTCCTTCACACACACACTTCTCGGTCGTCGTCTTCACTGTCTTCATCTCTAGCTCCCAAGGC'

        aligned_s1, aligned_s2, start_s1, start_s2 = ambivert.smith_waterman(s1,s2)
        self.assertEqual('GCGCCAGCCACTCAACTATA-TTTGGTTAATATCTCCTT--CACACACACACTTCTCG---GTCGTCGTCTTCACTGTCTTCATCTCTAGCTCCCAAGGC',aligned_s2)

        s1 = 'GCGCCAGCCACTCAACTATATTTTGGTTAATATCTCCTTCACACACACACACTTCTCGGTCGTCGTCGTCTTCACTGTCTTCATCTCTAGCTCCCAAGGC'
        s2 = 'GCGCCAGCCACTCAACTATATTTGGTTAATATCTCCTGCACACTTCTCGGTCGTCGTCTTCACTGTCTTCATCTCTAGCTCCCAAGGC'
        aligned_s1, aligned_s2, start_s1, start_s2 = ambivert.smith_waterman(s1,s2)
        
        self.assertEqual('GCGCCAGCCACTCAACTATA-TTTGGTTAATATCTCCT--------GCACACTTCTCG---GTCGTCGTCTTCACTGTCTTCATCTCTAGCTCCCAAGGC',aligned_s2)

    def test_AmpliconData_str(self):
        self.assertEqual(str(self.amplicons),'Data file names: \nReference file names: \nReadpairs: 0\nThreshold: 0\nUnique read groups: 0\nUnmerged pairs: 0\n')

    def test_AmpliconData_merge_overlaps(self):
        
        f_name = '@M00267:63:000000000-A5GL7:1:1104:11288:4465 1:N:0:15'
        f_seq = 'TTACCTTCCATGAGTTGTAGGTTTCTGCTGTGCCTGACTGGCATTTGGTTGTACTTTTTTTTCTTTATCTCTTCACTGCTAGAACAACTATCAATTTGCAATTCAGTACAATTAGGTGGGCTTAGATTTCTACTGACTACTAGTTCAAGCG'
        f_qual = 'BBBBBFFFFFFFGG5GGGGGGFHHHHFHHGHHHHHGFHHHHHGHHHGHGGFEFHHHHGHGGGHHHEFGHFDHHGFHHHHHHHHGHFHFHHEFFFHHHHGHHFHFHHHFFEGHHFHGBFFCGHHHFHFFHHHHHHHHHGHHHHHHHHEDHHD'
        r_name = '@M00267:63:000000000-A5GL7:1:1104:11288:4465 2:N:0:15'
        r_seq = 'GTTAAATATCCACAATTCAAAAGCACCTAAAAAGAATAGGCTGAGGAGGAAGTCTTCTACCAGGCATATTCATGCGCTTGAACTAGTAGTCAGTAGAAATCTAAGCCCACCTAATTGTACTGAATTGCAAATTGATAGTTGTTCTAGCAGT'
        r_qual = 'AABAAFFFBFFFGFGGGGGGGGHHHHHHHHGHHGGHGHGGHFHGGEHGGGHFHGHHHHHFFHHEHHFHHHHHHFHGGGGGHHHHHHGGHHHGFHHGHGGHHGGHHGHGGFGHHHGHGHHHHHHGHGHFHHHHHGFHHHGHHHHHGFFHHHB'
        self.amplicons.add_reads(f_name, f_seq, f_qual,r_name, r_seq, r_qual)
        
        f_name = 'dodgybrothers'
        f_seq =  'TTACCTTCCATGAGTTGTAGGTTTCTGCTG'
        f_qual = 'BBBBBFFFFFFFGG5GGGGGGFHHHHFHHG'
        r_name = 'dodgybrothers'                
        r_seq =  'GTTAAATATCCACAATTCAAAAGCACCTAA'
        r_qual = 'AABAAFFFBFFFGFGGGGGGGGHHHHHHHH'
        self.amplicons.add_reads(f_name, f_seq, f_qual,r_name, r_seq, r_qual)
        
        logfile = io.StringIO()
        self.amplicons.merge_overlaps()
        self.assertEqual(self.amplicons.merged,{md5('TTACCTTCCATGAGTTGTAGGTTTCTGCTGTGCCTGACTGGCATTTGGTTGTACTTTTTTTTCTTTATCTCTTCACTGCTAGAACAACTATCAATTTGCAATTCAGTACAATTAGGTGGGCTTAGATTTCTACTGACTACTAGTTCAAGCGGTTAAATATCCACAATTCAAAAGCACCTAAAAAGAATAGGCTGAGGAGGAAGTCTTCTACCAGGCATATTCATGCGCTTGAACTAGTAGTCAGTAGAAATCTAAGCCCACCTAATTGTACTGAATTGCAAATTGATAGTTGTTCTAGCAGT').hexdigest(): 'TTACCTTCCATGAGTTGTAGGTTTCTGCTGTGCCTGACTGGCATTTGGTTGTACTTTTTTTTCTTTATCTCTTCACTGCTAGAACAACTATCAATTTGCAATTCAGTACAATTAGGTGGGCTTAGATTTCTACTGACTACTAGTTCAAGCGCATGAATATGCCTGGTAGAAGACTTCCTCCTCAGCCTATTCTTTTTAGGTGCTTTTGAATTGTGGATATTTAAC'})
        self.assertEqual(self.amplicons.unmergable,[md5('TTACCTTCCATGAGTTGTAGGTTTCTGCTG'+'GTTAAATATCCACAATTCAAAAGCACCTAA').hexdigest(),])

    def test_add_references_from_manifest_and_fasta(self):
        manifest = io.TextIOWrapper(resource_stream(__name__, 'data/testdatamanifest.txt'))
        self.amplicons.add_references_from_manifest(manifest)
        self.assertEqual(md5(str(sorted(self.amplicons.reference_sequences.items()))).hexdigest(),'fc3c6701032dbd84c6f5731d344df060')
        fasta = io.TextIOWrapper(resource_stream(__name__, 'data/testdatareferences.fasta'))
        self.amplicons.reference_sequences = {}
        self.amplicons.add_references_from_fasta(fasta)
        self.assertEqual(md5(str(sorted(self.amplicons.reference_sequences.items()))).hexdigest(),'fc3c6701032dbd84c6f5731d344df060')
        pass
        
    def test_process_amplicon_data(self):
        forward_file = resource_stream(__name__, 'data/testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'data/testdata_R2.fastq')
        manifest = io.TextIOWrapper(resource_stream(__name__, 'data/testdatamanifest.txt'))
        amplicons = ambivert.process_amplicon_data(forward_file, reverse_file,
                                  manifest=manifest, fasta=None,
                                  threshold=50, overlap=20, 
                                  savehashtable=None, hashtable=None,
                                  )
        self.assertEqual(md5(str(amplicons.potential_variants)).hexdigest(),'f0a10362afd9d576c6082af8030f9cdf')
        pass
    
    def test_process_amplicon_data(self):
        forward_file = resource_stream(__name__, 'data/testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'data/testdata_R2.fastq')
        manifest = io.TextIOWrapper(resource_stream(__name__, 'data/testdatamanifest.txt'))
        amplicons = ambivert.process_amplicon_data(forward_file, reverse_file,
                                  manifest=manifest, fasta=None,
                                  threshold=50, overlap=20,
                                  savehashtable=None, hashtable=None,
                                  )
        self.assertEqual(md5(str(sorted(amplicons.potential_variants))).hexdigest(),'1a28f10a7a1e2ea430e79453e367a342')
    
    def test_AmpliconData_test_get_amplicon_counts(self):
        forward_file = resource_stream(__name__, 'data/testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'data/testdata_R2.fastq')
        manifest = io.TextIOWrapper(resource_stream(__name__, 'data/testdatamanifest.txt'))
        amplicons = ambivert.process_amplicon_data(forward_file, reverse_file,
                                  manifest=manifest, fasta=None,
                                  threshold=50, overlap=20, 
                                  savehashtable=None, hashtable=None,
                                  )
        self.assertEqual(md5(str(sorted(amplicons.get_amplicon_counts()))).hexdigest(),'564424f9a9e909cb9323d5ceac93a559')
        
        self.assertEqual(amplicons.get_amplicon_count('335eb2f760b98b4f25d281537400baf4'),144)
    
    def test_save_hash_table(self):
        forward_file = resource_stream(__name__, 'data/testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'data/testdata_R2.fastq')
        manifest = io.TextIOWrapper(resource_stream(__name__, 'data/testdatamanifest.txt'))
        amplicons = ambivert.process_amplicon_data(forward_file, reverse_file,
                                  manifest=manifest, fasta=None,
                                  threshold=50, overlap=20, 
                                  savehashtable=None, hashtable=None,
                                  )
        outfile = NamedTemporaryFile(delete=False,mode='wb')
        outfilename = outfile.name
        amplicons.save_hash_table(outfile)
        amplicons.reference = {}
        amplicons.load_hash_table(open(outfilename,mode='rb'))
        self.assertEqual(md5(str(sorted(amplicons.reference))).hexdigest(),'c5e174ffd7ee8cfbd47e162667c83796')
        outfile.close()
        os.unlink(outfile.name)
        pass
    
    def test_print_variants_as_alignments(self):
        forward_file = resource_stream(__name__, 'data/testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'data/testdata_R2.fastq')
        manifest = io.TextIOWrapper(resource_stream(__name__, 'data/testdatamanifest.txt'))
        amplicons = ambivert.process_amplicon_data(forward_file, reverse_file,
                                  manifest=manifest, fasta=None,
                                  threshold=50, overlap=20, 
                                  savehashtable=None, hashtable=None,
                                  )
        outfile = io.StringIO()
        amplicons.print_variants_as_alignments(outfile)
        #print(outfile.getvalue())
        self.assertEqual(md5(outfile.getvalue()).hexdigest(),'6905671f81cc776afecc8a13bc01b1c3')
        pass
    
    def test_call_amplicon_variants(self):
        forward_file = resource_stream(__name__, 'data/testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'data/testdata_R2.fastq')
        manifest = io.TextIOWrapper(resource_stream(__name__, 'data/testdatamanifest.txt'))
        amplicons = ambivert.process_amplicon_data(forward_file, reverse_file,
                                  manifest=manifest, fasta=None,
                                  threshold=50, overlap=20, 
                                  savehashtable=None, hashtable=None,
                                  )
        amplicons.call_amplicon_variants()
        #print(sorted(amplicons.called_variants.items()))
        self.assertEqual(md5(str(sorted(amplicons.called_variants.items()))).hexdigest(),'395ea8c7d8c55a4a752d9fa47dae8c8a')
        pass
    
    def test_get_variant_positions(self):
        forward_file = resource_stream(__name__, 'data/testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'data/testdata_R2.fastq')
        manifest = io.TextIOWrapper(resource_stream(__name__, 'data/testdatamanifest.txt'))
        amplicons = ambivert.process_amplicon_data(forward_file, reverse_file,
                                  manifest=manifest, fasta=None,
                                  threshold=50, overlap=20, 
                                  savehashtable=None, hashtable=None,
                                  )
        amplicons.call_amplicon_variants()
        #print(amplicons.get_variant_positions())
        self.assertEqual(md5(str(amplicons.get_variant_positions())).hexdigest(),'563614d673fc502e6a191b7d37dd2070')
        pass
    
    def test_get_amplicons_overlapping(self):
        forward_file = resource_stream(__name__, 'data/testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'data/testdata_R2.fastq')
        manifest = io.TextIOWrapper(resource_stream(__name__, 'data/testdatamanifest.txt'))
        amplicons = ambivert.process_amplicon_data(forward_file, reverse_file,
                                  manifest=manifest, fasta=None,
                                  threshold=50, overlap=20, 
                                  savehashtable=None, hashtable=None,
                                  )
        
        self.assertEqual(sorted(amplicons.get_amplicons_overlapping('chr17', '41245364')),['32ff3ea55601305b5e8b3bd266c4080a', '357199abf30c529a625b24916a341ff8', '4a2c50f29c387e9f5853d7e071c1d4ff'])
        self.assertEqual(sorted(amplicons.get_amplicons_overlapping_without_softmasking('chr17', '41245364')),['4a2c50f29c387e9f5853d7e071c1d4ff'])
        
        pass
    
    def test_consolidate_variants(self):
        forward_file = resource_stream(__name__, 'data/testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'data/testdata_R2.fastq')
        manifest = io.TextIOWrapper(resource_stream(__name__, 'data/testdatamanifest.txt'))
        amplicons = ambivert.process_amplicon_data(forward_file, reverse_file,
                                  manifest=manifest, fasta=None,
                                  threshold=50, overlap=20, 
                                  savehashtable=None, hashtable=None,
                                  )
        amplicons.call_amplicon_variants()
        amplicons.consolidate_variants()
        #print(sorted(amplicons.consolidated_variants))
        self.assertEqual(md5(str(sorted(amplicons.consolidated_variants))).hexdigest(),'a38f22f05f6c7d91f6615c404b206b2a')
        pass
    
    def test_get_filtered_variants(self):
        forward_file = resource_stream(__name__, 'data/testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'data/testdata_R2.fastq')
        manifest = io.TextIOWrapper(resource_stream(__name__, 'data/testdatamanifest.txt'))
        amplicons = ambivert.process_amplicon_data(forward_file, reverse_file,
                                  manifest=manifest, fasta=None,
                                  threshold=50, overlap=20, 
                                  savehashtable=None, hashtable=None,
                                  )
        amplicons.call_amplicon_variants()
        amplicons.consolidate_variants()
        #print(list(amplicons.get_filtered_variants(min_freq=0.5)))
        self.assertEqual(md5(str(list(amplicons.get_filtered_variants(min_freq=0.5)))).hexdigest(),'870ffb1202aa1f3acde75c18a6b4ff1f')
        pass

    def test_get_filtered_variants_filtering_ambiguous(self):
        forward_file = resource_stream(__name__, 'data/testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'data/testdata_R2.fastq')
        manifest = io.TextIOWrapper(resource_stream(__name__, 'data/testdatamanifest.txt'))
        amplicons = ambivert.process_amplicon_data(forward_file, reverse_file,
                                  manifest=manifest, fasta=None,
                                  threshold=0, overlap=10, 
                                  savehashtable=None, hashtable=None,
                                  )
        amplicons.call_amplicon_variants()
        amplicons.consolidate_variants()
        #print(list(amplicons.get_filtered_variants(min_reads=0, min_cover=0, min_freq=0.00001)))
        self.assertEqual(md5(str(list(amplicons.get_filtered_variants(min_reads=0, min_cover=0, min_freq=0.00001)))).hexdigest(),'4c5d87513811af17fe80bb7029490edf')
        pass

    
    def test_print_to_fastq(self):
        forward_file = resource_stream(__name__, 'data/testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'data/testdata_R2.fastq')
        manifest = io.TextIOWrapper(resource_stream(__name__, 'data/testdatamanifest.txt'))
        amplicons = ambivert.process_amplicon_data(forward_file, reverse_file,
                                  manifest=manifest, fasta=None,
                                  threshold=50, overlap=20, 
                                  savehashtable=None, hashtable=None,
                                  )
        outfile = io.StringIO()
        amplicons.print_to_fastq('335eb2f760b98b4f25d281537400baf4',forwardfile=outfile,reversefile=outfile)
        #print(outfile.getvalue())
        self.assertEqual(md5(outfile.getvalue()).hexdigest(),'f7208f5b8a2e6a33d48e76f861ed35c9')
        pass

    def test_is_homopolymer_at_position(self):
        forward_file = resource_stream(__name__, 'data/testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'data/testdata_R2.fastq')
        manifest = io.TextIOWrapper(resource_stream(__name__, 'data/testdatamanifest.txt'))
        amplicons = ambivert.process_amplicon_data(forward_file, reverse_file,
                                  manifest=manifest, fasta=None,
                                  threshold=50, overlap=20, 
                                  savehashtable=None, hashtable=None,
                                  )
        self.assertEqual(sorted(amplicons.get_reference_overlapping('chr17',41243290)),
                        [('BRCA1_Exon9_UserDefined_(9825051)_7473609_chr17_41243125_41243349', 'chr17', '41243125', '41243349', '+'),
                        ('BRCA1_Exon9_UserDefined_(9825051)_7473610_chr17_41243267_41243491', 'chr17', '41243267', '41243491', '-')])
        
        self.assertTrue(amplicons.is_homopolymer_at_position('chr17',41243290))
        self.assertTrue(amplicons.is_homopolymer_at_position('chr17',41243290, minimum=6))
        self.assertFalse(amplicons.is_homopolymer_at_position('chr17',41243290, minimum=7))
        pass
        
        
    
    def test_print_consolidated_vcf(self):
        forward_file = resource_stream(__name__, 'data/testdata_R1.fastq')
        reverse_file = resource_stream(__name__, 'data/testdata_R2.fastq')
        manifest = io.TextIOWrapper(resource_stream(__name__, 'data/testdatamanifest.txt'))
        amplicons = ambivert.process_amplicon_data(forward_file, reverse_file,
                                  manifest=manifest, fasta=None,
                                  threshold=50, overlap=20, 
                                  savehashtable=None, hashtable=None,
                                  )
        outfile = io.StringIO()
        amplicons.print_consolidated_vcf(outfile=outfile)
        #print(outfile.getvalue())
        self.assertEqual(md5(outfile.getvalue()).hexdigest(),'8c132f329a65a01fc7c31c456fe7f7c9')
        pass
    
    
if __name__ == '__main__':
    unittest.main()