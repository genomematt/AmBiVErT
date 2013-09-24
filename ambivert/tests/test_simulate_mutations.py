#!/usr/bin/env python
# encoding: utf-8
"""
test_simulate_mutations.py

Created by Matthew Wakefield on 2013-09-24.
Copyright (c) 2013 Matthew Wakefield and The University of Melbourne. All rights reserved.
"""

import unittest
from ambivert.simulate_mutations import *
from hashlib import md5

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2013,  Matthew Wakefield and The University of Melbourne"
__credits__ = ["Matthew Wakefield",]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Matthew Wakefield"
__email__ = "matthew.wakefield@unimelb.edu.au"
__status__ = "Development"



class test_simulate_mutations(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_point_mutate_sequence(self):
        sequence = 'cagtGATCGATCacgt'
        self.assertEqual(list(point_mutate_sequence('cagtGATCGATCacgt',position_excludes_softmasked=False,one_based_site=5)),
                        [('None_1_16M_4A11', 'cagtAATCGATCacgt'), ('None_1_16M_4C11', 'cagtCATCGATCacgt'), ('None_1_16M_4T11', 'cagtTATCGATCacgt')])
        self.assertEqual(list(point_mutate_sequence('cagtGATCGATCacgt',chromosome='chrX', one_based_start=123456, position_excludes_softmasked=False,one_based_site=7)),
                        [('chrX_123456_16M_6A9', 'cagtGAACGATCacgt'), ('chrX_123456_16M_6C9', 'cagtGACCGATCacgt'), ('chrX_123456_16M_6G9', 'cagtGAGCGATCacgt')])
        self.assertEqual(list(point_mutate_sequence('cagtGATCGATCacgt',position_excludes_softmasked=True,one_based_site=5)),
                        [('None_5_8M_A7', 'cagtAATCGATCacgt'), ('None_5_8M_C7', 'cagtCATCGATCacgt'), ('None_5_8M_T7', 'cagtTATCGATCacgt')])
        self.assertEqual(list(point_mutate_sequence('cagtGATCGATCacgt',chromosome='chrX', one_based_start=123456, position_excludes_softmasked=True,one_based_site=7)),
                        [('chrX_123460_8M_2A5', 'cagtGAACGATCacgt'), ('chrX_123460_8M_2C5', 'cagtGACCGATCacgt'), ('chrX_123460_8M_2G5', 'cagtGAGCGATCacgt')])
        pass
    
    def test_make_all_point_mutations(self):
        sequence = 'cagtGATCGATCacgt'
        #print(list(make_all_point_mutations(sequence)))
        self.assertEqual('63d7c766aebd4de2d4f6ac93ca0c2550',md5(str(list(make_all_point_mutations(sequence)))).hexdigest())
        
        #print(list(make_all_point_mutations(sequence, skip_softmasked=False, position_excludes_softmasked=False)))
        self.assertEqual('e39f6ab49cf9fa517368a8dcba416c97',md5(str(list(make_all_point_mutations(sequence, skip_softmasked=False, position_excludes_softmasked=False)))).hexdigest())
        pass
    
    def test_format_fastq_entry(self):
        #name,sequence,q_generator=constant_quality
        self.assertEqual('@None_1_16M_4A11\ncagtAATCGATCacgt\n+\nIIIIIIIIIIIIIIII\n',
                        format_fastq_entry(*('None_1_16M_4A11', 'cagtAATCGATCacgt'))
                        )
        pass
    
    def test_expand_cigar(self):
        self.assertEqual(expand_cigar('80M5D2M2I10M'),'M'*80+'D'*5+'M'*2+'I'*2+'M'*10)
        self.assertEqual(expand_cigar('2S4M5D2M2I10M'),'SSMMMMDDDDDMMIIMMMMMMMMMM')
        pass
    
    def test_engap(self):
        self.assertEqual('MMMM-----MMIIMMMMMMMMMM',
                        engap(seq = 'MMMMMMIIMMMMMMMMMM', cigar='4M5D2M2I10M', delete='D', insert='I', match='M', gap='-'))
        self.assertEqual('MMMMDDDDDMM--MMMMMMMMMM',
                        engap(seq = 'MMMMDDDDDMMMMMMMMMMMM', cigar='4M5D2M2I10M', delete='I', insert='D', match='M', gap='-'))
        pass
    
    def test_cigar_trimmer(self):
        self.assertEqual('4M5D2M2I10M',cigar_trimmer(cigar='4M5D2M2I10M',trim_from_start=0,trim_from_end=0))
        self.assertEqual('1M2I10M',cigar_trimmer(cigar='4M5D2M2I10M',trim_from_start=5,trim_from_end=0))
        self.assertEqual('3M',cigar_trimmer(cigar='4M5D2M2I10M',trim_from_start=0,trim_from_end=15))
        self.assertEqual('2M',cigar_trimmer(cigar='4M5D2M2I10M',trim_from_start=1,trim_from_end=15))
        pass
    
    def test_compact_cigar(self):
        self.assertEqual('80M5D2M2I10M',compact_cigar('M'*80+'D'*5+'M'*2+'I'*2+'M'*10))
        self.assertEqual('2S4M5D2M2I10M',compact_cigar('SSMMMMDDDDDMMIIMMMMMMMMMM'))
        pass
    
    def test_expand_mdtag_tokens(self):
        self.assertEqual(['', '', '', '', '', '', '', 'A', '', '', '', '', '', '', '', ''],expand_mdtag_tokens('7A8'))
        self.assertEqual(['', '', 'G', '', '', 'A', '', ''],expand_mdtag_tokens('2G2A2'))
        self.assertEqual(['G', '', '', 'A'],expand_mdtag_tokens('G2A'))
        pass
        
    def test_compact_expanded_mdtag_tokens(self):
        self.assertEqual(compact_expanded_mdtag_tokens(['', '', '', '', '', '', '', 'A', '', '', '', '', '', '', '', '']),'7A8')
        self.assertEqual(compact_expanded_mdtag_tokens(['', '', 'G', '', '', 'A', '', '']),'2G2A2')
        self.assertEqual(compact_expanded_mdtag_tokens(['G', '', '', 'A']),'G2A')
        pass
        
    def test_mutation_detection_tag_trimmer(self):
        self.assertEqual(mutation_detection_tag_trimmer('7A8',trim_from_start=0,trim_from_end=0),'7A8')
        self.assertEqual(mutation_detection_tag_trimmer('7A8',trim_from_start=3,trim_from_end=0),'4A8')
        self.assertEqual(mutation_detection_tag_trimmer('7A8',trim_from_start=0,trim_from_end=5),'7A3')
        self.assertEqual(mutation_detection_tag_trimmer('7A8',trim_from_start=3,trim_from_end=5),'4A3')
        self.assertEqual(mutation_detection_tag_trimmer('7A8G',trim_from_start=3,trim_from_end=5),'4A4')
        pass
    
    def test_mutated_amplicon_to_paired_reads(self):
        self.assertEqual(mutated_amplicon_to_paired_reads('cagtGATCGATCacgt','chrX','12345','16M','7C8',readlength=10),
                        (('cagtGATCGA', '', 'chrX', '12345', '10M', '7C2'), ('acgtGATCGA', '', 'chrX', '12351', '10M', '1C8')))
        self.assertEqual(mutated_amplicon_to_paired_reads('cagtGATCGATCacgt','chrX','12345','16M','7C8'),
                        (('cagtGATCGATCacgt', '', 'chrX', '12345', '16M', '7C8'), ('acgtGATCGATCactg', '', 'chrX', '12345', '16M', '7C8')))
        self.assertEqual(mutated_amplicon_to_paired_reads('cagtGATCGATCacgt','chrX','12345','16M','7C8',quality='ABCGEFGHIJKLMNOP'),
                         (('cagtGATCGATCacgt', 'ABCGEFGHIJKLMNOP', 'chrX', '12345', '16M', '7C8'),
                          ('acgtGATCGATCactg', 'PONMLKJIHGFEGCBA', 'chrX', '12345', '16M', '7C8')))
        pass
    

if __name__ == '__main__':
    unittest.main()