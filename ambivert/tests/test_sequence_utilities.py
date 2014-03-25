#!/usr/bin/env python3
# encoding: utf-8
"""
test_truseq_manifest.py

Created by Matthew Wakefield.
Copyright (c) 2013-2014  Matthew Wakefield and The University of Melbourne. All rights reserved.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

"""

import unittest
from ambivert.sequence_utilities import *

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2013-2014,  Matthew Wakefield and The University of Melbourne"
__credits__ = ["Matthew Wakefield",]
__license__ = "GPLv3"
__version__ = "0.1.11"
__maintainer__ = "Matthew Wakefield"
__email__ = "matthew.wakefield@unimelb.edu.au"
__status__ = "Development"

def md5(data):
    #python3 compatibility
    return hashlib.md5(bytes(data,'ascii'))

class test_sequence_utilities(unittest.TestCase):
    def setUp(self):
        pass
    def test_fix_softmasked_expanded_cigar(self):
        self.assertEqual(fix_softmasked_expanded_cigar(['S', 'S', 'S', 'S', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'S', 'S', 'S', 'S']),
                        ('SSSSMMMMMMMMSSSS', 4, 8))
        self.assertEqual(fix_softmasked_expanded_cigar(['S', 'S', 'D', 'S', 'M', 'M', 'M', 'I', 'M', 'D', 'M', 'M', 'S', 'I', 'S', 'S']),
                        ('SSSMMMIMDMMSSSS', 3, 8))
        self.assertEqual(fix_softmasked_expanded_cigar(['S', 'S', 'D', 'S', 'M', 'M', 'M', 'I', 'M', 'D', 'M', 'M']),
                        ('SSSMMMIMDMM', 3, 8))
        self.assertEqual(fix_softmasked_expanded_cigar(['M', 'M', 'M', 'I', 'M', 'D', 'M', 'M', 'S', 'I', 'S', 'S']),
                        ('MMMIMDMMSSSS', 0, 8))
        pass
        
    def test_gapped_alignment_to_cigar(self):
        self.assertEqual(gapped_alignment_to_cigar('gtacACGTACGTgtac','GTACACGTACGTGTAC'),
                        ('4S8M4S', 4, 8)
                        )
        self.assertEqual(gapped_alignment_to_cigar('gtacACG-ACGTg-ac','GT-CACGTA-GTGTAC'),
                        ('3S3M1I1M1D2M4S', 3, 8)
                        )
        pass
    
    def test_expand_cigar(self):
        self.assertEqual(expand_cigar('80M5D2M2I10M'),'M'*80+'D'*5+'M'*2+'I'*2+'M'*10)
        self.assertEqual(expand_cigar('2S4M5D2M2I10M'),'SSMMMMDDDDDMMIIMMMMMMMMMM')
        pass
    
    def test_compact_cigar(self):
        self.assertEqual('80M5D2M2I10M',compact_cigar('M'*80+'D'*5+'M'*2+'I'*2+'M'*10))
        self.assertEqual('2S4M5D2M2I10M',compact_cigar('SSMMMMDDDDDMMIIMMMMMMMMMM'))
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
    
    def test_expand_mdtag(self):
        self.assertEqual(['', '', '', '', '', '', '', 'A', '', '', '', '', '', '', '', ''],expand_mdtag('7A8'))
        self.assertEqual(['', '', 'G', '', '', 'A', '', ''],expand_mdtag('2G2A2'))
        self.assertEqual(['G', '', '', 'A'],expand_mdtag('G2A'))
        self.assertEqual(['', '', '', '', '', '', '', '^' ,'c','a','t', '', '', '', '', '', '', '', ''],expand_mdtag('7^CAT8'))
        self.assertEqual(['', '', '', '', '', '', '',],expand_mdtag('7'))
        self.assertEqual(['', '', '', '', '', '', '', '^' ,'c','a','t', 'G', '', '', '', '', '', '', '', ''],expand_mdtag('7^CAT0G8'))
        pass
    
    def test_compact_expanded_mdtag_tokens(self):
        self.assertEqual(compact_expanded_mdtag_tokens(['', '', '', '', '', '', '', 'A', '', '', '', '', '', '', '', '']),'7A8')
        self.assertEqual(compact_expanded_mdtag_tokens(['', '', 'G', '', '', 'A', '', '']),'2G2A2')
        self.assertEqual(compact_expanded_mdtag_tokens(['G', '', '', 'A']),'G2A')
        self.assertEqual(compact_expanded_mdtag_tokens(['', '', '', '', '', '', '',]),'7')
        self.assertEqual(compact_expanded_mdtag_tokens(['', '', '', '', '', '', '', '^' ,'c','a','t', 'G', '', '', '', '', '', '', '']),'7^CAT0G7')
        self.assertEqual(compact_expanded_mdtag_tokens(['', '', '', '', '', '', '', '^' ,'c','a','t', '^', 'g', '', '', '', '', '', '', '']),'7^CATG7')
        pass
        
    
    
    
    
if __name__ == '__main__':
    unittest.main()