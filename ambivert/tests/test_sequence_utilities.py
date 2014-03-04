#!/usr/bin/env python3
# encoding: utf-8
"""
test_truseq_manifest.py

Created by Matthew Wakefield on 2013-09-23.
Copyright (c) 2013 The University of Melbourne. All rights reserved.
"""

import unittest
from ambivert.sequence_utilities import *

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
    
if __name__ == '__main__':
    unittest.main()