#!/usr/bin/env python
# encoding: utf-8
"""
test_call_variants.py

Created by Matthew Wakefield and Graham Taylor.
Copyright (c) 2013-2015  Matthew Wakefield and The University of Melbourne. All rights reserved.

   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

"""

from __future__ import print_function
import unittest
from io import StringIO
from ambivert import call_variants

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2013-2015,  Matthew Wakefield and The University of Melbourne"
__credits__ = ["Matthew Wakefield","Graham Taylor"]
__license__ = "GPLv3"
__version__ = "0.5.1"
__maintainer__ = "Matthew Wakefield"
__email__ = "matthew.wakefield@unimelb.edu.au"
__status__ = "Development"

class test_call_variants(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_caller(self):
        """Caller takes aligned sequence and returns a list of tuples of:
            Muation type: 'I': Insertion, 'D': Deletion, 'X': Substitition
            Start: Zero based position in gapped sequence
            Pos:   One based position in ungapped reference sequence
            Bases: Substitute, inserted or deleted sequence
        """
        self.assertEqual(list(call_variants.caller('AAAA','AAGA')),[('X', 2, 3, 'A')])
        self.assertEqual(list(call_variants.caller('AAAA','AGGA')),[('X', 1, 2, 'A'), ('X', 2, 3, 'A')])
        self.assertEqual(list(call_variants.caller('AAAA','AA-A')),[('I', 2, 3, 'A')])
        self.assertEqual(list(call_variants.caller('AA-A','AAAA')),[('D', 2, 3, 'A')])
        self.assertEqual(list(call_variants.caller('AAAA','ACGA')),[('X', 1, 2, 'A'), ('X', 2, 3, 'A')])
        self.assertEqual(list(call_variants.caller('AATA','AC-A')),[('X', 1, 1, 'A'), ('I', 2, 3, 'T')])
        self.assertEqual(list(call_variants.caller('AA-A','ACAA')),[('X', 1, 2, 'A'), ('D', 2, 3, 'A')])
        self.assertEqual(list(call_variants.caller('AAAAA','AAGCA')),[('X', 2, 3, 'A'), ('X', 3, 4, 'A')])
        self.assertEqual(list(call_variants.caller('AAAAA','AA-CA')),[('I', 2, 3, 'A'), ('X', 3, 3, 'A')])
        self.assertEqual(list(call_variants.caller('AA-AA','AAACA')),[('D', 2, 3, 'A'), ('X', 3, 4, 'A')])
        self.assertEqual(list(call_variants.caller('GAAAA','gAAGA')),[('X', 3, 4, 'A')])
        self.assertEqual(list(call_variants.caller('GAAAA','gAA-A')),[('I', 3, 4, 'A')])
        self.assertEqual(list(call_variants.caller('GAA-A','gAAAA')),[('D', 3, 4, 'A')])
        self.assertEqual(list(call_variants.caller('CAA-A','gAAAA')),[('D', 3, 4, 'A')])
        self.assertEqual(list(call_variants.caller('CAA-A','gAAAA',softmask=False)),[('X', 0, 1, 'C'), ('D', 3, 4, 'A')])
        self.assertEqual(list(call_variants.caller('CAA-ATCT','gAAAAttt')),[('D', 3, 4, 'A')])
        self.assertEqual(list(call_variants.caller('CAA-ATCT','gAAAAttt',softmask=False)),[('X', 0, 1, 'C'), ('D', 3, 4, 'A'), ('X', 6, 7, 'C')])
        self.assertEqual(list(call_variants.caller('CAA-AT-T','gAAAAttt')),[('D', 3, 4, 'A')])
        self.assertEqual(list(call_variants.caller('CAA-AT-T','gAAAAttt',softmask=False)),[('X', 0, 1, 'C'), ('D', 3, 4, 'A'), ('D', 6, 7, 't')])
        self.assertEqual(list(call_variants.caller('CAA-ATTT','gAAAAt-t')),[('D', 3, 4, 'A')])
        self.assertEqual(list(call_variants.caller('CAA-ATTT','gAAAAt-t',softmask=False)),[('X', 0, 1, 'C'), ('D', 3, 4, 'A'), ('I', 6, 7, 'T')])
        
        #softmasked insertions
        self.assertEqual(list(call_variants.caller('CGGGGGCCATCAT','c-ggggcCTTCAT',softmask=False)),[('I', 1, 2, 'G'), ('X', 8, 8, 'A')])
        self.assertEqual(list(call_variants.caller('CGGGGGCCATCAT','c-ggggcCTTCAT',softmask=True)),[('X', 8, 8, 'A')])
        self.assertEqual(list(call_variants.caller('CGGGGGCCATCAT','cggggg-CTTCAT',softmask=False)),[('I', 6, 7, 'C'), ('X', 8, 8, 'A')])
        self.assertEqual(list(call_variants.caller('CGGGGGCCATCAT','cggggg-CTTCAT',softmask=True)),[('I', 6, 7, 'C'), ('X', 8, 8, 'A')])
        self.assertEqual(list(call_variants.caller('CGGGGGCCATCAT','c--gggcCTTCAT',softmask=False)),[('I', 1, 2, 'GG'), ('X', 8, 7, 'A')])
        self.assertEqual(list(call_variants.caller('CGGGGGCCATCAT','c--gggcCTTCAT',softmask=True)),[('X', 8, 7, 'A')])
        self.assertEqual(list(call_variants.caller('CGGGGGCCATCAT','cggggg--TTCAT',softmask=False)),[('I', 6, 7, 'CC'), ('X', 8, 7, 'A')])
        self.assertEqual(list(call_variants.caller('CGGGGGCCATCAT','cggggg--TTCAT',softmask=True)),[('I', 6, 7, 'CC'), ('X', 8, 7, 'A')])
        
        #adjacent insertion and deletions
        self.assertEqual(list(call_variants.caller('CGGGGGC-TCAT','CGGGGG-ATCAT',softmask=True)),[('I', 6, 7, 'C'), ('D', 7, 7, 'A')])
        self.assertEqual(list(call_variants.caller('CGGGGGCC-TCAT','CGGGGG--ATCAT',softmask=True)),[('I', 6, 7, 'CC'), ('D', 8, 7, 'A')])
        self.assertEqual(list(call_variants.caller('CGGGGG-ATCAT','CGGGGGC-TCAT',softmask=True)),[('D', 6, 7, 'C'), ('I', 7, 8, 'A')])
        self.assertEqual(list(call_variants.caller('CGGGGG--ATCAT','CGGGGGCC-TCAT',softmask=True)),[('D', 6, 7, 'CC'), ('I', 8, 9, 'A')])
        self.assertEqual(list(call_variants.caller('CGGGGG---ATCAT','CGGGGGCCC-TCAT',softmask=True)),[('D', 6, 7, 'CCC'), ('I', 9, 10, 'A')])
        
        #terminal mutation states
        self.assertEqual(list(call_variants.caller('CGGGGGCCATCAT','CGGGGGCCATCAG',softmask=True)),[('X', 12, 13, 'T')])
        self.assertEqual(list(call_variants.caller('CGGGGGCCATCA-','CGGGGGCCATCAG',softmask=True)),[('D', 12, 13, 'G')])
        self.assertEqual(list(call_variants.caller('CGGGGGCCATCAT','CGGGGGCCATCA-',softmask=True)),[('I', 12, 13, 'T')])
        
        
        #unequal sequence lengths
        #self.assertEqual(list(call_variants.caller('CGGGGGCCATCAT','cggg-TTCAT',softmask=True)),[]) #as expected this line raises a runtime error
        #however the unit test below fails.  This in a testing error, not a code error so marked as #pragma no cover
        #self.assertRaises(RuntimeError, call_variants.caller,'CGGGGGCCATCAT','cggg-TTCAT',softmask=True)
        
        pass

    def test_call_variants(self):
        self.maxDiff = None
        #                        555555555566666666667777777777888888888899999999990000000000111111111122222222222333333333344444444444444
        #                        012345678901234567890123456789012345678901234567890123456789012345678901234567899012345678901234555567890
        #                                                                                 v              v       v   vvv          vvv
        result = call_variants.call_variants( 'CAGGACTGGCTGCCGGCCCTTCTCTCCAGGTACTGGCCCCACGGCCTGAAGACTTCACGCGGCCCAGACGTG-TCAGCGGGCAG---GTACCCCGGGCATGTGCA',#mutant
                                'CAGGACTGGCTGCCGGCCCTTCTCTCCAGGTACTGGCCCCACGGCCTGAAGACTTCATGCGGCCCAGACGTGTTCAGCGG-CAGCTCGTACCCCGGG---GTGCA',#ref
                                'X',
                                ref_start = 153457150)
        self.assertEqual(repr(result),"[AmpliconVariant(chromosome='X', start=153457207, end=153457207, vcf_start=153457207, ref_allele='T', alt_allele='C'), " +\
                                "AmpliconVariant(chromosome='X', start=153457222, end=153457222, vcf_start=153457221, ref_allele='GT', alt_allele='G'), " +\
                                "AmpliconVariant(chromosome='X', start=153457230, end=153457230, vcf_start=153457229, ref_allele='G', alt_allele='GG'), " +\
                                "AmpliconVariant(chromosome='X', start=153457233, end=153457235, vcf_start=153457232, ref_allele='GCTC', alt_allele='G'), " +\
                                "AmpliconVariant(chromosome='X', start=153457246, end=153457246, vcf_start=153457245, ref_allele='G', alt_allele='GCAT')]")

        #                        555555555566666666667777777777888888888899999999990000000000111111111122222222222333333333344444444444444
        #                        012345678901234567890123456789012345678901234567890123456789012345678901234567899012345678901234555567890
        #                                                                                 v              v       v   vvv          vvv
        result = call_variants.call_variants( 'CAGGACTGGCTGCCGGCCCTTCTCTCCAGGTACTGGCCCCACGGCCTGAAGACTTCACGCGGCCCAGACGTG-TCAGCGGGCAG---GTACCCCGGGCATGTGCA',#mutant
                                'caggactGGCTGCCGGCCCTTCTCTCCAGGTACTGGCCCCACGGCCTGAAGACTTCATGCGGCCCAGACGTGTTCAGCGG-CAGCTCGTACCCCGGG---GTGCA',#ref
                                'X',
                                ref_start = 153457150)
        self.assertEqual(repr(result),"[AmpliconVariant(chromosome='X', start=153457207, end=153457207, vcf_start=153457207, ref_allele='T', alt_allele='C'), " +\
                                "AmpliconVariant(chromosome='X', start=153457222, end=153457222, vcf_start=153457221, ref_allele='GT', alt_allele='G'), " +\
                                "AmpliconVariant(chromosome='X', start=153457230, end=153457230, vcf_start=153457229, ref_allele='G', alt_allele='GG'), " +\
                                "AmpliconVariant(chromosome='X', start=153457233, end=153457235, vcf_start=153457232, ref_allele='GCTC', alt_allele='G'), " +\
                                "AmpliconVariant(chromosome='X', start=153457246, end=153457246, vcf_start=153457245, ref_allele='G', alt_allele='GCAT')]")

        result = call_variants.call_variants( 'CAGGGGGACTGGCTGTCGGCC',#mutant
                                'cagG---ACTGGCTGCCGGCC',#ref
                                'X',
                                ref_start = 153457150)
        self.assertEqual(repr(result),"[AmpliconVariant(chromosome='X', start=153457154, end=153457154, vcf_start=153457153, ref_allele='G', alt_allele='GGGG'), " +\
                                "AmpliconVariant(chromosome='X', start=153457162, end=153457162, vcf_start=153457162, ref_allele='C', alt_allele='T')]")
        
        result = call_variants.call_variants( 'CAGGGACTGGCTGTCGGCC',#mutant
                                'cagg-actGGCTGCCGGCC',#ref
                                'X',
                                ref_start = 153457150)
        self.assertEqual(repr(result),"[AmpliconVariant(chromosome='X', start=153457162, end=153457162, vcf_start=153457162, ref_allele='C', alt_allele='T')]")
        
        result = call_variants.call_variants( 'CAGGGGGACTGGCTGTCGGCC',#mutant
                                'cagg---actGGCTGCCGGCC',#ref
                                'X',
                                ref_start = 153457150)
        self.assertEqual(repr(result),"[AmpliconVariant(chromosome='X', start=153457162, end=153457162, vcf_start=153457162, ref_allele='C', alt_allele='T')]")

        pass


    
    def test_call_variants_to_vcf(self):
        resultfile = StringIO()
        expected_result = "\n".join([
        "X\t153457207\t.\tT\tC\t.\tPASS\t.",
        "X\t153457221\t.\tGT\tG\t.\tPASS\t.",
        "X\t153457229\t.\tG\tGG\t.\tPASS\t.",
        "X\t153457232\t.\tGCTC\tG\t.\tPASS\t.",
        "X\t153457245\t.\tG\tGCAT\t.\tPASS\t.",
        ])+"\n"
        call_variants.call_variants_to_vcf('CAGGACTGGCTGCCGGCCCTTCTCTCCAGGTACTGGCCCCACGGCCTGAAGACTTCACGCGGCCCAGACGTG-TCAGCGGGCAG---GTACCCCGGGCATGTGCA',#mutant
                              'CAGGACTGGCTGCCGGCCCTTCTCTCCAGGTACTGGCCCCACGGCCTGAAGACTTCATGCGGCCCAGACGTGTTCAGCGG-CAGCTCGTACCCCGGG---GTGCA',#ref
                              'X',
                              ref_start = 153457150,
                              outfile = resultfile)
        self.assertEqual(resultfile.getvalue(), expected_result)
        pass
        

    
if __name__ == '__main__':
    unittest.main()