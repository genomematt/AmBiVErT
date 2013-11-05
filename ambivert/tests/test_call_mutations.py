#!/usr/bin/env python
# encoding: utf-8
"""
test_call_mutations.py

Created by Matthew Wakefield on 2013-11-04.
Copyright (c) 2013 Matthew Wakefield and The University of Melbourne. All rights reserved.
"""

from __future__ import print_function
import unittest
from StringIO import StringIO
from ambivert.call_mutations import *


class test_call_mutations(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_caller(self):
        self.assertEqual(list(caller('AAAA','AAGA')),[('X', 2, 3, 'A')])
        self.assertEqual(list(caller('AAAA','AGGA')),[('X', 1, 2, 'A'), ('X', 2, 3, 'A')])
        self.assertEqual(list(caller('AAAA','AA-A')),[('I', 2, 3, 'A')])
        self.assertEqual(list(caller('AA-A','AAAA')),[('D', 2, 3, 'A')])
        self.assertEqual(list(caller('AAAA','ACGA')),[('X', 1, 2, 'A'), ('X', 2, 3, 'A')])
        self.assertEqual(list(caller('AATA','AC-A')),[('X', 1, 1, 'A'), ('I', 2, 3, 'T')])
        self.assertEqual(list(caller('AA-A','ACAA')),[('X', 1, 2, 'A'), ('D', 2, 3, 'A')])
        self.assertEqual(list(caller('AAAAA','AAGCA')),[('X', 2, 3, 'A'), ('X', 3, 4, 'A')])
        self.assertEqual(list(caller('AAAAA','AA-CA')),[('I', 2, 3, 'A'), ('X', 3, 3, 'A')])
        self.assertEqual(list(caller('AA-AA','AAACA')),[('D', 2, 3, 'A'), ('X', 3, 4, 'A')])
        self.assertEqual(list(caller('GAAAA','gAAGA')),[('X', 3, 4, 'A')])
        self.assertEqual(list(caller('GAAAA','gAA-A')),[('I', 3, 4, 'A')])
        self.assertEqual(list(caller('GAA-A','gAAAA')),[('D', 3, 4, 'A')])
        self.assertEqual(list(caller('CAA-A','gAAAA')),[('D', 3, 4, 'A')])
        self.assertEqual(list(caller('CAA-A','gAAAA',softmask=False)),[('X', 0, 1, 'C'), ('D', 3, 4, 'A')])
        self.assertEqual(list(caller('CAA-ATCT','gAAAAttt')),[('D', 3, 4, 'A')])
        self.assertEqual(list(caller('CAA-ATCT','gAAAAttt',softmask=False)),[('X', 0, 1, 'C'), ('D', 3, 4, 'A'), ('X', 6, 7, 'C')])
        self.assertEqual(list(caller('CAA-AT-T','gAAAAttt')),[('D', 3, 4, 'A')])
        self.assertEqual(list(caller('CAA-AT-T','gAAAAttt',softmask=False)),[('X', 0, 1, 'C'), ('D', 3, 4, 'A'), ('D', 6, 7, 't')])
        self.assertEqual(list(caller('CAA-ATTT','gAAAAt-t')),[('D', 3, 4, 'A')])
        self.assertEqual(list(caller('CAA-ATTT','gAAAAt-t',softmask=False)),[('X', 0, 1, 'C'), ('D', 3, 4, 'A'), ('I', 6, 7, 'T')])
        pass
    
    def test_call_mutations_to_vcf(self):
        resultfile = StringIO()
        expected_result = "\n".join([
        "X\t153457207\t.\tT\tC\t.\tPASS\t.",
        "X\t153457221\t.\tGT\tG\t.\tPASS\t.",
        "X\t153457229\t.\tG\tGG\t.\tPASS\t.",
        "X\t153457232\t.\tGCTC\tG\t.\tPASS\t.",
        "X\t153457245\t.\tG\tGCAT\t.\tPASS\t.",
        ])+"\n"
        #                      555555555566666666667777777777888888888899999999990000000000111111111122222222222333333333344444444444444
        #                      012345678901234567890123456789012345678901234567890123456789012345678901234567899012345678901234555567890
        #                                                                               v              v       v   vvv          vvv
        call_mutations_to_vcf('CAGGACTGGCTGCCGGCCCTTCTCTCCAGGTACTGGCCCCACGGCCTGAAGACTTCACGCGGCCCAGACGTG-TCAGCGGGCAG---GTACCCCGGGCATGTGCA',#mutant
                              'CAGGACTGGCTGCCGGCCCTTCTCTCCAGGTACTGGCCCCACGGCCTGAAGACTTCATGCGGCCCAGACGTGTTCAGCGG-CAGCTCGTACCCCGGG---GTGCA',#ref
                              'X',
                              ref_start = 153457150,
                              outfile = resultfile)
        print(resultfile.getvalue())
        self.assertEqual(resultfile.getvalue(), expected_result)
        pass
        

    
if __name__ == '__main__':
    unittest.main()