#!/usr/bin/env python3
# encoding: utf-8
"""
test_truseq_manifest.py

Created by Matthew Wakefield.
Copyright (c) 2013-2015  Matthew Wakefield and The University of Melbourne. All rights reserved.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

"""

import io
import unittest
from string import ascii_uppercase, ascii_lowercase
from ambivert.sequence_utilities import *

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2013-2015,  Matthew Wakefield and The University of Melbourne"
__credits__ = ["Matthew Wakefield",]
__license__ = "GPLv3"
__version__ = "0.5.1"
__maintainer__ = "Matthew Wakefield"
__email__ = "matthew.wakefield@unimelb.edu.au"
__status__ = "Development"

SAMDATA = """@HQ\tVN:1.5\tSO:unsorted
@SQ\tSN:chr1\tLN:249250621
@SQ\tSN:chr2\tLN:243199373
@PG	ID:AmBiVErT	VN:0.1.11
1186cf52bcce24d9e5446675e896b90c_1	0	chr17	41244429	40	3S16M2S	*	0	16	CTGGAACCTACTTCATTAATA	IIIIIIIIIIIIIIIIIIIII
1186cf52bcce24d9e5446675e896b90c_2	0	chr17	41244429	40	3S16M2S	*	0	16	CTGGAACCTACTTCATTAATA	IIIIIIIIIIIIIIIIIIIII
1186cf52bcce24d9e5446675e896b90c_3	0	chr17	41244429	40	3S16M2S	*	0	16	CTGGAACCTACTTCATTAATA	IIIIIIIIIIIIIIIIIIIII
HWI-ST960:63:D0CYJACXX:4:1201:14988:16838	16	MT	1	3	21S29M	*	0	0	CTTATTTAAGGGGAACGTGTGGATCACAGGTCTATCACCCTATTAACCAC	AIGEHGDF>IIHEHGAIGIIIGGEEFHIIGHHFAFEAHHDHHFFDDFCC@	AS:i:58	XS:i:44XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:29	YT:Z:UU
"""

class test_sequence_utilities(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_open_potentially_gzipped_file(self):
        pass
        
    def test_parse_fastq(self):
        testdata = io.BytesIO(b'\n'.join([b"@firstreadname",b"GACAT",b'+',b'IAB+@',
                               b"@secondreadname:1_2 1:2",b"GNATCAT",b'+',b'IAB+@EE'])+b'\n')
        result = list(parse_fastq(testdata))
        self.assertEqual(result, [('firstreadname', 'GACAT', 'IAB+@'), ('secondreadname:1_2 1:2', 'GNATCAT', 'IAB+@EE')])
    
    def test_parse_fasta(self):
        testdata = io.StringIO('>firstsequence\nGACAT\n>secondsequence\nGNATCAT')
        self.assertEqual(list(parse_fasta(testdata)),[('firstsequence', 'GACAT'), ('secondsequence', 'GNATCAT')])
    
    def test_complement(self):
        self.assertEqual(complement('ACBDGHKMmNSRUTWVYnacbdghk-.srutwvy'),'TGVHCDMKkNSYAAWBRntgvhcdm-.syaawbr')
        self.assertEqual(complement('ACBDGHKMmNSRUTWVYInacbdghk-.srutwvy'),'TGVHCDMKkNSYAAWBRnntgvhcdm-.syaawbr')

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement('ACBDGHKMmNSRUTWVYnacbdghk-.srutwvy'),'rbwaays.-mdchvgtnRBWAAYSNkKMDCHVGT')
    
    def test_encode_ambiguous(self):
        self.assertEqual([encode_ambiguous(x) for x in ['CG','AGT','ACG','GT','CGT','AT','ACT','ACGT','AC','AG','CT','ANT','A','C','T','G','N']],
                        ['S', 'D', 'V', 'K', 'B', 'W', 'H', 'N', 'M', 'R', 'Y', 'N','A','C','T','G','N'])
    
    def test_flatten_paired_alignment(self):
        testdata = [('GATC---',
                     '---GATC',
                     'GATSATC'),
                    ('CAGGACTGGCTGCCGGCCCTTCTCTCCAGGTACTGGCCCCACGGCCTGAAGA----------------------------',
                     '--------------------------------CTGGCCCCACGGCCTGAAGACTTCATGCGGCCCAGACGTGTTCAGCGG',
                     'CAGGACTGGCTGCCGGCCCTTCTCTCCAGGTACTGGCCCCACGGCCTGAAGACTTCATGCGGCCCAGACGTGTTCAGCGG'),
                    ('--------------------------------CTGGCCCCACGGCCTGAAGACTTCATGCGGCCCAGACGTGTTCAGCGG',
                     'CAGGACTGGCTGCCGGCCCTTCTCTCCAGGTACTGGCCCCACGGCCTGAAGA----------------------------',
                     'CAGGACTGGCTGCCGGCCCTTCTCTCCAGGTACTGGCCCCACGGCCTGAAGACTTCATGCGGCCCAGACGTGTTCAGCGG'),
                    ('CAGGACTGGCTGCCGGCCCTTCTCTCCAGGTACTGGCC---CGGCCTGAAGA----------------------------',
                     '--------------------------------CTGGCCCCACGGCCTGAAGACTTCATGCGGCCCAGACGTGTTCAGCGG',
                     'CAGGACTGGCTGCCGGCCCTTCTCTCCAGGTACTGGCCCCACGGCCTGAAGACTTCATGCGGCCCAGACGTGTTCAGCGG'),
                    ]
        for (seq1,seq2,flat) in testdata:
            self.assertEqual(flatten_paired_alignment(seq1,seq2),flat)
    
    def test_make_blocklist(self):
        canned_result = ['abcdefghij', 'klmnopqrst', 'uvwxyzABCD',
                        'EFGHIJKLMN', 'OPQRSTUVWX', 'YZ']
        result = make_blocklist(ascii_lowercase+ascii_uppercase, 10)
        self.assertEqual(canned_result, result)
    
    def test_slice_string_in_blocks(self):
        canned_result = 'abcdefghij\nklmnopqrst\nuvwxyzABCD\nEFGHIJKLMN\nOPQRSTUVWX\nYZ\n'
        result = slice_string_in_blocks(ascii_lowercase+ascii_uppercase, 10)
        self.assertEqual(canned_result, result)
    
    def test_format_fasta(self):
        canned_result = '>fooGene\nabcdefghij\nklmnopqrst\nuvwxyzABCD\nEFGHIJKLMN\nOPQRSTUVWX\nYZ\n'
        result = format_fasta('fooGene',ascii_lowercase+ascii_uppercase, 10)
        self.assertEqual(canned_result, result)
    
    def test_get_sam_header(self):
        self.assertEqual(get_sam_header(io.StringIO(SAMDATA)),
                ['@HQ\tVN:1.5\tSO:unsorted', '@SQ\tSN:chr1\tLN:249250621', '@SQ\tSN:chr2\tLN:243199373', '@PG\tID:AmBiVErT\tVN:0.1.11'])
    
    def test_get_tag(self):
        data = SAMDATA.split('\n')[5:]
        self.assertEqual([get_tag(x.split('\t'),'MD') for x in data],[None, None, '29', None])
    
    def test_fix_softmasked_expanded_cigar(self):
        self.assertEqual(fix_softmasked_expanded_cigar(['S', 'S', 'S', 'S', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'S', 'S', 'S', 'S']),
                        ('SSSSMMMMMMMMSSSS', 4, 8))
        self.assertEqual(fix_softmasked_expanded_cigar(['S', 'S', 'D', 'S', 'M', 'M', 'M', 'I', 'M', 'D', 'M', 'M', 'S', 'I', 'S', 'S']),
                        ('SSSMMMIMDMMSSSS', 3, 8))
        self.assertEqual(fix_softmasked_expanded_cigar(['S', 'S', 'D', 'S', 'M', 'M', 'M', 'I', 'M', 'D', 'M', 'M']),
                        ('SSSMMMIMDMM', 3, 8))
        self.assertEqual(fix_softmasked_expanded_cigar(['M', 'M', 'M', 'I', 'M', 'D', 'M', 'M', 'S', 'I', 'S', 'S']),
                        ('MMMIMDMMSSSS', 0, 8))
        self.assertEqual(fix_softmasked_expanded_cigar(['M', 'M', 'M', 'I', 'M', 'D', 'M', 'M', 'M', 'M', 'M', 'M']),
                        ('MMMIMDMMMMMM', 0, 11))
        pass
    
    def test_gapped_alignment_to_cigar(self):
        self.assertEqual(gapped_alignment_to_cigar('gtacACGTACGTgtac','GTACACGTACGTGTAC'),
                        ('4S8M4S', 4, 8)
                        )
        self.assertEqual(gapped_alignment_to_cigar('gtacACG-ACGTg-ac','GT-CACGTA-GTGTAC'),
                        ('3S3M1I1M1D2M4S', 3, 8)
                        )
        self.assertEqual(gapped_alignment_to_cigar('gtacACGAACGTg-ac','GT-CACGTA-GTGTAC',snv='X'),
                        ('3S3M1X1M1D2M4S', 3, 8)
                        )
        
        self.assertRaises(RuntimeError,gapped_alignment_to_cigar,'gtacACG-ACGTg-ac','GT-CACGTA-GTG')
        self.assertRaises(RuntimeError,gapped_alignment_to_cigar,'gtacACG-ACGTg-ac','gtacACG-ACGTg-ac')
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