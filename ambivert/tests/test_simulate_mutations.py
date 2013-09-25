#!/usr/bin/env python
# encoding: utf-8
"""
test_simulate_mutations.py

Created by Matthew Wakefield on 2013-09-24.
Copyright (c) 2013 Matthew Wakefield and The University of Melbourne. All rights reserved.
"""

import unittest
import io
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

SAMFILE = u"""@SQ	SN:chrM	LN:16571
@SQ	SN:chr1	LN:249250621
@SQ	SN:chr2	LN:243199373
@SQ	SN:chr3	LN:198022430
@SQ	SN:chr4	LN:191154276
@SQ	SN:chr5	LN:180915260
@SQ	SN:chr6	LN:171115067
@SQ	SN:chr7	LN:159138663
@SQ	SN:chr8	LN:146364022
@SQ	SN:chr9	LN:141213431
@SQ	SN:chr10	LN:135534747
@SQ	SN:chr11	LN:135006516
@SQ	SN:chr12	LN:133851895
@SQ	SN:chr13	LN:115169878
@SQ	SN:chr14	LN:107349540
@SQ	SN:chr15	LN:102531392
@SQ	SN:chr16	LN:90354753
@SQ	SN:chr17	LN:81195210
@SQ	SN:chr18	LN:78077248
@SQ	SN:chr19	LN:59128983
@SQ	SN:chr20	LN:63025520
@SQ	SN:chr21	LN:48129895
@SQ	SN:chr22	LN:51304566
@SQ	SN:chrX	LN:155270560
@SQ	SN:chrY	LN:59373566
chr17_41243125_150M_27A122_chr17_41243200_150M_150	99	chr17	41243125	60	150M	=	41243200	225	CATTGGCAGGACTGGATTTACTTTCATGACACACAAAATGATTAAATTCCTTGCTTTGGGACACCTGGATTTGCTTTTATAAAATGAAACCAGAAGTAAGTCCACCAGTAATTAGGATGTTAAAGCTCATTCAGTCAAAGATGACGTCCT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:1	AS:i:145	XS:i:19
chr17_41243125_150M_27A122_chr17_41243200_150M_150	147	chr17	41243200	60	150M	=	41243125	-225	TTTATAAAATGAAACCAGAAGTAAGTCCACCAGTAATTAGGATGTTAAAGCTCATTCAGTCAAAGATGACGTCCTAGCTGTGTGAAGGACTTTTTTCTATGAAAAGCACCTTAGGAGGAACATGTTTCAAGTTTAAGAAGCAGTTCCTTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:0	AS:i:150	XS:i:0
chr13_32900016_120M_20T99_chr13_32900109_121M_121	99	chr13	32899987	60	150M	=	32900109	272	AATTTTTGAGTTTAAAATACACGGTTTCCAGCAGCTGAAATTTGTGAGTTCATATGTGTTGGCATTTTAAACATCACTTGATGATTATTTAATGCTTCATGAGAGATTTACTTTTTAAAATGTAATATAAAATATCTAAAAGTAGTATTC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:1	AS:i:145	XS:i:0
chr13_32900016_120M_20T99_chr13_32900109_121M_121	147	chr13	32900109	60	150M	=	32899987	-272	TAATATAAAATATCTAAAAGTAGTATTCCAACAATTTATATGAATGAGAATCTTCTTTTAAAAATAAGATAAACTAGTTTTTGCCAGTTTTTTAAAATAACCTAAGGGATTTGCTTTGTTTTATTTTAGTCCTGTTGTTCTACAATGTAC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:0	AS:i:150	XS:i:0
chr17_41243125_150M_27C122_chr17_41243200_150M_150	99	chr13	32899987	60	150M	=	32900109	272	AATTTTTGAGTTTAAAATACACGGTTTCCAGCAGCTGAAATTTGTGAGTAAATATGTGTTGGCATTTTAAACATCACTTGATGATTATTTAATGCTTCATGAGAGATTTACTTTTTAAAATGTAATATAAAATATCTAAAAGTAGTATTC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:1	AS:i:145	XS:i:20
chr17_41243125_150M_27C122_chr17_41243200_150M_150	147	chr13	32900109	60	150M	=	32899987	-272	TAATATAAAATATCTAAAAGTAGTATTCCAACAATTTATATGAATGAGAATCTTCTTTTAAAAATAAGATAAACTAGTTTTTGCCAGTTTTTTAAAATAACCTAAGGGATTTGCTTTGTTTTATTTTAGTCCTGTTGTTCTACAATGTAC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:0	AS:i:150	XS:i:0
chr17_41243125_150M_27C122_chr17_41243200_150M_150	0	.	.	.	.	.	.	.	TAATATAAAATATCTAAAAGTAGTATTCCAACAATTTATATGAATGAGAATCTTCTTTTAAAAATAAGATAAACTAGTTTTTGCCAGTTTTTTAAAATAACCTAAGGGATTTGCTTTGTTTTATTTTAGTCCTGTTGTTCTACAATGTAC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""

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
        self.assertEqual(mutated_amplicon_to_paired_reads('cagtGATCGATCacgt','chrX','12349','8M','3C4', readlength=10, position_excludes_softmasked=True),
                        (('cagtGATCGA', '', 'chrX', '12349', '6M', '3C2'), ('acgtGATCGA', '', 'chrX', '12351', '6M', '1C4')))
        pass
    
    def test_check_point_mutate_sequence(self):
        samfile = io.StringIO(SAMFILE)
        outfile = io.StringIO()
        check_point_mutate_sequence(samfile, outfile=outfile)
        print(outfile.getvalue())
        self.assertEqual(md5(outfile.getvalue()).hexdigest(),'de119c1e94dbf6e1664b419f09500c07')
        pass
    

if __name__ == '__main__':
    unittest.main()