#!/usr/bin/env python
# encoding: utf-8
"""
test_simulate_variants.py

Created by Matthew Wakefield.
Copyright (c) 2013-2014  Matthew Wakefield and The University of Melbourne. All rights reserved.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

"""

import unittest
import io
from ambivert.simulate_variants import *
import hashlib

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2013-2015,  Matthew Wakefield and The University of Melbourne"
__credits__ = ["Matthew Wakefield",]
__license__ = "GPLv3"
__version__ = "0.5.1"
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
def md5(data):
    #python3 compatibility
    return hashlib.md5(bytes(data,'ascii'))

class test_simulate_variants(unittest.TestCase):
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
    
    def test_make_all_point_variants(self):
        sequence = 'cagtGATCGATCacgt'
        #print(list(make_all_point_variants(sequence)))
        self.assertEqual('63d7c766aebd4de2d4f6ac93ca0c2550',md5(str(list(make_all_point_variants(sequence)))).hexdigest())
        
        #print(list(make_all_point_variants(sequence, skip_softmasked=False, position_excludes_softmasked=False)))
        self.assertEqual('e39f6ab49cf9fa517368a8dcba416c97',md5(str(list(make_all_point_variants(sequence, skip_softmasked=False, position_excludes_softmasked=False)))).hexdigest())
        pass
    
    def test_format_fastq_entry(self):
        #name,sequence,q_generator=constant_quality
        self.assertEqual('@None_1_16M_4A11\ncagtAATCGATCacgt\n+\nIIIIIIIIIIIIIIII\n',
                        format_fastq_entry(*('None_1_16M_4A11', 'cagtAATCGATCacgt'))
                        )
        pass
    
    def test_cigar_add_deletion(self):
        self.assertEqual(cigar_add_deletion('10M',one_based_start=5),'4M1D5M')
        self.assertEqual(cigar_add_deletion('10M',one_based_start=5, length=5),'4M5D1M')
        self.assertEqual(cigar_add_deletion('2S4M5D2M2I10M',one_based_start=5),'2S4M6D1M2I10M')
        self.assertRaises(AssertionError,cigar_add_deletion,'2S4M5D2M2I10M',one_based_start=6,length=3)
        self.assertRaises(AssertionError,cigar_add_deletion,'10M',one_based_start=6,length=10)
        
    def test_variant_detection_tag_trimmer(self):
        self.assertEqual(variant_detection_tag_trimmer('7A8',trim_from_start=0,trim_from_end=0),'7A8')
        self.assertEqual(variant_detection_tag_trimmer('7A8',trim_from_start=3,trim_from_end=0),'4A8')
        self.assertEqual(variant_detection_tag_trimmer('7A8',trim_from_start=0,trim_from_end=5),'7A3')
        self.assertEqual(variant_detection_tag_trimmer('7A8',trim_from_start=3,trim_from_end=5),'4A3')
        self.assertEqual(variant_detection_tag_trimmer('7A8G',trim_from_start=3,trim_from_end=5),'4A4')
        self.assertEqual(variant_detection_tag_trimmer('17',trim_from_start=3,trim_from_end=5),'9')
        pass
    
    def test_mdtag_add_deletion(self):
        self.assertEqual(mdtag_add_deletion('10',seq='GGGGGGGGGG',one_based_start=3,length=1),'2^G7')
        self.assertEqual(mdtag_add_deletion('10',seq='GGTGGGGGGG',one_based_start=3,length=1),'2^T7')
        self.assertEqual(mdtag_add_deletion('7A8G',seq='TTTTTTTTTTTTTTTTT',one_based_start=3,length=1),'2^T4A8G')
        self.assertEqual(mdtag_add_deletion('7A8G',seq='TTTTTTTTTTTTTTTTT',one_based_start=7,length=1),'6^T0A8G')
        self.assertEqual(mdtag_add_deletion('7A8G',seq='TTTTTTTTTTTTTTTTT',one_based_start=7,length=2),'6^TT8G')
        self.assertEqual(mdtag_add_deletion('2^CAT8',seq='GGGGGGGGGG',one_based_start=8,length=1),'2^CAT5^G2')
        self.assertEqual(mdtag_add_deletion('10',seq='GGGGGGGGGG',one_based_start=3,length=1),'2^G7')
        self.assertEqual(mdtag_add_deletion('10',seq='GGGGGGGGGG',one_based_start=3,length=1),'2^G7')
        self.assertEqual(mdtag_add_deletion('10',seq='GGGGGGGGGG',one_based_start=3,length=1),'2^G7')
        pass
    
    def test_deletion_mutate_sequence(self):
        #sequence,chromosome=None,one_based_start=1,one_based_site=1, length=1, position_excludes_softmasked=True
        self.assertEqual(deletion_mutate_sequence('acgtCATGacgt',one_based_start=1,one_based_site=6, length=1),('None_5_1M1D2M_1^A2', 'acgtCTGacgt'))############# TODO WORKING ON THIS
        pass
    
    def test_make_all_deletion_variants(self):
        self.assertEqual(list(make_all_deletion_variants('acgtCATGacgt',length=2)),[('None_5_2D2M_^CA2', 'acgtTGacgt'),
                                                                                    ('None_5_1M2D1M_1^AT1', 'acgtCGacgt'),
                                                                                    ('None_5_2M2D_2^TG', 'acgtCAacgt'),])
        pass

    
    def test_mutated_amplicon_to_paired_reads(self):
        self.assertEqual(mutated_amplicon_to_paired_reads('cagtGATCGATCacgt','chrX','12345','16M','7C8',readlength=10),
                        (('cagtGATCGA', '', 'chrX', '12345', '10M', '7C2'), ('acgtGATCGA', '', 'chrX', '12351', '10M', '1C8')))
        self.assertEqual(mutated_amplicon_to_paired_reads('cagtGATCGATCacgt','chrX','12345','16M','7C8'),
                        (('cagtGATCGATCacgt', '', 'chrX', '12345', '16M', '7C8'), ('acgtGATCGATCactg', '', 'chrX', '12345', '16M', '7C8')))
        self.assertEqual(mutated_amplicon_to_paired_reads('cagtGATCGATCacgt','chrX','12345','16M','7C8',quality='ABCGEFGHIJKLMNOP'),
                         (('cagtGATCGATCacgt', 'ABCGEFGHIJKLMNOP', 'chrX', '12345', '16M', '7C8'),
                          ('acgtGATCGATCactg', 'PONMLKJIHGFEGCBA', 'chrX', '12345', '16M', '7C8')))
        self.assertEqual(mutated_amplicon_to_paired_reads('cagtGATCGATCacgtt','chrX','12349','8M','3C4', readlength=10, position_excludes_softmasked=True),
                        (('cagtGATCGA', '', 'chrX', '12349', '6M', '3C2'), ('aacgtGATCG', '', 'chrX', '12352', '5M', 'C4')))
        pass
    
    def test_check_point_mutated_sequence(self):
        samfile = io.StringIO(SAMFILE)
        outfile = io.StringIO()
        check_point_mutated_sequence(samfile, outfile=outfile, verbose=False)
        print(outfile.getvalue())
        self.assertEqual(md5(outfile.getvalue()).hexdigest(),'17c5c2f058c4ff56a173d0255314ba6a')
        pass
    

if __name__ == '__main__':
    unittest.main()