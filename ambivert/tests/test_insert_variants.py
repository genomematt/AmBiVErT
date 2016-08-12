#!/usr/bin/env python
# encoding: utf-8
"""
test_insert_variants.py

Created by Matthew Wakefield on 2013-09-26.
Copyright (c) 2013-2015  Matthew Wakefield and The University of Melbourne. All rights reserved.
"""

import unittest
import io
from ambivert.insert_variants import *
from ambivert.sequence_utilities import reverse_complement

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2013-2015,  Matthew Wakefield and The University of Melbourne"
__credits__ = ["Matthew Wakefield",]
__license__ = "GPLv3"
__version__ = "0.5.1"
__maintainer__ = "Matthew Wakefield"
__email__ = "matthew.wakefield@unimelb.edu.au"
__status__ = "Development"

SAMFILE = u"""@SQ	SN:chr13	LN:115169878
@SQ	SN:chr17	LN:81195210
chr17_41243153_122M_15G106_chr17_41243200_123M_123	99	chr17	41243153	60	122M	=	41243200	170	TCACACAAAATGATTGAATTCCTTGCTTTGGGACACCTGGATTTGCTTTTATAAAATGAAACCAGAAGTAAGTCCACCAGTAATTAGGATGTTAAAGCTCATTCAGTCAAAGATGACGTCCT	0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK	NM:i:1	AS:i:117	XS:i:19
chr17_41243153_122M_15G106_chr17_41243200_123M_123	147	chr17	41243200	60	123M	=	41243153	-170	TTTATAAAATGAAACCAGAAGTAAGTCCACCAGTAATTAGGATGTTAAAGCTCATTCAGTCAAAGATGACGTCCTAGCTGTGTGAAGGACTTTTTTCTATGAAAAGCACCTTAGGAGGAACAT	0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK	NM:i:0	AS:i:123	XS:i:0
chr17_41243153_122M_27G94_chr17_41243200_123M_123	147	chr17	41243200	60	123M	=	41243153	-170	TTTATAAAATGAAACCAGAAGTAAGTCCACCAGTAATTAGGATGTTAAAGCTCATTCAGTCAAAGATGACGTCCTAGCTGTGTGAAGGACTTTTTTCTATGAAAAGCACCTTAGGAGGAACAT	0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK	NM:i:0	AS:i:123	XS:i:0
chr13_32900227_120M_79G40_chr13_32900272_123M_34G88	99	chr13	32900227	60	120M	=	32900272	168	TTTTATTTTAGTCCTGTTGTTCTACAATGTACACATGTAACACCACAAAGAGATAAGTCAGGTATGATTAAAAACAATGGTTTTTATTCTTAGAATACTAGAAATGTTAATAAAAATAAA	0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK	NM:i:1	AS:i:115	XS:i:21
chr13_32900227_120M_79G40_chr13_32900272_123M_34G88	147	chr13	32900272	60	123M	=	32900227	-168	CAAAGAGATAAGTCAGGTATGATTAAAAACAATGGTTTTTATTCTTAGAATACTAGAAATGTTAATAAAAATAAAACTTAACAATTTTCCCCTTTTTTTACCCCCAGTGGTATGTGGGAGTTT	0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK	NM:i:1	AS:i:118	XS:i:20
chr13_32900227_120M_120_chr13_32900272_123M_122C	147	chr13	32900272	60	123M	=	32900227	-168	CAAAGAGATAAGTCAGGTATGATTAAAAACAATGCTTTTTATTCTTAGAATACTAGAAATGTTAATAAAAATAAAACTTAACAATTTTCCCCTTTTTTTACCCCCAGTGGTATGTGGGAGTTC	0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK	NM:i:1	AS:i:122	XS:i:0
"""

SAMFILE_SNP_CHR13 = u"""@SQ	SN:chr13	LN:115169878
@SQ	SN:chr17	LN:81195210
chr17_41243153_122M_15G106_chr17_41243200_123M_123	99	chr17	41243153	60	122M	=	41243200	170	TCACACAAAATGATTGAATTCCTTGCTTTGGGACACCTGGATTTGCTTTTATAAAATGAAACCAGAAGTAAGTCCACCAGTAATTAGGATGTTAAAGCTCATTCAGTCAAAGATGACGTCCT	0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK	NM:i:1	AS:i:117	XS:i:19
chr17_41243153_122M_15G106_chr17_41243200_123M_123	147	chr17	41243200	60	123M	=	41243153	-170	TTTATAAAATGAAACCAGAAGTAAGTCCACCAGTAATTAGGATGTTAAAGCTCATTCAGTCAAAGATGACGTCCTAGCTGTGTGAAGGACTTTTTTCTATGAAAAGCACCTTAGGAGGAACAT	0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK	NM:i:0	AS:i:123	XS:i:0
chr17_41243153_122M_27G94_chr17_41243200_123M_123	147	chr17	41243200	60	123M	=	41243153	-170	TTTATAAAATGAAACCAGAAGTAAGTCCACCAGTAATTAGGATGTTAAAGCTCATTCAGTCAAAGATGACGTCCTAGCTGTGTGAAGGACTTTTTTCTATGAAAAGCACCTTAGGAGGAACAT	0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK	NM:i:0	AS:i:123	XS:i:0
chr13_32900227_120M_79G40_chr13_32900272_123M_34G88	99	chr13	32900227	60	120M	=	32900272	168	TTTTATTTTAGTCCTGTTGTTCTACAATGTACACATGTAACACCACAAAGAGATAAGTCAGGTATGATTAAAAACAATGGTTTTTATTCTTAGAATACTAGAAATGTTAATAAAAATAAA	0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK	NM:i:1	AS:i:115	XS:i:21
chr13_32900227_120M_79G40_chr13_32900272_123M_34G88	147	chr13	32900272	60	123M	=	32900227	-168	CAAAGAGATAAGTCAGGTATGATTAAAAACAATGGTTTTTATTCTTAGAATACTAGAAATGTTAATAAAAATAAAACTTAACAATTTTCCCCTTTTTTTACCCCCAGTGGTATGTGGGAGTTT	0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK	NM:i:1	AS:i:118	XS:i:20
chr13_32900227_120M_120_chr13_32900272_123M_122C	147	chr13	32900272	60	123M	=	32900227	-168	CAAAGAGATAAGTCAGGTATGATTAAAAACAATGCTTTTTATTCTTAGAATACTAGAAATGTTAATAAAAATAAAACTTAACAATTTTCCCCTTTTTTTACCCCCAGTGGTATGTGGGAGTTC	0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK	NM:i:1	AS:i:122	XS:i:0
"""

SAMFILE_INDEL_CHR13 = u"""@SQ	SN:chr13	LN:115169878
@SQ	SN:chr17	LN:81195210
chr17_41243153_122M_15G106_chr17_41243200_123M_123	99	chr17	41243153	60	122M	=	41243200	170	TCACACAAAATGATTGAATTCCTTGCTTTGGGACACCTGGATTTGCTTTTATAAAATGAAACCAGAAGTAAGTCCACCAGTAATTAGGATGTTAAAGCTCATTCAGTCAAAGATGACGTCCT	0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK	NM:i:1	AS:i:117	XS:i:19
chr17_41243153_122M_15G106_chr17_41243200_123M_123	147	chr17	41243200	60	123M	=	41243153	-170	TTTATAAAATGAAACCAGAAGTAAGTCCACCAGTAATTAGGATGTTAAAGCTCATTCAGTCAAAGATGACGTCCTAGCTGTGTGAAGGACTTTTTTCTATGAAAAGCACCTTAGGAGGAACAT	0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK	NM:i:0	AS:i:123	XS:i:0
chr17_41243153_122M_27G94_chr17_41243200_123M_123	147	chr17	41243200	60	123M	=	41243153	-170	TTTATAAAATGAAACCAGAAGTAAGTCCACCAGTAATTAGGATGTTAAAGCTCATTCAGTCAAAGATGACGTCCTAGCTGTGTGAAGGACTTTTTTCTATGAAAAGCACCTTAGGAGGAACAT	0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK	NM:i:0	AS:i:123	XS:i:0
chr13_32900227_120M_79G40_chr13_32900272_123M_34G88	99	chr13	32900227	60	120M	=	32900272	168	TTTTATTTTAGTCCTGTTGTTCTACAATGTACACATGTAACACCACAAAGAGATAAGTCAGGTATGATTAAAAACAATGGTTTTTATTCTTAGAATACTAGAAATGTTAATAAAAATAAA	0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK	NM:i:1	AS:i:115	XS:i:21
chr13_32900227_120M_79G40_chr13_32900272_123M_34G88	147	chr13	32900272	60	123M	=	32900227	-168	CAAAGAGATAAGTCAGGTATGATTAAAAACAATGGTTTTTATTCTTAGAATACTAGAAATGTTAATAAAAATAAAACTTAACAATTTTCCCCTTTTTTTACCCCCAGTGGTATGTGGGAGTTT	0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK	NM:i:1	AS:i:118	XS:i:20
chr13_32900227_120M_120_chr13_32900272_123M_122C	147	chr13	32900272	60	123M	=	32900227	-168	CAAAGAGATAAGTCAGGTATGATTAAAAACAATGCTTTTTATTCTTAGAATACTAGAAATGTTAATAAAAATAAAACTTAACAATTTTCCCCTTTTTTTACCCCCAGTGGTATGTGGGAGTTC	0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK	NM:i:1	AS:i:122	XS:i:0
"""


#111111111111111111111111111111111111111111111111122222===22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
#555555555666666666677777777778888888888999999999900000===00000111111111122222222223333333333444444444455555555556666666666777777777788888888889999999999
#123456789012345678901234567890123456789012345678901234===56789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
#  TCACACAAAATGATTgAATTCCTTGCTTTGGGACACCTGGATTTGCTTTTATAAAAAAATGAAACCAGAAGTAAGTCCACCAGTAATTAGGATGTTAAAGCTCATTCAGTCAAAGATGACGT
#                                                    ATAAAAAAATGAAACCAGAAGTAAGTCCACCAGTAATTAGGATGTTAAAGCTCATTCAGTCAAAGATGACGTCCTAGCTGTGTGAAGGACTTTTTTCTATGAAAAGCACCTTAGGAGGAACAT
#
SAMFILE_INDEL_CHR17 = [
("chr17_41243153_52M3I67M_15G36^AAA67_chr17_41243200_2M3I118M_2^AAA118", "CACACAAAATGATTGAATTCCTTGCTTTGGGACACCTGGATTTGCTTTTATAAAAAAATGAAACCAGAAGTAAGTCCACCAGTAATTAGGATGTTAAAGCTCATTCAGTCAAAGATGACGT", """0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPIIIQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGH"""),
("chr17_41243153_52M3I67M_15G36^AAA67_chr17_41243200_2M3I118M_2^AAA118", "ATAAAAAAATGAAACCAGAAGTAAGTCCACCAGTAATTAGGATGTTAAAGCTCATTCAGTCAAAGATGACGTCCTAGCTGTGTGAAGGACTTTTTTCTATGAAAAGCACCTTAGGAGGAACAT", """34III56789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK"""),
("chr17_41243153_52M3I67M_27G24^AAA67_chr17_41243200_2M3I118M_2^AAA118", "ATAAAAAAATGAAACCAGAAGTAAGTCCACCAGTAATTAGGATGTTAAAGCTCATTCAGTCAAAGATGACGTCCTAGCTGTGTGAAGGACTTTTTTCTATGAAAAGCACCTTAGGAGGAACAT", """34III56789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK"""),
("chr13_32900227_120M_79G40_chr13_32900272_123M_34G88", "TTTTATTTTAGTCCTGTTGTTCTACAATGTACACATGTAACACCACAAAGAGATAAGTCAGGTATGATTAAAAACAATGGTTTTTATTCTTAGAATACTAGAAATGTTAATAAAAATAAA"	"""0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK"""),
("chr13_32900227_120M_79G40_chr13_32900272_123M_34G88", "CAAAGAGATAAGTCAGGTATGATTAAAAACAATGGTTTTTATTCTTAGAATACTAGAAATGTTAATAAAAATAAAACTTAACAATTTTCCCCTTTTTTTACCCCCAGTGGTATGTGGGAGTTT"	"""0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK"""),
("chr13_32900227_120M_120_chr13_32900272_123M_122C", "CAAAGAGATAAGTCAGGTATGATTAAAAACAATGCTTTTTATTCTTAGAATACTAGAAATGTTAATAAAAATAAAACTTAACAATTTTCCCCTTTTTTTACCCCCAGTGGTATGTGGGAGTTC"	"""0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK"""),
]

READS = [
        (
         ("chr17_41243153_122M_15G106_chr17_41243200_123M_123","41243153","122M","TCACACAAAATGATTGAATTCCTTGCTTTGGGACACCTGGATTTGCTTTTATAAAATGAAACCAGAAGTAAGTCCACCAGTAATTAGGATGTTAAAGCTCATTCAGTCAAAGATGACGTCCT","""0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK"""),
         ("chr17_41243153_122M_15G106_chr17_41243200_123M_123","41243200","123M","TTTATAAAATGAAACCAGAAGTAAGTCCACCAGTAATTAGGATGTTAAAGCTCATTCAGTCAAAGATGACGTCCTAGCTGTGTGAAGGACTTTTTTCTATGAAAAGCACCTTAGGAGGAACAT","""0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK"""),
        ),
        (
         None,
         ("chr17_41243153_122M_27G94_chr17_41243200_123M_123","41243200","123M","TTTATAAAATGAAACCAGAAGTAAGTCCACCAGTAATTAGGATGTTAAAGCTCATTCAGTCAAAGATGACGTCCTAGCTGTGTGAAGGACTTTTTTCTATGAAAAGCACCTTAGGAGGAACAT","""0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK""")),
        (
         ("chr13_32900227_120M_79G40_chr13_32900272_123M_34G88","32900227","120M","TTTTATTTTAGTCCTGTTGTTCTACAATGTACACATGTAACACCACAAAGAGATAAGTCAGGTATGATTAAAAACAATGGTTTTTATTCTTAGAATACTAGAAATGTTAATAAAAATAAA","""0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK"""),
         ("chr13_32900227_120M_79G40_chr13_32900272_123M_34G88","32900272","123M","CAAAGAGATAAGTCAGGTATGATTAAAAACAATGGTTTTTATTCTTAGAATACTAGAAATGTTAATAAAAATAAAACTTAACAATTTTCCCCTTTTTTTACCCCCAGTGGTATGTGGGAGTTT","""0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK""")
        ),
        (
         None,
         ("chr13_32900227_120M_120_chr13_32900272_123M_122C","32900272","123M","CAAAGAGATAAGTCAGGTATGATTAAAAACAATGCTTTTTATTCTTAGAATACTAGAAATGTTAATAAAAATAAAACTTAACAATTTTCCCCTTTTTTTACCCCCAGTGGTATGTGGGAGTTC","""0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK"""),
        ),
        ]







class test_insert_variants(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_overlaps(self):
        self.assertEqual(overlaps(10,1,10,10),True)
        self.assertEqual(overlaps(9,1,10,10),False)
        self.assertEqual(overlaps(10,5,10,20),True)
        self.assertEqual(overlaps(21,1,10,20),False)
        self.assertEqual(overlaps(13,2,10,20),True)
        self.assertEqual(overlaps(5,25,10,20),True)
        self.assertEqual(overlaps('41243153','122','41243153','41243153'),True)
        pass

    def test_get_readpairs_overlapping_position(self):
        self.assertEqual(get_readpairs_overlapping_position(io.StringIO(SAMFILE),*split_location_string("chr17:41243153-41243153")),READS[0:1])
        self.assertEqual(get_readpairs_overlapping_position(io.StringIO(SAMFILE),*split_location_string("chr17:41243322-41243322")),READS[0:2])
        self.assertEqual(get_readpairs_overlapping_position(io.StringIO(SAMFILE),*split_location_string("chr17:41243323-41243323")),[])
        self.assertEqual(get_readpairs_overlapping_position(io.StringIO(SAMFILE),*split_location_string("chr17:41243160-41243161")),READS[0:1])
        self.assertEqual(get_readpairs_overlapping_position(io.StringIO(SAMFILE),*split_location_string("chr17:41243250-41243250")),READS[0:2])
        self.assertEqual(get_readpairs_overlapping_position(io.StringIO(SAMFILE),*split_location_string("chr13:32900230-32900231")),READS[2:3])
        self.assertEqual(get_readpairs_overlapping_position(io.StringIO(SAMFILE),*split_location_string("chr13:32900351-32900351")),READS[2:])
        self.assertEqual(get_readpairs_overlapping_position(io.StringIO(SAMFILE),*split_location_string("chr17:41243152-41243152")),[])
        self.assertEqual(get_readpairs_overlapping_position(io.StringIO(SAMFILE),*split_location_string("chrX:100000000-100000001")),[])
        pass
    
    def test_get_readpairs_not_overlapping_position(self):
        self.assertEqual(list(get_reads_not_overlapping_position(io.StringIO(SAMFILE),*split_location_string("chr17:41243160-41243161"))),
                        [(READS[1][1][0],reverse_complement(READS[1][1][3]),READS[1][1][4][::-1]),
                        (READS[2][0][0],READS[2][0][3],READS[2][0][4]),
                        (READS[2][1][0],reverse_complement(READS[2][1][3]),READS[2][1][4][::-1]),
                        (READS[3][1][0],reverse_complement(READS[3][1][3]),READS[3][1][4][::-1]),
                        ])
        self.assertEqual(list(get_reads_not_overlapping_position(io.StringIO(SAMFILE),*split_location_string("chr17:41243250-41243250"))),
                        [(READS[2][0][0],READS[2][0][3],READS[2][0][4]),
                        (READS[2][1][0],reverse_complement(READS[2][1][3]),READS[2][1][4][::-1]),
                        (READS[3][1][0],reverse_complement(READS[3][1][3]),READS[3][1][4][::-1]),
                        ])
        self.assertEqual(list(get_reads_not_overlapping_position(io.StringIO(SAMFILE),*split_location_string("chr13:32900230-32900231"))),
                        [(READS[0][0][0],READS[0][0][3],READS[0][0][4]),
                        (READS[0][1][0],reverse_complement(READS[0][1][3]),READS[0][1][4][::-1]),
                        (READS[1][1][0],reverse_complement(READS[1][1][3]),READS[1][1][4][::-1]),
                        (READS[3][1][0],reverse_complement(READS[3][1][3]),READS[3][1][4][::-1]),
                        ])
        self.assertEqual(list(get_reads_not_overlapping_position(io.StringIO(SAMFILE),*split_location_string("chr13:32900351-32900351"))),
                        [(READS[0][0][0],READS[0][0][3],READS[0][0][4]),
                        (READS[0][1][0],reverse_complement(READS[0][1][3]),READS[0][1][4][::-1]),
                        (READS[1][1][0],reverse_complement(READS[1][1][3]),READS[1][1][4][::-1]),
                        ])
        self.assertEqual(list(get_reads_not_overlapping_position(io.StringIO(SAMFILE),*split_location_string("chrX:100000000-100000001"))),
                        [(READS[0][0][0],READS[0][0][3],READS[0][0][4]),
                        (READS[0][1][0],reverse_complement(READS[0][1][3]),READS[0][1][4][::-1]),
                        (READS[1][1][0],reverse_complement(READS[1][1][3]),READS[1][1][4][::-1]),
                        (READS[2][0][0],READS[2][0][3],READS[2][0][4]),
                        (READS[2][1][0],reverse_complement(READS[2][1][3]),READS[2][1][4][::-1]),
                        (READS[3][1][0],reverse_complement(READS[3][1][3]),READS[3][1][4][::-1]),
                        ])
        pass
    
    def test_point_mutate_read(self):
        #point_mutate_read(name,pos,cigar,seq,qual,one_based_variant_site,variant_base)
        self.assertEqual(point_mutate_read(*READS[0][0], one_based_variant_site=41243163,variant_base='g'),
                        ('chr17_41243153_122M_15G106_chr17_41243200_123M_123',
                        'TCACACAAAAgGATTGAATTCCTTGCTTTGGGACACCTGGATTTGCTTTTATAAAATGAAACCAGAAGTAAGTCCACCAGTAATTAGGATGTTAAAGCTCATTCAGTCAAAGATGACGTCCT',
                        '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJK')
                        )
        self.assertEqual(point_mutate_read('','10','12M','GCTCGCTCGCTC','I'*12, one_based_variant_site=10,variant_base='A'),
                        ('', 'ACTCGCTCGCTC', 'IIIIIIIIIIII'))
        self.assertEqual(point_mutate_read('','10','12M','GCTCGCTCGCTC','I'*12, one_based_variant_site=21,variant_base='T'),
                        ('', 'GCTCGCTCGCTT', 'IIIIIIIIIIII'))
        self.assertEqual(point_mutate_read('','12','2S8M','GCTCGCTCGCTC','I'*12, one_based_variant_site=15,variant_base='A'),
                        ('', 'GCTCGATCGCTC', 'IIIIIIIIIIII'))
        self.assertEqual(point_mutate_read('','10','2M1D10M','GCTCGCTCGCTC','I'*12, one_based_variant_site=15,variant_base='A'),
                        ('', 'GCTCACTCGCTC', 'IIIIIIIIIIII'))
        self.assertEqual(point_mutate_read('','10','2M1I9M','GCTCGCTCGCTC','I'*12, one_based_variant_site=15,variant_base='A'),
                        ('', 'GCTCGCACGCTC', 'IIIIIIIIIIII'))
        self.assertEqual(point_mutate_read('','10','4M1I7M','GCTCGCTCGCTC','I'*12, one_based_variant_site=15,variant_base='A'),
                        ('', 'GCTCGCACGCTC', 'IIIIIIIIIIII'))
        self.assertEqual(point_mutate_read('','10','3M3I6M','GCTCGCTCGCTC','I'*12, one_based_variant_site=15,variant_base='A'),
                        ('', 'GCTCGCTCACTC', 'IIIIIIIIIIII'))
        self.assertRaises(IndexError,point_mutate_read,'','10','12M','GCTCGCTCGCTC','I'*12, one_based_variant_site=22,variant_base='T')
        self.assertRaises(IndexError,point_mutate_read,'','10','6M1I5M','GCTCGCTCGCTC','I'*12, one_based_variant_site=21,variant_base='T')
        self.assertRaises(IndexError,point_mutate_read,'','10','2S10M','GCTCGCTCGCTC','I'*12, one_based_variant_site=20,variant_base='T')
        self.assertEqual(point_mutate_read('','10','4M1D8M','GCTCGCTCGCTC','I'*12, one_based_variant_site=15,variant_base='A'),
                        ('', 'GCTCACTCGCTC', 'IIIIIIIIIIII')) #corner case - not currently correct behaviour
        #self.assertEqual(point_mutate_read('','10','4M1D8M','GCTCGCTCGCTC','I'*12, one_based_variant_site=15,variant_base='A'),
        #                ('', 'GCTCAGCTCGCTC', 'IIIIIIIIIIIII')) #should put mismatch to right of deletion and replace a deletion
        #self.assertEqual(point_mutate_read('','10','4M3D8M','GCTCGCTCGCTC','I'*12, one_based_variant_site=15,variant_base='A'),
        #                ('', 'GCTCAGCTCGCTC', 'IIIIIIIIIIIII')) #string is the same no matter size or location in deletion.
    
    def test_deletion_mutate_read(self):
        self.assertEqual(deletion_mutate_read('','10','12M','GCTCGCTCGCTC','0123456789!@', one_based_variant_site=15,length=1),
                        ('', 'GCTCGTCGCTC', '012346789!@'))
        self.assertEqual(deletion_mutate_read('','10','12M','GCTCGCTCGCTC','0123456789!@', one_based_variant_site=15,length=3),
                        ('', 'GCTCGGCTC', '0123489!@'))
        self.assertEqual(deletion_mutate_read('','10','3M1D9M','GCTCGCACGCTC','0123456789!@', one_based_variant_site=12,length=5),
                        ('', 'GCACGCTC', '016789!@'))
        self.assertEqual(deletion_mutate_read('','10','3M1I8M','GCTCGCTCGCTC','0123456789!@', one_based_variant_site=12,length=5),
                        ('', 'GCGCTC', '0189!@'))
        self.assertRaises(IndexError,deletion_mutate_read,'','10','2S10M','GCTCGCTCGCTC','I'*12, one_based_variant_site=25,length=1)
           
if __name__ == '__main__':
    unittest.main()