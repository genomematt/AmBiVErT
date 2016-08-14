#!/usr/bin/env python
# encoding: utf-8
"""
test_insilicoPCR_to_amplicons.py

Copyright (c) 2016 Matthew Wakefield and The University of Melbourne. All rights reserved.
"""

import unittest, io, os, hashlib, logging
from tempfile import NamedTemporaryFile


__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2016,  Matthew Wakefield and The University of Melbourne"
__credits__ = ["Matthew Wakefield",]
__license__ = "GPLv3"
__version__ = "0.5.1"
__maintainer__ = "Matthew Wakefield"
__email__ = "matthew.wakefield@unimelb.edu.au"
__status__ = "Development"

PCRFASTA = u""">chr17:7674820-7675000 181bp GCCTCTGATTCCTCACTGAT TTAACCCCTCCTCCCAGAGA
GCCTCTGATTCCTCACTGATtgctcttaggtctggcccctcctcagcatc
ttatccgagtggaaggaaatttgcgtgtggagtatttggatgacagaaac
acttttcgacatagtgtggtggtgccctatgagccgcctgaggtctggtt
tgcaactggggTCTCTGGGAGGAGGGGTTAA

>chr22:33908517+33908966 450bp GGACAGATTGATGATGCATGAAATGGG CCCATGAGTGGCTCCTAAAGCAGCTGC
ttACAGATTGATGATGCATGAAATGGGgggtggccaggggtggggggtga
gactgcagagaaaggcagggctggttcataacaagctttgtgcgtcccaa
tatgacagctgaagttttccaggggctgatggtgagccagtgagggtaag
tacacagaacatcctagagaaaccctcattccttaaagattaaaaataaa
gacttgctgtctgtaagggattggattatcctatttgagaaattctgtta
tccagaatggcttaccccacaatgctgaaaagtgtgtaccgtaatctcaa
agcaagctcctcctcagacagagaaacaccagccgtcacaggaagcaaag
aaattggcttcacttttaaggtgaatccagaacccagatgtcagagctcc
aagcactttgctctcagctccacGCAGCTGCTTTAGGAGCCACTCATGaG
"""

def md5(data):
    #python3 compatibility
    return hashlib.md5(bytes(data,'ascii'))

class test_insilicoPCR_to_amplicons(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_insilicoPCR_to_amplicons(self):
        from ambivert.insilicoPCR_to_amplicons import make_fasta
        pcrfasta = io.StringIO(PCRFASTA)
        outfile = NamedTemporaryFile(delete=False,mode='wt')
        outfilename = outfile.name
        output = io.StringIO('')
        make_fasta(pcrfasta, outfile)
        #print(open(outfilename,mode='rt').read())
        outfile = open(outfilename,mode='rt')
        self.assertEqual(md5(outfile.read()).hexdigest(), '0332d90379d62155d9a0a82417f02a9f')
        outfile.close()
        os.unlink(outfile.name)
        pass
        
    
if __name__ == '__main__':
    #logging.disable(logging.CRITICAL)
    unittest.main()
    #logging.disable(logging.NOTSET)