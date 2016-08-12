#!/usr/bin/env python
# encoding: utf-8
"""
test_truseq_manifest.py

Created by Matthew Wakefield on 2013-09-23.
Copyright (c) 2013-2015 The Walter and Eliza Hall Institute. All rights reserved.
"""

from __future__ import print_function
import unittest
import io
from tempfile import NamedTemporaryFile
import hashlib
from ambivert.truseq_manifest import *
from pkg_resources import resource_stream

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2013-2015,  Matthew Wakefield and The University of Melbourne"
__credits__ = ["Matthew Wakefield",]
__license__ = "GPLv3"
__version__ = "0.5.1"
__maintainer__ = "Matthew Wakefield"
__email__ = "matthew.wakefield@unimelb.edu.au"
__status__ = "Development"

MANIFEST = u"""[Header]
Customer Name	Joseph Bloggs
Product Type	15026626
Date Manufactured	25/12/2012
Lot	WO0001234567890
DesignStudio ID	99999
Target Plexity	2

[Probes]
Target Region Name	Target Region ID	Target ID	Species	Build ID	Chromosome	Start Position	End Position	Submitted Target Region Strand	ULSO Sequence	ULSO Genomic Hits	DLSO Sequence	DLSO Genomic Hits	Probe Strand	Designer	Design Score	Expected Amplifed Region Size	SNP Masking	Labels
BRCA1_Exon9_UserDefined (9825051)	9825051	BRCA1_Exon9_UserDefined (9825051)_7473609	Homo sapiens	hg19	chr17	41243252	41247077	-	AAAGGAACTGCTTCTTAAACTTGAAAC	0	CATGAAAGTAAATCCAGTCCTGCCAATG	0	-	ILLUMINA	0.937	225	false	
BRCA1_Exon9_UserDefined (9825051)	9825051	BRCA1_Exon9_UserDefined (9825051)_7473610	Homo sapiens	hg19	chr17	41243252	41247077	+	GACGTCCTAGCTGTGTGAAGGACTTTT	0	TCCATGCTTTGCTCTTCTTGATTATTTTC	0	+	ILLUMINA	0.937	225	false	
BRCA2_exon22_23_UserDefined (9825052)	9825052	BRCA2_exon22_23_UserDefined (9825052)_7473650	Homo sapiens	hg19	chr13	32900038	32900619	-	GTACATTGTAGAACAACAGGACTAAAATAA	0	GGAAACCGTGTATTTTAAACTCAAAAATT	0	-	ILLUMINA	0.791	272	false	
BRCA2_exon22_23_UserDefined (9825052)	9825052	BRCA2_exon22_23_UserDefined (9825052)_7473651	Homo sapiens	hg19	chr13	32900038	32900619	+	TTTTTAAAATAACCTAAGGGATTTGCTTTG	0	GTTTCATACACCAAAGTTTGTGAAGGT	0	+	ILLUMINA	0.791	225	false	
[Targets]
TargetA	TargetB	Target Number	Chromosome	Start Position	End Position	Probe Strand	Sequence	Species	Build ID
BRCA1_Exon9_UserDefined (9825051)_7473609	BRCA1_Exon9_UserDefined (9825051)_7473609	1	chr17	41243125	41243349	+	TCACACAAAATGATTAAATTCCTTGCTTTGGGACACCTGGATTTGCTTTTATAAAATGAAACCAGAAGTAAGTCCACCAGTAATTAGGATGTTAAAGCTCATTCAGTCAAAGATGACGTCCTAGCTGTGTGAAGGACTTTTTTCTATGAAAAGCACCTTAGGAGGAACAT	Homo sapiens	hg19
BRCA1_Exon9_UserDefined (9825051)_7473610	BRCA1_Exon9_UserDefined (9825051)_7473610	1	chr17	41243267	41243491	-	TTCAAACTTAGGTATTGGAACCAGGTTTTTGTGTTTGCCCCAGTCTATTTATAGAAGTGAGCTAAATGTTTATGCTTTTGGGGAGCACATTTTACAAATTTCCAAGTATAGTTAAAGGAACTGCTTCTTAAACTTGAAACATGTTCCTCCTAAGGTGCTTTTCATAGAA	Homo sapiens	hg19
BRCA2_exon22_23_UserDefined (9825052)_7473650	BRCA2_exon22_23_UserDefined (9825052)_7473650	1	chr13	32899987	32900258	+	AGCAGCTGAAATTTGTGAGTACATATGTGTTGGCATTTTAAACATCACTTGATGATTATTTAATGCTTCATGAGAGATTTACTTTTTAAAATGTAATATAAAATATCTAAAAGTAGTATTCCAACAATTTATATGAATGAGAATCTTCTTTTAAAAATAAGATAAACTAGTTTTTGCCAGTTTTTTAAAATAACCTAAGGGATTTGCTTTGTT	Homo sapiens	hg19
BRCA2_exon22_23_UserDefined (9825052)_7473651	BRCA2_exon22_23_UserDefined (9825052)_7473651	1	chr13	32900197	32900421	-	AAACTCCCACATACCACTGGGGGTAAAAAAAGGGGAAAATTGTTAAGTTTTATTTTTATTAACATTTCTAGTATTCTAAGAATAAAAAGCATTGTTTTTAATCATACCTGACTTATCTCTTTGTGGTGTTACATGTGTACATTGTAGAACAACAGGACTAAAATAAAA	Homo sapiens	hg19

"""

def md5(data):
    #python3 compatibility
    return hashlib.md5(bytes(data,'ascii'))

class test_truseq_manifest(unittest.TestCase):
    def setUp(self):
        manifest = io.TextIOWrapper(resource_stream(__name__, 'data/testdatamanifest.txt'))
        self.header, self.probes, self.targets = parse_truseq_manifest(manifest)
        pass
    
    def test_parse_truseq_manifest(self):        
        #print(self.header, self.probes, self.targets)
        self.assertEqual(md5(str(sorted(self.header.items()))).hexdigest(), '9442a2a72901b006f8ea16640cf5d418')
        self.assertEqual(md5(str(self.probes)).hexdigest(), '7b3526b5867201bf8e2b848235a39343')
        self.assertEqual(md5(str(self.targets)).hexdigest(), '6a96cfc57eeec0c981302c7770a525fb')
        pass
        
    def test_make_probes(self):
        #without adaptors
        outfile = NamedTemporaryFile(delete=False,mode='wt')
        make_probes(self.header, self.probes, self.targets, output=outfile)
        outfile = open(outfile.name)
        #print(outfile.read())
        #outfile.seek(0)
        self.assertEqual(md5(outfile.read()).hexdigest(),'a9d739c4ba4e7ecb1c2b10b69171e8dd')
        outfile.close()
        os.unlink(outfile.name)

        #with adaptors
        outfile = NamedTemporaryFile(delete=False,mode='wt')
        make_probes(self.header, self.probes, self.targets, adaptors=True, output=outfile)
        outfile = open(outfile.name)
        #print(outfile.read())
        #outfile.seek(0)
        self.assertEqual(md5(outfile.read()).hexdigest(),'d4f7115e442204fd668b0980fb23d305')
        outfile.close()
        os.unlink(outfile.name)
        pass
        
    def test_make_fasta(self):
        # with_probes=False, softmask_probes=False, all_plus=True
        outfile = NamedTemporaryFile(delete=False,mode='wt')
        make_fasta(self.header, self.probes, self.targets, output=outfile, all_plus=True)
        outfile = open(outfile.name)
        #print(outfile.read())
        #outfile.seek(0)
        self.assertEqual(md5(outfile.read()).hexdigest(),'721ea92b3bc94910f0f1396077694161')
        outfile.close()
        os.unlink(outfile.name)
        
        # with_probes=False, softmask_probes=False, all_plus=False
        outfile = NamedTemporaryFile(delete=False,mode='wt')
        make_fasta(self.header, self.probes, self.targets, output=outfile, all_plus=False)
        outfile = open(outfile.name)
        #print(outfile.read())
        #outfile.seek(0)
        self.assertEqual(md5(outfile.read()).hexdigest(),'4d5112f709c03dc63019cfa7d1f07941')
        outfile.close()
        os.unlink(outfile.name)
        
        # with_probes=True, softmask_probes=False, all_plus=True
        outfile = NamedTemporaryFile(delete=False,mode='wt')
        make_fasta(self.header, self.probes, self.targets, output=outfile, with_probes=True, softmask_probes=False, all_plus=True)
        outfile = open(outfile.name)
        #print(outfile.read())
        #outfile.seek(0)
        self.assertEqual(md5(outfile.read()).hexdigest(),'cfcec95a220139e574132b94b690aa7e')
        outfile.close()
        os.unlink(outfile.name)
        
        # with_probes=True, softmask_probes=True, all_plus=True
        outfile = NamedTemporaryFile(delete=False,mode='wt')
        make_fasta(self.header, self.probes, self.targets, output=outfile, with_probes=True, softmask_probes=True, all_plus=True)
        outfile = open(outfile.name)
        #print(outfile.read())
        #outfile.seek(0)
        self.assertEqual(md5(outfile.read()).hexdigest(),'1a2ea63012b713fc9f1dce9eddf36a9c')
        outfile.close()
        os.unlink(outfile.name)
        
        # with_probes=True, softmask_probes=True, all_plus=False
        outfile = NamedTemporaryFile(delete=False,mode='wt')
        make_fasta(self.header, self.probes, self.targets, output=outfile, with_probes=True, softmask_probes=True, all_plus=False)
        outfile = open(outfile.name)
        #print(outfile.read())
        #outfile.seek(0)
        self.assertEqual(md5(outfile.read()).hexdigest(),'c09802e71fb028768388b6574539e606')
        outfile.close()
        os.unlink(outfile.name)

    
if __name__ == '__main__':
    unittest.main()