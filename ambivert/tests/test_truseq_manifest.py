#!/usr/bin/env python
# encoding: utf-8
"""
test_truseq_manifest.py

Created by Matthew Wakefield on 2013-09-23.
Copyright (c) 2013 The Walter and Eliza Hall Institute. All rights reserved.
"""

from __future__ import print_function
import unittest
import io
from tempfile import NamedTemporaryFile
from hashlib import md5
from ambivert.truseq_manifest import *
from pkg_resources import resource_stream

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

class test_truseq_manifest(unittest.TestCase):
    def setUp(self):
        self.header, self.probes, self.targets = parse_truseq_manifest(io.StringIO(MANIFEST))
        pass
    
    def test_parse_truseq_manifest(self):        
        self.assertEqual(md5(str(self.header)).hexdigest(), '8fe60036dd4cc7d22c731338678a94e0')
        self.assertEqual(md5(str(self.probes)).hexdigest(), '34c79046f6822cd1ff585eee928e9088')
        self.assertEqual(md5(str(self.targets)).hexdigest(), 'fbab15899042e8e751bdad0d3997f3af')
        #print(self.header, self.probes, self.targets)
        
        manifest = resource_stream(__name__, 'testdatamanifest.txt')
        self.header, self.probes, self.targets = parse_truseq_manifest(manifest)
        self.assertEqual(md5(str(self.header)).hexdigest(), '2f981c324cea3f412abee2c4c306d74d')
        self.assertEqual(md5(str(self.probes)).hexdigest(), '7b3526b5867201bf8e2b848235a39343')
        self.assertEqual(md5(str(self.targets)).hexdigest(), '6a96cfc57eeec0c981302c7770a525fb')
        
        pass
        
    def test_make_probes(self):
        #without adaptors
        outfile = NamedTemporaryFile(delete=False)
        make_probes(self.header, self.probes, self.targets, output=outfile)
        outfile = open(outfile.name)
        #print(outfile.read())
        #outfile.seek(0)
        self.assertEqual(md5(outfile.read()).hexdigest(),'0958aee1ba65c510320ea41fdde12359')
        outfile.close()
        os.unlink(outfile.name)

        #with adaptors
        outfile = NamedTemporaryFile(delete=False)
        make_probes(self.header, self.probes, self.targets, adaptors=True, output=outfile)
        outfile = open(outfile.name)
        #print(outfile.read())
        #outfile.seek(0)
        self.assertEqual(md5(outfile.read()).hexdigest(),'e43c1292d8792bbab8f4b99a9b55db59')
        outfile.close()
        os.unlink(outfile.name)
        pass
        
    def test_make_fasta(self):
        # with_probes=False, softmask_probes=False, all_plus=True
        outfile = NamedTemporaryFile(delete=False)
        make_fasta(self.header, self.probes, self.targets, output=outfile, all_plus=True)
        outfile = open(outfile.name)
        #print(outfile.read())
        #outfile.seek(0)
        self.assertEqual(md5(outfile.read()).hexdigest(),'b2dde8e4110ccb20bbb2c474c77c4712')
        outfile.close()
        os.unlink(outfile.name)
        
        # with_probes=False, softmask_probes=False, all_plus=False
        outfile = NamedTemporaryFile(delete=False)
        make_fasta(self.header, self.probes, self.targets, output=outfile, all_plus=False)
        outfile = open(outfile.name)
        #print(outfile.read())
        #outfile.seek(0)
        self.assertEqual(md5(outfile.read()).hexdigest(),'d253080b5236b77a4bafce88f4b7e329')
        outfile.close()
        os.unlink(outfile.name)
        
        # with_probes=True, softmask_probes=False, all_plus=True
        outfile = NamedTemporaryFile(delete=False)
        make_fasta(self.header, self.probes, self.targets, output=outfile, with_probes=True, softmask_probes=False, all_plus=True)
        outfile = open(outfile.name)
        #print(outfile.read())
        #outfile.seek(0)
        self.assertEqual(md5(outfile.read()).hexdigest(),'6742cc3a47e4cb1e29eec03244e8028b')
        outfile.close()
        os.unlink(outfile.name)
        
        # with_probes=True, softmask_probes=True, all_plus=True
        outfile = NamedTemporaryFile(delete=False)
        make_fasta(self.header, self.probes, self.targets, output=outfile, with_probes=True, softmask_probes=True, all_plus=True)
        outfile = open(outfile.name)
        #print(outfile.read())
        #outfile.seek(0)
        self.assertEqual(md5(outfile.read()).hexdigest(),'1c1414efef85de5c96be8965704f658a')
        outfile.close()
        os.unlink(outfile.name)
        
        # with_probes=True, softmask_probes=True, all_plus=False
        outfile = NamedTemporaryFile(delete=False)
        make_fasta(self.header, self.probes, self.targets, output=outfile, with_probes=True, softmask_probes=True, all_plus=False)
        outfile = open(outfile.name)
        #print(outfile.read())
        #outfile.seek(0)
        self.assertEqual(md5(outfile.read()).hexdigest(),'1bfd9f13e175f65380d8b2c6f7d17a74')
        outfile.close()
        os.unlink(outfile.name)

    
if __name__ == '__main__':
    unittest.main()