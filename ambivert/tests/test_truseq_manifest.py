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
        self.assertEqual(md5(outfile.read()).hexdigest(),'8fe08cd60ed7c306f7fb63218dba4a84')
        outfile.close()
        os.unlink(outfile.name)
        
        # with_probes=False, softmask_probes=False, all_plus=False
        outfile = NamedTemporaryFile(delete=False)
        make_fasta(self.header, self.probes, self.targets, output=outfile, all_plus=False)
        outfile = open(outfile.name)
        #print(outfile.read())
        #outfile.seek(0)
        self.assertEqual(md5(outfile.read()).hexdigest(),'96c2e8ea6fe8f0a357f5365a80af4937')
        outfile.close()
        os.unlink(outfile.name)
        
        # with_probes=True, softmask_probes=False, all_plus=True
        outfile = NamedTemporaryFile(delete=False)
        make_fasta(self.header, self.probes, self.targets, output=outfile, with_probes=True, softmask_probes=False, all_plus=True)
        outfile = open(outfile.name)
        #print(outfile.read())
        #outfile.seek(0)
        self.assertEqual(md5(outfile.read()).hexdigest(),'bfdae1f47a817b667b2a34e6e1c32fc0')
        outfile.close()
        os.unlink(outfile.name)
        
        # with_probes=True, softmask_probes=True, all_plus=True
        outfile = NamedTemporaryFile(delete=False)
        make_fasta(self.header, self.probes, self.targets, output=outfile, with_probes=True, softmask_probes=True, all_plus=True)
        outfile = open(outfile.name)
        #print(outfile.read())
        #outfile.seek(0)
        self.assertEqual(md5(outfile.read()).hexdigest(),'e455c9435915c8dcf9b24acb70e38feb')
        outfile.close()
        os.unlink(outfile.name)
        
        # with_probes=True, softmask_probes=True, all_plus=False
        outfile = NamedTemporaryFile(delete=False)
        make_fasta(self.header, self.probes, self.targets, output=outfile, with_probes=True, softmask_probes=True, all_plus=False)
        outfile = open(outfile.name)
        #print(outfile.read())
        #outfile.seek(0)
        self.assertEqual(md5(outfile.read()).hexdigest(),'486df188e89ed92f3277f255d83d91bc')
        outfile.close()
        os.unlink(outfile.name)

    
if __name__ == '__main__':
    unittest.main()