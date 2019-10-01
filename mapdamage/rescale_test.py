import unittest
import optparse
from . import rescale
import pysam
import filecmp

def mock_options(filename,rescale_out,folder):
    """Make the options object with nice values for testing"""
    return optparse.Values({
        "filename":filename,
        "rescale_out":rescale_out,
        "verbose":True,
        "folder":folder,
        "quiet":True
        })


class testPairedFile(unittest.TestCase):
    """Tests if rescaling of a paired end file"""
    def test_paired_end_file(self):
        """Test, paired end SAM file"""
        # here is the example
        # >ref
        # 5p CGA AAA CGA 3p
        # 3p GCT TTT GCT 5p
        # >r001/1
        # CGA
        # >r001/2
        # GCA
        # The sam file looks like this
        #@HD VN:1.5 SO:coordinate
        #@SQ SN:ref LN:9
        #r001 163 ref 1 30 3M = 7 9 CGA III MR:f:0 #the normal ones
        #r001 83 ref 7 30 3M = 1 -9 TCG III MR:f:0
        #r002 163 ref 1 30 3M = 7 9 TGA III MR:f:0.9 #With one dam subs
        #r002 83 ref 7 30 3M = 1 -9 CAA III MR:f:0.009
        #r003 83 ref 1 30 3M = 7 9 TGA III #The reverse complement, should not rescale (thus no flags)
        #r003 163 ref 7 30 3M = 1 -9 CAA III
        # 
        #hand calculated the correct rescaled sam file in pe_rescaled_correct.sam
        options = mock_options("rescale_test/pe_test/pe.sam","rescale_test/pe_test/pe_rescaled.sam","rescale_test/pe_test/pe_output/")
        ref = pysam.Fastafile("rescale_test/pe_test/ref.fa")
        rescale.rescale_qual(ref,options,debug=True)
        self.assertTrue(filecmp.cmp("rescale_test/pe_test/pe_rescaled.sam","rescale_test/pe_test/pe_rescaled_correct.sam"))

class testCases(unittest.TestCase):
    """Various cases that failed"""
    def test_longalignshortread(self):
        """Check if fails on an aligment longer than the read"""
        prefix="rescale_test/long_align/"
        options = mock_options(prefix+"pe.sam",prefix+"pe_out.sam",prefix+"pe_output/")
        ref = pysam.Fastafile("rescale_test/pe_test/ref.fa")
        rescale.rescale_qual(ref,options,debug=True)
        #Should run without an error

if  __name__=='__main__':
    unittest.main()
