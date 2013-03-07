import unittest

import rpy2.robjects as robjects
import bioc.biostrings as bs

from test_utils import RS4Test_Type


class XStringTestCase(unittest.TestCase):
    __metaclass__ = RS4Test_Type
    __targetclass__ = bs.XString
    obj = None

    def setUp(self):
        self.obj = bs.DNAString.new("ATTACG")

    def test_nchar(self):
        xs = self.obj
        nchar = xs.nchar
        self.assertEquals(1, len(nchar))
        self.assertEquals(6, nchar[0])

class DNAStringTestCase(unittest.TestCase):
    
    def test___init__(self):
        ds = bs.DNAString.new("ATTACG")
        self.assertEquals(6, ds.nchar[0])

class RNAStringTestCase(unittest.TestCase):
    
    def test___init__(self):
        ds = bs.RNAString.new("AUUACG")
        self.assertEquals(6, ds.nchar[0])

class AAStringTestCase(unittest.TestCase):
    
    def test___init__(self):
        ds = bs.AAString.new("AULVCG")
        self.assertEquals(6, ds.nchar[0])


class DNAStringSetTestCase(unittest.TestCase):
    
    def test___init__(self):
        strvect = robjects.StrVector(("ATTACG", "ATTAAT"))
        dss = bs.DNAStringSet.new(strvect)
        self.assertEquals(2, len(dss))


class AlignedXStringSet0TestCase(unittest.TestCase):
    __metaclass__ = RS4Test_Type
    __targetclass__ = bs.AlignedXStringSet0
    obj = None

    def setUp(self):
        self.obj = robjects.r('''
        local({
        pattern <- AAString("LAND")
        subject <- AAString("LEAVES")
        nw1 <- pairwiseAlignment(pattern, subject,
                                 substitutionMatrix = "BLOSUM50",
                                 gapOpening = -3, gapExtension = -1)
        pattern(nw1)
        })''')
    
    def test_nchar(self):
        xs = self.obj
        nchar = xs.nchar
        self.assertEquals(1, len(nchar))
        self.assertEquals(6, nchar[0])


class InDelTestCase(unittest.TestCase):
    __metaclass__ = RS4Test_Type
    __targetclass__ = bs.InDel
    obj = None

    def setUp(self):
        self.obj = robjects.r('''
        local({
        pattern <- AAString("LAND")
        subject <- AAString("LEAVES")
        nw1 <- pairwiseAlignment(pattern, subject,
                                 substitutionMatrix = "BLOSUM50",
                                 gapOpening = -3, gapExtension = -1)
        nindel(nw1)
        })''')
    

class PairwiseAlignedXStringSetTestCase(unittest.TestCase):
    __metaclass__ = RS4Test_Type
    __targetclass__ = bs.PairwiseAlignedXStringSet
    obj = None

    def setUp(self):
        self.obj = robjects.r('''
        local({
        pattern <- AAString("LAND")
        subject <- AAString("LEAVES")
        nw1 <- pairwiseAlignment(pattern, subject,
                                 substitutionMatrix = "BLOSUM50",
                                 gapOpening = -3, gapExtension = -1)
        nw1
        })''')
    

def suite():
    tl = unittest.TestLoader()
    suite = tl.loadTestsFromTestCase(XStringTestCase)
    suite.addTest(tl.loadTestsFromTestCase(DNAStringTestCase))
    suite.addTest(tl.loadTestsFromTestCase(RNAStringTestCase))
    suite.addTest(tl.loadTestsFromTestCase(AAStringTestCase))
    suite.addTest(tl.loadTestsFromTestCase(DNAStringSetTestCase))
    suite.addTest(tl.loadTestsFromTestCase(AlignedXStringSet0TestCase))
    suite.addTest(tl.loadTestsFromTestCase(InDelTestCase))
    suite.addTest(tl.loadTestsFromTestCase(PairwiseAlignedXStringSetTestCase))
    return suite

if __name__ == '__main__':
     unittest.main()
