import unittest

import rpy2.robjects as robjects
import bioc.biobase as bb

ridentical = robjects.r["identical"]

class AnnotatedDataFrameTestCase(unittest.TestCase):

    def setUp(self):
        self.adf = robjects.r["new"]("AnnotatedDataFrame")

    def test_pData(self):
        adf = self.adf
        pdata = adf.pdata
        self.assertTrue(robjects.r["is.data.frame"](pdata)[0])
        adf.pdata = pdata
        rpdata = robjects.r["pData"](adf)
        self.assertTrue(adf.pdata.rsame(rpdata))

    def test_combine(self):
        adf = self.adf
        cadf = adf.combine(adf)
        self.assertFalse(cadf.pdata.rsame(adf))
        self.assertTrue(ridentical(cadf, adf)[0])

    def test_featurenames(self):
        adf = self.adf
        featurenames = adf.featurenames
        self.assertEquals(len(featurenames), 
                          robjects.r["nrow"](adf.pdata)[0]) 

    def test_new(self):
        adf = bb.AnnotatedDataFrame.new()
        self.assertTrue(False) # test not properly implemented (yet)
    
class ExpressionSetTestCase(unittest.TestCase):

    def setUp(self):
        self.eset = robjects.r["new"]("ExpressionSet")

    def test_exprs(self):
        eset = self.eset
        exprs = eset.exprs
        self.assertTrue(robjects.r["is.matrix"](exprs)[0])
        eset.exprs = exprs
        rexprs = robjects.r["exprs"](eset)
        self.assertTrue(eset.exprs.rsame(rexprs))

    def test_new(self):
        
        eset = bb.ExpressionSet.new()

def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(BiobaseTestCase)
    suite.addTest(unittest.TestLoader().loadTestsFromTestCase(ExpressionSetTestCase))
    return suite

if __name__ == '__main__':
     unittest.main()
