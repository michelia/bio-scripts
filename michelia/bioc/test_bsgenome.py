import unittest

import rpy2.robjects as robjects
import bioc.bsgenome as bg
from test_utils import RS4Test_Type

ridentical = robjects.r["identical"]


class GenomeDescriptionTestCase(unittest.TestCase):
    __metaclass__ = RS4Test_Type
    __targetclass__ = bg.GenomeDescription
    obj = None
    
    def setUp(self):
        robjects.r('library(BSgenome.Celegans.UCSC.ce2)')
        self.obj = robjects.r['Celegans']



def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(GenomeDescriptionTestCase)
    #suite.addTest(unittest.TestLoader().loadTestsFromTestCase(FooTestCase))
    return suite

if __name__ == '__main__':
     unittest.main()
