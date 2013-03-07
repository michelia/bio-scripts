#    map bioconductor classes to Python classes 
#    Copyright (C) 2009-2010  Laurent Gautier
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
A module to model the ShortRead library in Bioconductor

Copyright 2009 - Laurent Gautier

"""


import rpy2.robjects.methods
import rpy2.robjects as robjects
import rpy2.robjects.packages
from rpy2.robjects.packages import importr
conversion = robjects.conversion

getmethod = robjects.baseenv.get("getMethod")

import bioc.biobase

__rname__ = 'ShortRead'

__rpackage__ = importr(__rname__)

StrVector = robjects.StrVector

class _ShortReadBase(robjects.methods.RS4):
    pass

class ShortRead(_ShortReadBase):
    __metaclass__ = rpy2.robjects.methods.RS4_Type

    __accessors__ = (('length', 'ShortRead', '__len__', True, None),
                     ('width', 'ShortRead', None, True, None),
                     ('srorder', 'ShortRead', None, False, None),
                     ('srrank', 'ShortRead', None, False, None),
                     ('srsort', 'ShortRead', None, False, None),
                     ('srduplicated', 'ShortRead', None, False, None),
                     ('detail', 'ShortRead', None, False, None),
                     ('clean', 'ShortRead', None, False, None),
                     ('alphabetByCycle', 'ShortRead', 'alphabetbycycle', False, None)
                     )
    _tables = getmethod("tables", 
                        signature = StrVector(["ShortRead", ]), 
                        where="package:ShortRead")
    def tables(self, n = 50):
        res = self._tables(self, n = n)
        return res

    _narrow = getmethod("narrow", 
                        signature = StrVector(["ShortRead", ]), 
                        where="package:ShortRead")
    def narrow(self, start = None, end = None, width = None, use_names = True):
        res = self._narrow(self, start = start, end = end, 
                           width = width, use_names = use_names)
        return res

class ShortReadQ(ShortRead):
    pass

class AlignedRead(ShortReadQ):
    _coverage_coords = robjects.StrVector(('leftmost', 'fiveprime'))

    _coverage = getmethod("coverage", 
                          signature = StrVector(["AlignedRead", ]), 
                          where = "package:ShortRead")
    _extract = getmethod("[",
                         signature = StrVector(["AlignedRead", 
                                                "ANY", "missing"])
                         )

    def coverage(self, 
                 start = robjects.rinterface.NA_Logical,
                 end = robjects.rinterface.NA_Logical,
                 coords = _coverage_coords, 
                 extend = 0):
        res = self._coverage(self, start=start, end=end, coords=coords, extend=extend)
        return res

    def rx(self, i, drop = True):
        res = self._extract(self, i)
        res = conversion.ri2py(res)
        return res

class Intensity(_ShortReadBase):
    # __metaclass__ = rpy2.robjects.methods.RS4_Type

    # __accessors__ = (('readInfo', 'ShortRead', 'readinfo', True, None),
    #                  ('intensity', 'ShortRead', None, True, None),
    #                  ('measurementError', 'ShortRead', 'measurementerror', True, None),
    #                  ('dim', 'ShortRead', None, True, None),
    #                  )
    pass

class IntensityInfo(_ShortReadBase):
    pass
class IntensityMeasure(_ShortReadBase):
    pass
class ArrayIntensity(IntensityMeasure):
    pass

class QualityScore(_ShortReadBase):
    __metaclass__ = rpy2.robjects.methods.RS4_Type

    __accessors__ = (('length', 'ShortRead', '__len__', True, None),
                     ('width', 'ShortRead', None, True, None),
                     ('detail', 'ShortRead', None, False, None),
                     )

class NumericQuality(QualityScore):
    pass
class MatrixQuality(QualityScore):
    pass
class FastqQuality(QualityScore):
    __metaclass__ = rpy2.robjects.methods.RS4_Type

    __accessors__ = (
                     ('alphabetByCycle', 'ShortRead', 'alphabetbycycle', False, None),
                     ('alphabetFrequency', 'ShortRead', 'alphabetfrequency', False, None),
                     ('srorder', 'ShortRead', None, False, None),
                     ('srrank', 'ShortRead', None, False, None),
                     ('srduplicated', 'ShortRead', None, False, None),
                     )

class IntegerQuality(NumericQuality):
    pass

class SFastqQuality(FastqQuality):
    pass

class ExperimentPath(_ShortReadBase):
    __metaclass__ = rpy2.robjects.methods.RS4_Type

    __accessors__ = (
                     ('detail', 'ShortRead', None, False, None),
                     )

class _Solexa(_ShortReadBase):
    pass
class SolexaSet(_Solexa):
    pass
class _Roche(_ShortReadBase):
    pass
class RocheSet(_Roche):
    pass

class SolexaPath(ExperimentPath, _Solexa):
    __metaclass__ = rpy2.robjects.methods.RS4_Type

    __accessors__ = (
                     ('detail', 'ShortRead', None, False, None),
                     ('SolexaSet', 'ShortRead', 'solexaset', False, None),
                     )

    _readintensities =  getmethod("readIntensities", 
                                  signature = StrVector(["SolexaPath", ]), 
                                  where="package:ShortRead")

    _readprb =  getmethod("readPrb", 
                          signature = StrVector(["SolexaPath", ]), 
                          where="package:ShortRead")

    _readfastq =  getmethod("readFastq", 
                            signature = StrVector(["SolexaPath", ]), 
                            where="package:ShortRead")

    _readbasequality =  getmethod("readBaseQuality", 
                                  signature = StrVector(["SolexaPath", ]), 
                                  where="package:ShortRead")

    _readqseq =  getmethod("readQseq", 
                           signature = StrVector(["SolexaPath", ]), 
                           where="package:ShortRead")

    _readaligned =  getmethod("readAligned", 
                              signature = StrVector(["SolexaPath", ]), 
                              where="package:ShortRead")

    def readintensities(self, **kwargs):
        res = self._readintensities(self, **kwargs)
        return res

    def readprb(self, **kwargs):
        res = self._readprb(self, **kwargs)
        return res
    
    def readfastq(self, **kwargs):
        res = self._readfastq(self, **kwargs)
        return res
    
    def readqseq(self, **kwargs):
        res = self._readqseq(self, **kwargs)
        return res
    
    def readaligned(self, **kwargs):
        res = self._readaligned(self, **kwargs)
        return res

class RochePath(ExperimentPath, _Roche):
    __metaclass__ = rpy2.robjects.methods.RS4_Type

    __accessors__ = (
                     ('detail', 'ShortRead', None, False, None),
                     ('read454', 'ShortRead', None, False, None),
                     ('runNames', 'ShortRead', 'runnames', False, None),
                     ('RocheSet', 'ShortRead', 'rocheset', False, None),
                     )

    _readfasta =  getmethod("readFasta", 
                            signature = StrVector(["RochePath", ]), 
                            where="package:ShortRead")

    _readqual =  getmethod("readQual",
                            signature = StrVector(["RochePath", ]), 
                            where="package:ShortRead")

    def readfasta(self, **kwargs):
        res = self._readfastq(self, **kwargs)
        return res
    def readqual(self, **kwargs):
        res = self._readqual(self, **kwargs)
        return res
    def readpath(self):
        res = self.do_slot("readPath")
        return res

class SRSet(_ShortReadBase):
    pass

class SolexaIntensity(Intensity):
    pass

class SolexaIntensityInfo(IntensityInfo):
    pass

class AlignedDataFrame(bioc.biobase.AnnotatedDataFrame):
    __append = getmethod("append",
                         signature = StrVector(["AlignedDataFrame", "AlignedDataFrame", "missing"]),
                         where = "package:ShortRead")

    def append(self, **kwargs):
        res = self.__append(self, **kwargs)
        return res


_shortread_dict = {
    '.ShortReadBase': _ShortReadBase,
    'ShortRead': ShortRead,
    'ShortReadQ': ShortReadQ,
    'AlignedRead': AlignedRead,
    'AlignedDataFrame': AlignedDataFrame,
    'Intensity': Intensity,
    'IntensityInfo': IntensityInfo,
    'IntensityMeasure': IntensityMeasure,
    'ArrayIntensity': ArrayIntensity,
    'QualityScore': QualityScore,
    'NumericQuality': NumericQuality,
    'MatrixQuality': MatrixQuality,
    'FastqQuality': FastqQuality,
    'IntegerQuality': IntegerQuality,
    'SFastqQuality': FastqQuality,
    'ExperimentPath': ExperimentPath,
    'SRSet': SRSet,
    'SolexaIntensity': SolexaIntensity,
    'SolexaIntensityInfo': SolexaIntensityInfo,
    'SolexaPath': SolexaPath,
    'SolexaSet': SolexaSet,
    'RochePath': RochePath,
    'RocheSet': RocheSet
    }

original_conversion = conversion.ri2py
def shortread_conversion(robj):

    pyobj = original_conversion(robj)
    if isinstance(pyobj, robjects.RS4):
        rclass = [x for x in pyobj.rclass][0]
        try:
            cls = _shortread_dict[rclass]
            pyobj = cls(pyobj)
        except KeyError, ke:
            pass

    return pyobj

conversion.ri2py = shortread_conversion

