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
A module to model the IRanges library in Bioconductor

Copyright 2009 - Laurent Gautier

"""



import rpy2.robjects as robjects
import rpy2.robjects.packages
conversion = robjects.conversion

getmethod = robjects.baseenv.get("getMethod")

__rname__ = 'IRanges'
rpy2.robjects.packages.quiet_require(__rname__)

iranges_env = robjects.baseenv['as.environment']('package:IRanges')
__rpackage__ = robjects.packages.SignatureTranslatedPackage(iranges_env, __rname__)

StrVector = robjects.StrVector

class Sequence(robjects.methods.RS4):
    """ Abstract class """

class Rle(Sequence):
    _length = getmethod("length", 
                        signature = StrVector(["Rle", ]), 
                        where="package:IRanges")
    def __len__(self):
        res = self._length(self)
        return res[0]



class AtomicList(Sequence):
    pass


class LogicalList(AtomicList):
    pass

class IntegerList(AtomicList):
    pass

class NumericList(AtomicList):
    pass

class ComplexList(AtomicList):
    pass

class CharacterList(AtomicList):
    pass

class RawList(AtomicList):
    pass


class CompressedAtomicList(AtomicList):
    pass

class SimpleAtomicList(AtomicList):
    pass


class Ranges(robjects.RS4):
    _length = getmethod("length", 
                        signature = StrVector(["Ranges", ]), 
                        where="package:IRanges")
    def __len__(self):
        res = self._length(self)
        return res[0]

acs = (
    ('start', None, True, None),
    ('width', None, True, None),
    ('end', None, True, None),
    ('mid', None, False, None),
    ('isEmpty', 'isempty', False, None),
    ('isDisjoint', 'isdisjoint', False, None),
    ('isNormal', 'isnormal', False, None),
    ('whichFirstNotNormal', 'whichfirst_notnormal', False, None),
    ('as.matrix', 'as_matrix', False, None),
    ('as.integer', 'as_integer', False, None),
    ('as.data.frame', 'as_DataFrame', False, None)
    )
robjects.methods.set_accessors(Ranges,
                               "Ranges", "IRanges",
                               acs)
del(acs)

class XRanges(Ranges):
    pass

class IntervalTree(XRanges):
    _length = getmethod("length", 
                        signature = StrVector(["IntervalTree", ]), 
                        where="package:IRanges")
    def __len__(self):
        res = self._length(self)
        return res[0]



class IRanges(Ranges):
    pass

class NormalIRanges(IRanges):
    pass
acs = (
    ('max', None, True, None),
    ('min', None, True, None)
    )
robjects.methods.set_accessors(NormalIRanges,
                               "NormalIRanges", "IRanges",
                               acs)
del(acs)


class SimpleList(Sequence):
    _length = getmethod("length", 
                        signature = StrVector(["SimpleList", ]),
                        where="package:IRanges")
    def __len__(self):
        res = self._length(self)
        return res[0]

    _extract = getmethod("[",
                         signature = StrVector(["SimpleList", 
                                                "ANY", "ANY"]),
                         where="package:IRanges")
    def __getitem__(self, i):
        res = self._extract(self, i)
        return res

class AnnotatedList(robjects.RS4):
    pass

class DataFrame(SimpleList):
    _extract = getmethod("[[",
                         signature = StrVector(["DataFrame", ]),
                         where="package:IRanges")

    def __getitem__(self, i):
        res = self._extract(self, i)
        return res


class XDataFrameList(AnnotatedList):
    pass

class RangesList(AnnotatedList):
    pass

class RangesDataList(AnnotatedList):
    pass

class IRangesList(RangesList):
    pass

class Views(IRanges):
    _extract = getmethod("[[",
                         signature = StrVector(["Views", ]),
                         where="package:IRanges")
    def __getitem__(self, i):
        res = self._extract(self, i)
        return res

    _subject = getmethod("subject", 
                         signature = StrVector(["Views", ]), 
                         where="package:IRanges")

    def get_subject(self):
        """ Get the 'subject' """
        res = self._subject(self)
        return res
    subject = property(get_subject, None)

class RleViews(Views):
    pass

class MaskCollection(robjects.RS4):
    _length = getmethod("length", 
                        signature = "MaskCollection",
                        where = "package:IRanges")
    _mask_constructor = iranges_env['Mask']
    _width = getmethod("width",
                       signature = "MaskCollection",
                       where = "package:IRanges")

    @classmethod
    def new(cls, mask_width):
        res = cls(cls._mask_constructor(mask_width))
        return cls
        
    def __len__(self):
        res = self._length(self)[0]
        return res

    def get_width(self):
        res = self._length(self)
        return res
        
class RangedData(robjects.RS4):
    _length = getmethod("length", 
                        signature = "RangedData",
                        where = "package:IRanges")
    _names = getmethod("names", 
                        signature = "RangedData",
                        where = "package:IRanges")
    _extract = getmethod("[[",
                         signature = StrVector(["RangedData", ]),
                         where="package:IRanges")
    def __len__(self):
        res = self._length(self)[0]
        return res

    def get_names(self):
        res = self._names(self)
        return res
    names = property(get_names, None)

    def __getitem__(self, i):
        res = self._extract(self, i)
        return res


acs = (('start', None, True, None),
       ('width', None, True, None),
       ('end', None, True, None),
       )
robjects.methods.set_accessors(RangedData,
                               "RangedData", "IRanges",
                               acs)
del(acs)


_iranges_dict = {
    'Sequence': Sequence,
    'Rle': Rle,
    'AtomicList': AtomicList,
    'LogicalList': LogicalList,
    'IntegerList': IntegerList,
    'NumericList': NumericList,
    'CharacterList': CharacterList,
    'ComplexList': ComplexList,
    'RawList': RawList,
    'RangedData': RangedData,
    'Ranges': Ranges,
    'IRanges': IRanges,
    'XRanges': XRanges,
    'IntervalTree': IntervalTree,
    'NormalIRanges': NormalIRanges,
    'Views': Views,
    'RangesList': RangesList,
    'RangesDataList': RangesDataList,
    'IRangesList': IRangesList,
    'RleViews': RleViews,
    'DataFrame': DataFrame,
    'XDataFrameList': XDataFrameList
    }

original_conversion = conversion.ri2py
def iranges_conversion(robj):

    pyobj = original_conversion(robj)
    if isinstance(pyobj, robjects.RS4):
        rclass = [x for x in pyobj.rclass][0]
        try:
            cls = _iranges_dict[rclass]
            pyobj = cls(pyobj)
        except KeyError, ke:
            pass

    return pyobj

conversion.ri2py = iranges_conversion

