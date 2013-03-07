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
A module to model the Biostrings library in Bioconductor

Copyright 2009-2010 - Laurent Gautier

"""


import rpy2.robjects.methods
import rpy2.robjects as robjects
import rpy2.robjects.packages
conversion = robjects.conversion
import iranges
from rpy2.robjects.packages import importr

base = importr('base')
methods = importr('methods')
getmethod = methods.getMethod
selectmethod = methods.selectMethod

__rname__ = 'Biostrings'

biostrings = importr(__rname__)
__rpackage__ = biostrings

StrVector = robjects.StrVector



def _setExtractDelegators(self):
    self.rx = robjects.vectors.ExtractDelegator(self)
    self.rx2 = robjects.vectors.DoubleExtractDelegator(self)        

   
class XString(rpy2.robjects.methods.RS4):
    """ Arbitrary string """

    __metaclass__ = rpy2.robjects.methods.RS4_Type

    __accessors__ = (('nchar', 'Biostrings', None, True, None),
                     ('reverse', 'Biostrings', None, False, None))

    _alphabet = getmethod("alphabet", 
                          signature = StrVector(["ANY", ]), 
                          where="package:Biostrings")


    def __init__(self, *args, **kwargs):
        super(XString, self).__init__(*args, **kwargs)
        _setExtractDelegators(self)
        
    def get_alphabet(self):
        res = self._alphabet(self)
        return res
    alphabet = property(get_alphabet, None, "Get the 'alphabet'")




class BString(XString):
    """ Biological string """
    _bstring_constructor = biostrings.BString

    @classmethod
    def new(cls, x):
        """ :param x: a (biological) string """
        res = cls(cls._bstring_constructor(conversion.py2ri(x)))
        _setExtractDelegators(res)
        return res


class DNAString(XString):
    ''' DNA string '''
    _dnastring_constructor = biostrings.DNAString
    _reverse_complement = getmethod("reverseComplement", 
                                    signature = StrVector(["DNAString", ]), 
                                    where="package:Biostrings")

    @classmethod
    def new(cls, x):
        """ :param x: a DNA string """
        res = cls(cls._dnastring_constructor(conversion.py2ri(x)))
        _setExtractDelegators(res)
        return res

    def reverse_complement(self):
        ''' Return the reverse complement '''
        res = self._reverse_complement(self)
        return res


    

class RNAString(XString):
    """ RNA string """
    _rnastring_constructor = biostrings.RNAString
    _reverse_complement = getmethod("reverseComplement", 
                                    signature = StrVector(["RNAString", ]), 
                                    where="package:Biostrings")

    @classmethod
    def new(cls, x):
        """ :param x: an RNA string """
        res = cls(cls._rnastring_constructor(conversion.py2ri(x)))
        _setExtractDelegators(res)
        return res

    def reverse_complement(self):
        ''' Return the reverse complement '''
        res = self._reverse_complement(self)
        return res


class AAString(XString):
    _aastring_constructor = biostrings.AAString

    @classmethod
    def new(cls, x):
        """ :param x: a string of amino-acids """
        res = cls(cls._aastring_constructor(conversion.py2ri(x)))
        _setExtractDelegators(res)
        return res


class XStringSet(rpy2.robjects.methods.RS4):
    """ An abstract class, parent to *StringSet classes """
    _subset = selectmethod("[", 
                           signature = "XStringSet")
    _length = base.length

    r = None
    
    def subset(self, x):
        res = self._subset(self, x)
        return res

    #FIXME: should implement an rpy2-defined sequence interface ?
    def __len__(self):
        res = self._length(self)[0]
        return res

class BStringSet(XStringSet):
    ''' Set of biological strings '''
    _bstringset_constructor = biostrings.BStringSet
    @classmethod
    def new(cls, x):
        res = cls(cls._bstringset_constructor(conversion.py2ri(x)))
        _setExtractDelegators(res)
        return res


class DNAStringSet(XStringSet):
    ''' Set of DNA strings '''
    _dnastringset_constructor = biostrings.DNAStringSet

    @classmethod
    def new(cls, x):
        res = cls(cls._dnastringset_constructor(conversion.py2ri(x)))
        _setExtractDelegators(res)
        return res


class RNAStringSet(XStringSet):
    ''' Set of RNA strings '''
    _rnastringset_constructor = biostrings.RNAStringSet
    @classmethod
    def new(cls, x):
        res = cls(cls._rnastringset_constructor(conversion.py2ri(x)))
        _setExtractDelegators(res)
        return res


class AAStringSet(XStringSet):
    ''' Set of amino-acid strings '''
    _aastringset_constructor = biostrings.AAStringSet
    @classmethod
    def new(cls, x):
        res = cls(cls._aastringset_constructor(conversion.py2ri(x)))
        _setExtractDelegators(res)
        return res

class AlignedXStringSet0(rpy2.robjects.methods.RS4):
    __metaclass__ = rpy2.robjects.methods.RS4_Type
    __accessors__ = (('unaligned', 'Biostrings', None, True, None),
                     ('start', 'Biostrings', None, True, None),
                     ('end', 'Biostrings', None, True, None),
                     ('width', 'Biostrings', None, True, None),
                     ('indel', 'Biostrings', None, False, None),
                     ('length', 'Biostrings', '__len__', False, None),
                     ('nchar', 'Biostrings', None, True, None),
                     #('alphabet', None, True, ''),
                     ('as.character', 'Biostrings', 'as_character',
                      False, None))

    _alphabet = getmethod("alphabet", 
                          signature = StrVector(["ANY", ]), 
                          where="package:Biostrings")
    def get_alphabet(self):
        res = self._alphabet(self)
        return res
    alphabet = property(get_alphabet, None, "Get the 'alphabet'")


class AlignedXStringSet(AlignedXStringSet0):
    pass


class XStringViews(iranges.Views):
    """ View on an arbitrary string """

    
    _nchar = getmethod("nchar", 
                       signature = "XStringViews",
                       where="package:Biostrings")
    _reverse_complement = getmethod("reverseComplement", 
                                    signature = StrVector(["XStringViews", ]), 
                                    where="package:Biostrings")
    _as_matrix = getmethod("as.matrix", 
                           signature = StrVector(["XStringViews", ]), 
                           where="package:Biostrings")

        

    def get_nchar(self):
        res = self._nchar(self)
        return res
    nchar = property(get_nchar, None,
                     "Number of characters in the view")
    width = property(get_nchar, None,
                     "Number of characters in the view (identical to nchar)")

    
    def reverse_complement(self):
        """ Return the reverse and complement for the view """
        res = self._reverse_complement(self)
        return res

    def as_matrix(self, **kwargs):
        res = self._as_matrix(self, **kwargs)
        return res

class MaskedXString(rpy2.robjects.methods.RS4):
    ''' "Masked" arbitrary string '''
    _unmasked = getmethod("unmasked", 
                          signature = "MaskedXString",
                          where="package:Biostrings")
    _masks = getmethod("masks", 
                       signature = "MaskedXString",
                       where="package:Biostrings")
    _masks_set_maskcollection = getmethod("masks<-", 
                                          signature = robjects.StrVector(("MaskedXString", 
                                                                          "MaskCollection")),
                                          where="package:Biostrings")
    

    def __length__(self):
        res = self._length(self)
        return res

    def get_unmasked(self):
        """ Return the strings without its 'mask' """
        res = self._unmasked(self)
        return res
    unmasked = property(get_unmasked, None)

    def get_masks(self):
        res = res._masks(self)
        return res

    def set_masks(self, value):
        if isinstance(value, MaskCollection):
            res = res._masks_set_maskcollection(self, value)
        else:
            raise ValueError('Value should be a MaskCollection')
        return res
    masks = property(get_masks, set_masks,
                     None,
                     "Property for both R's 'masks' and 'masks<-'")


class MaskedBString(MaskedXString):
    """ Masked biological string """
    pass

class MaskedDNAString(MaskedXString):
    """ Masked DNA string """
    pass

class MaskedRNAString(MaskedXString):
    """ Masked RNA string """
    pass

class MaskedAAString(MaskedXString):
    """ Masked string of amino-acids """
    pass

class XStringQuality(rpy2.robjects.methods.RS4):
    pass

class PhredQuality(XStringQuality, XStringSet):
    _constructor = biostrings.PhredQuality

    @classmethod
    def new(cls, x):
        res = cls(cls._constructor(conversion.py2ri(x)))
        return res

    
class SolexaQuality(XStringQuality, XStringSet):
    _constructor = biostrings.SolexaQuality

    @classmethod
    def new(cls, x):
        res = cls(cls._constructor(conversion.py2ri(x)))
        return res

class QualityScaledDNAStringSet(XStringSet):
    _constructor = biostrings.QualityScaledDNAStringSet
    @classmethod
    def new(cls, x, quality):
        res = cls(cls._constructor(conversion.py2ri(x)))
        return res


class QualityScaledXStringSet(XStringSet):
    pass

# acs = (
#        ('quality', None, True, None),
#       )
# rpy2.robjects.methods.set_accessors(QualityScaledXStringSet,
#                                     "QualityScaledXStringSet", "Biostrings",
#                                     acs)
# del(acs)

class QualityScaledDNAStringSet(DNAStringSet, QualityScaledXStringSet):
    _constructor = biostrings.QualityScaledDNAStringSet
    @classmethod
    def new(cls, x, quality):
        res = cls(cls._constructor(conversion.py2ri(x)))
        return res

class QualityScaledRNAStringSet(RNAStringSet, QualityScaledXStringSet):
    _constructor = biostrings.QualityScaledRNAStringSet
    @classmethod
    def new(cls, x, quality):
        res = cls(cls._constructor(conversion.py2ri(x)))
        return res

class QualityScaledBStringSet(BStringSet, QualityScaledXStringSet):
    _constructor = biostrings.QualityScaledBStringSet
    @classmethod
    def new(cls, x, quality):
        res = cls(cls._constructor(conversion.py2ri(x)))
        return res

class QualityScaledAAStringSet(AAStringSet, QualityScaledXStringSet):
    _constructor = biostrings.QualityScaledAAStringSet
    @classmethod
    def new(cls, x, quality):
        res = cls(cls._constructor(conversion.py2ri(x)))
        return res


    
    
class InDel(rpy2.robjects.methods.RS4):
    __metaclass__ = rpy2.robjects.methods.RS4_Type
    __accessors__ = (('insertion', 'Biostrings', None, True, None),
                     ('deletion', 'Biostrings', None, True, None))


class PairwiseAlignedXStringSet(rpy2.robjects.methods.RS4):
    __metaclass__ = rpy2.robjects.methods.RS4_Type
    __accessors__ = (('score', 'Biostrings', None, True, None),
                     ('nchar', 'Biostrings', None, True, None),
                     ('nindel', 'Biostrings', None, True, None),
                     ('pattern', 'Biostrings', None, True, None),
                     )
        

    _fromXString_XString = \
            getmethod("PairwiseAlignedXStringSet", 
                      signature = robjects.StrVector(("XString",
                                                      "XString")),
                      where="package:Biostrings")
    _fromCharacter_Character = \
            getmethod("PairwiseAlignedXStringSet", 
                      signature = robjects.StrVector(("character",
                                                      "character")),
                                       where="package:Biostrings")

    _fromXStringSet_missing = \
            getmethod("PairwiseAlignedXStringSet", 
                      signature = robjects.StrVector(("XStringSet",
                                                      "missing")),
                                   where="package:Biostrings")
    _fromCharacter_missing = \
            getmethod("PairwiseAlignedXStringSet", 
                      signature = robjects.StrVector(("character",
                                                      "missing")),
                      where="package:Biostrings")
    

    @staticmethod
    def fromXString_XString(pattern, target, **kwargs):
        res = PairwiseAlignedXStringSet._fromXString_XString(pattern,
                                                             target)
        return res

    @staticmethod
    def fromCharacter_Character(pattern, target, **kwargs):
        res = PairwiseAlignedXStringSet._fromCharacter_Character(pattern,
                                                                 target,
                                                                 **kwargs)
        return res

    @staticmethod
    def fromCharacter_missing(pattern, **kwargs):
        res = PairwiseAlignedXStringSet._fromCharacter_missing(pattern,
                                                               **kwargs)
        return res


class PairwiseAlignedFixedSubject(PairwiseAlignedXStringSet):
    pass

class PairwiseAlignedFixedSubjectSummary(rpy2.robjects.methods.RS4):
    pass

class PDict(rpy2.robjects.methods.RS4):
    """ Dictionnary of probes, that is dictionary of
    of rather short strings. """
    
    _length = getmethod("length", 
                        signature = "PDict",
                        where = "package:Biostrings")
    _pdict_constructor = biostrings.PDict
    _match = biostrings._env['matchPDict']
    _count = biostrings._env['countPDict']
    _which = biostrings._env['whichPDict']
    _width = getmethod("width",
                       signature = "PDict",
                       where = "package:Biostrings")

    @classmethod
    def create_instance(self, x):
        """ Create a preprocessed dictionnary of genomic patterns.

        :param x: a string vector, and DNAStringSet, 
                  or an XStringViews with s DNAString subject"""
        res = self._pdict_constructor(conversion.py2ri(x))
        return res
        
    def __length__(self):
        res = self._length(self)
        return res

    width = property(__length__, None,
                     "Number of elements in the dictionary (identical to length)")

    
    def match(self, subject, algorithm = "auto", 
              max_mismatch = 0, fixed = True, verbose = False):
        """ Match subject sequence(s) to the dictionary """
        res = self._match(self, subject, 
                          **{'algorithm': algorithm, 
                             'max.mismatch': max_mismatch,
                             'fixed': fixed, 
                             'verbose': verbose})
        return res

    def count(self, subject, algorithm = "auto", 
              max_mismatch = 0, fixed = True, verbose = False):
        """ Count the number of matching subject sequences """
        res = self._match(self, subject, 
                          **{'algorithm': algorithm, 
                             'max.mismatch': max_mismatch,
                             'fixed': fixed, 
                             'verbose': verbose})
        return res
    def which(self, subject, algorithm = "auto", 
              max_mismatch = 0, fixed = True, verbose = False):
        res = self._match(self, subject, 
                          **{'algorithm': algorithm, 
                             'max.mismatch': max_mismatch,
                             'fixed': fixed, 
                             'verbose': verbose})
        return res

class TB_PDict(PDict):
    """ 'Trusted-band' (TB) probe dictionary """
    _tb = getmethod("tb", 
                    signature = "TB_PDict",
                    where = "package:Biostrings")
    _tb_width = getmethod("tb.width", 
                          signature = "TB_PDict",
                          where = "package:Biostrings")
    def get_tb(self):
        res = self._tb(self)
        return res
    tb = property(get_tb, None,
                  "Return the trusted band for the dictionary")

    def get_tb_width(self):
        res = self._tb_width(self)
        return res
    tb_width = property(get_tb_width, None,
                        "Return the width of the trusted band")
   

class MTB_PDict(PDict):
    pass


_biostrings_dict = {
    'BString': BString, 
    'DNAString': DNAString,
    'RNAString': RNAString,
    #
    #'Views': Views,
    'XStringViews': XStringViews,
    #
    #'MaskCollection': MaskCollection,
    'MaskedXString': MaskedXString,
    'MaskedBString': MaskedBString,
    'MaskedDNAString': MaskedDNAString,
    'MaskedRNAString': MaskedRNAString,
    #
    'PDict': PDict,
    'TB_PDict': TB_PDict,
    'MTB_PDict': MTB_PDict,
    #
    'XStringSet': XStringSet,
    'BStringSet': BStringSet,
    'DNAStringSet': DNAStringSet,
    'RNAStringSet': RNAStringSet,
    'AAStringSet': AAStringSet,
    #
    'AlignedXStringSet0': AlignedXStringSet0,
    'AlignedXStringSet': AlignedXStringSet,
    #
    'XStringQuality': XStringQuality,
    'PhredQuality': PhredQuality,
    'SolexaQuality': SolexaQuality,
    #
    'QualityScaledXStringSet': QualityScaledXStringSet,
    'QualityScaledBStringSet': QualityScaledBStringSet,
    'QualityScaledDNAStringSet': QualityScaledDNAStringSet,
    'QualityScaledRNAStringSet': QualityScaledRNAStringSet,
    'QualityScaledAAStringSet': QualityScaledAAStringSet,
    #
    'InDel': InDel,
    #
    'PairwiseAlignedXStringSet': PairwiseAlignedXStringSet,
    'PairwiseAlignedFixedSubject': PairwiseAlignedFixedSubject,
    'PairwiseAlignedFixedSubjectSummary': PairwiseAlignedFixedSubjectSummary
    }

original_conversion = conversion.ri2py
def biostrings_conversion(robj):

    pyobj = original_conversion(robj)
    if isinstance(pyobj, robjects.RS4):
        rclass = [x for x in pyobj.rclass][0]
        try:
            cls = _biostrings_dict[rclass]
            pyobj = cls(pyobj)
        except KeyError, ke:
            pass

    return pyobj

conversion.ri2py = biostrings_conversion

