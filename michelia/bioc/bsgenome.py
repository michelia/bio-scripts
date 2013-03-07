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
A module to model the BSgenome library in Bioconductor

Copyright 2009 - Laurent Gautier

"""


import rpy2.robjects.methods
import rpy2.robjects as robjects
import rpy2.robjects.packages
import bioc.iranges, bioc.biostrings


getmethod = robjects.baseenv.get("getMethod")


__rname__ = 'BSgenome'
rpy2.robjects.packages.quiet_require(__rname__)

bsgenome_env = robjects.baseenv['as.environment']('package:BSgenome')
__rpackage__ = robjects.packages.SignatureTranslatedPackage(bsgenome_env, __rname__)

StrVector = robjects.StrVector


class GenomeDescription(rpy2.robjects.methods.RS4):
    __metaclass__ = rpy2.robjects.methods.RS4_Type
    
    __accessors__ = (('organism', 'BSgenome', None, True, None),
                     ('species', 'BSgenome', None, True, None),
                     ('provider', 'BSgenome', None, True, None),
                     ('providerVersion', 'BSgenome', 'provider_version',
                      True, None),
                     ('releaseDate', 'BSgenome', 'release_date', True, None),
                     ('releaseName', 'BSgenome', 'release_name', True, None),
                     )



class BSgenome(GenomeDescription):
    """ Arbitrary string """

    __accessors__ = (('seqnames', 'BSgenome', None, True, None),
                     ('mseqnames', 'BSgenome', None, True, None),
                     ('seqlengths', 'BSgenome', None, True, None),
                     )
    
    _extract = getmethod("[[", 
                         signature = StrVector(["BSgenome", "ANY", "ANY"]),
                         where="package:BSgenome")


    def __getitem__(self, x):
        if isinstance(x, int):
            x = x - 1
        elif not (isinstance(x, str) or isinstance(x, unicode)):
            raise ValueError('x should either be an integer or a string/unicode')
        res = self._extract(self, x)
        return res

class GenomeData(bioc.iranges.AnnotatedList):
    __metaclass__ = rpy2.robjects.methods.RS4_Type
    
    __accessors__ = (('organism', 'BSgenome', None, True, None),
                     ('provider', 'BSgenome', None, True, None),
                     ('providerVersion', 'BSgenome', 'provider_version',
                      True, None),
                     )

class GenomeDataList(bioc.iranges.AnnotatedList):
    pass

_bsgenome_dict = {
    'GenomeDescription': GenomeDescription,
    'BSgenome': BSgenome,
    'GenomeData': GenomeData,
    'GenomeDataList': GenomeDataList
    }

original_conversion = robjects.conversion.ri2py
def bsgenome_conversion(robj):

    pyobj = original_conversion(robj)
    if isinstance(pyobj, robjects.RS4):
        rclass = pyobj.rclass[0]
        try:
            cls = _bsgenome_dict[rclass]
            pyobj = cls(pyobj)
        except KeyError, ke:
            pass

    return pyobj

robjects.conversion.ri2py = bsgenome_conversion

