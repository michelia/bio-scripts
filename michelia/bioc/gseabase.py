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
A module to model the GSEABase library in Bioconductor

Copyright 2009 - Laurent Gautier

"""

__rname__ = 'GSEABase'

import rpy2.robjects.methods
import rpy2.robjects as robjects
import rpy2.robjects.packages

rinterface = robjects.rinterface
conversion = robjects.conversion
getmethod = robjects.baseenv.get("getMethod")

import bioc.biobase, bioc.annotationdbi
rpy2.robjects.packages.quiet_require(__rname__)

package_env = robjects.baseenv['as.environment']('package:%s' %__rname__)
__rpackage__ = robjects.packages.SignatureTranslatedPackage(package_env, 
                                                            __rname__)

StrVector = robjects.StrVector

class GeneSet(rpy2.robjects.methods.RS4):
    __rname__ = 'GeneSet'
    __metaclass__ = rpy2.robjects.methods.RS4_Type
    
    __accessors__ = (('collectionType', 'GSEABase', 'collectiontype', True, None),
                     ('contributor', 'GSEABase', 'contributor', True, None),
                     ('creationDate', 'GSEABase', 'creationdate', True, None),
                     ('description', 'GSEABase', 'description', True, None),
                     ('geneIds', 'GSEABase', 'geneids', True, None),
                     ('longDescription', 'GSEABase', 'longdescription', True, None),
                     ('organism', 'GSEABase', 'organism', True, None),
                     ('pubMedIds', 'GSEABase', 'pubmedids', True, None),
                     ('setName', 'GSEABase', 'name', True, 'maps GSEABase::setName'),
                     ('geneIdType', 'GSEABase', 'geneidtype', True, 'maps GSEABase::geneIdType'),
                     ('urls', 'GSEABase', 'urls', True, None),
                     )

    __extract_num = getmethod("[",
                              signature = StrVector(["GeneSet", 
                                                     "numeric", "ANY", "ANY"]))

    __or = getmethod("|",
                     signature = StrVector(["GeneSet", "GeneSet"]))
    __and = getmethod("&",
                      signature = StrVector(["GeneSet", "GeneSet"]))
    __setdiff = getmethod("setdiff",
                          signature = StrVector(["GeneSet", "GeneSet"]))
                              #where="package:GSEABase") #likely bug in R's "methods" 
    #FIXME: duck-punching for the ExpressionSet class as the bioconductor
    # package defines a new "[" method using a GeneSet as an index
    def rx(self, i):
        res = self.__extract_num(self, i)
        res = conversion.ri2py(res)
        return res

    def union(self, gs):
        res = self.__or(self, gs)
        return conversion.rpi2py(res)

    def intersection(self, gs):
        res = self.__and(self, gs)
        return conversion.rpi2py(res)

    def difference(self, gs):
        res = self.__setdiff(self, gs)
        return conversion.rpi2py(res)


class GeneColorSet(GeneSet):
    __rname__ = 'GeneColorSet'
    __metaclass__ = rpy2.robjects.methods.RS4_Type

    __accessors__ = (('coloring', 'GSEABase', 'coloring', True, None),
                     ('geneColor', 'GSEABase', 'genecolor', True, None),
                     ('phenotypeColor', 'GSEABase', 'phenotypecolor', False, None),
                     ('phenotype', 'GSEABase', 'phenotype', False, None),
                     )

    __extract_num = getmethod("[",
                              signature = StrVector(['GeneColorSet', 
                                                     "numeric","ANY"]))
                              #where="package:GSEABase")
    __extract_str = getmethod("[",
                              signature = StrVector(['GeneColorSet', 
                                                     "character", "ANY"]))
                              #where="package:GSEABase")
    
    def rx(self, i):
        if inherits(i, str) or \
                (inherits(i, rinterface.SexpVector) and \
                     i.typeof == rinterface.STRSXP):
            res = self.__extract_str(self, i)
        else:
            res = self.__extract_num(self, i)
        res = conversion.ri2py(res)
        return res


class GeneSetCollection(rpy2.robjects.methods.RS4, rpy2.robjects.Vector):
    __rname__ = 'GeneSetCollection'
    __metaclass__ = rpy2.robjects.methods.RS4_Type

    __accessors__ = (('geneIds', 'GSEABase', 'geneids', True, None),
                     ('names', 'GSEABase', 'names', True, None)
                     )

_gseabase_dict = {
    'GeneSet': GeneSet,
    'GeneColorSet': GeneColorSet,
    'GeneSetCollection': GeneSetCollection
    }

original_conversion = robjects.conversion.ri2py
def gseabase_conversion(robj):

    pyobj = original_conversion(robj)
    if isinstance(pyobj, robjects.RS4):
        rclass = [x for x in pyobj.rclass][0]
        try:
            cls = _gseabase_dict[rclass]
            pyobj = cls(pyobj)
        except KeyError, ke:
            pass

    return pyobj

robjects.conversion.ri2py = gseabase_conversion
