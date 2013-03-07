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

import rpy2.robjects.methods
import rpy2.robjects as robjects
import rpy2.robjects.packages

getmethod = robjects.baseenv.get("getMethod")


__rname__ = 'AnnotationDbi'
rpy2.robjects.packages.quiet_require(__rname__)

package_env = robjects.baseenv['as.environment']('package:AnnotationDbi')
__rpackage__ = robjects.packages.SignatureTranslatedPackage(package_env, __rname__)

StrVector = robjects.StrVector

class Bimap(rpy2.robjects.methods.RS4):
    ''' A Bimap as defined in the R package AnnotationDbi '''

    _length = getmethod("length", 
                        signature = StrVector(["Bimap", ]),
                        where="package:AnnotationDbi")
    _keys = getmethod("keys", 
                      signature = StrVector(["Bimap", ]),
                      where="package:AnnotationDbi")
    _mappedkeys = getmethod("mappedkeys", 
                            signature = StrVector(["Bimap", ]),
                            where="package:AnnotationDbi")
    _count_mappedkeys = getmethod("count.mappedkeys", 
                                  signature = StrVector(["Bimap", ]),
                                  where="package:AnnotationDbi")
    _subset = getmethod("subset", 
                        signature = StrVector(["Bimap", ]),
                        where="package:AnnotationDbi")


    def __length__(self):
        res = self._length(self)[0]
        return res

    def keys(self):
        res = self._keys(self)
        return res

    def mappedkeys(self):
        res = self._mappedkeys(self)
        return res

    def count_mappedkeys(self):
        res = self._count_mappedkeys(self)
        return res

    def subset(self, lkeys = None, rkeys = None):
        res = self._subset(self, Lkeys = lkeys, Rkeys = rkeys)
        return res


class AnnDbBimap(Bimap):
    _lkeys_get = getmethod("Lkeys", 
                       signature = StrVector(["AnnDbBimap", ]),
                       where="package:AnnotationDbi")
    _lkeys_set = getmethod("Lkeys<-", 
                           signature = StrVector(["AnnDbBimap", ]),
                           where="package:AnnotationDbi")
    _mappedlkeys = getmethod("mappedLkeys", 
                             signature = StrVector(["AnnDbBimap", ]),
                             where="package:AnnotationDbi")
    _rkeys_get = getmethod("Rkeys", 
                           signature = StrVector(["AnnDbBimap", ]),
                           where="package:AnnotationDbi")
    _rkeys_set = getmethod("Rkeys<-", 
                           signature = StrVector(["AnnDbBimap", ]),
                           where="package:AnnotationDbi")
    _mappedrkeys = getmethod("mappedRkeys", 
                             signature = StrVector(["AnnDbBimap", ]),
                             where="package:AnnotationDbi")
    _subset = getmethod("subset", 
                        signature = StrVector(["AnnDbBimap", ]),
                        where="package:AnnotationDbi")

    def get_lkeys(self):
        res = self._lkeys_get(self)
        return res
    def set_lkeys(self, value):
        sexp = self._lkeys_set(self, value)
        self.__sexp__ = sexp
    lkeys = property(get_lkeys, set_lkeys, None)

    def mappedlkeys(self):
        res = self._mappedlkeys(self)
        return res

    def get_rkeys(self):
        res = self._rkeys_get(self)
        return res
    def set_rkeys(self, value):
        sexp = self._rkeys_set(self, value)
        self.__sexp__ = sexp
    rkeys = property(get_rkeys, set_rkeys, None)


    def mappedrkeys(self):
        res = self._mappedrkeys(self)
        return res

    def subset(self, lkeys = None, rkeys = None, objname = None):
        res = self._subset(self, Lkeys = lkeys, Rkeys = rkeys, 
                           objName = objname)
        return res


class AnnDbMap(AnnDbBimap):
    pass

class GOTerms(rpy2.robjects.methods.RS4):
    __metaclass__ = rpy2.robjects.methods.RS4_Type
    
    __accessors__ = (('GOID', 'AnnotationDbi', 'GOID', True,
                      'maps AnnotationDbi::GOID'),
                     ('Term', 'AnnotationDbi', 'term', True,
                      'maps AnnotationDbi::Term'),
                     ('Ontology', 'AnnotationDbi', 'ontology', True,
                      'maps AnnotationDbi::Ontology'),
                     ('Definition', 'AnnotationDbi', 'definition', True,
                      'maps AnnotationDbi::Definition'),
                     ('Synonym', 'AnnotationDbi', 'synonym', True,
                      'maps AnnotationDbi::Synonym'),
                     ('Secondary', 'AnnotationDbi', 'secondary', True,
                      'maps AnnotationDbi::Secondary'),
                     )


annotationdbi_dict = {
    'Bimap': Bimap,
    'AnnDbBimap': AnnDbBimap,
    'AnnDbMap': AnnDbMap
    }

original_conversion = robjects.conversion.ri2py
def annotationdbi_conversion(robj):
    pyobj = original_conversion(robj)
    if isinstance(pyobj, robjects.RS4):
        rclass = [x for x in pyobj.rclass][0]
        try:
            cls = annotationdbi_dict[rclass]
            pyobj = cls(pyobj)
        except KeyError, ke:
            pass

    return pyobj

robjects.conversion.ri2py = annotationdbi_conversion

