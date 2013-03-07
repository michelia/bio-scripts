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
A module to model the GGtools library in Bioconductor

Copyright 2009 - Laurent Gautier

"""

__rname__ = 'GGBase'

import rpy2.robjects.methods
import rpy2.robjects as robjects
import rpy2.robjects.packages
from rpy2.robjects.packages import importr

conversion = robjects.conversion
getmethod = robjects.baseenv.get("getMethod")
import bioc.biobase

__rpackage__ = importr(__rname__)

StrVector = robjects.StrVector

class SmlSet(bioc.biobase.ESet):
    __rname__ = 'smlSet'
    __metaclass__ = rpy2.robjects.methods.RS4_Type

    __accessors__ = (('smList', 'GGBase', 'smlist', True, None),
                     ('smlEnv', 'GGBase', 'smlenv', True, None),
                     ('exprs', 'GGBase', 'exprs', True, None)
                     )
    _getalleles = getmethod("getAlleles", 
                            signature = StrVector(["smlSet", "rsid"]), 
                            where="package:GGBase")
    _snps = getmethod("snps", 
                      signature = StrVector(["smlSet", "chrnum"]), 
                      where="package:GGBase")
    def getalleles(self, rsid):
        res = self._getalleles(self, rsid)
        res = conversion.ri2py(res)
        return res
    def snps(self, chrnum):
        res = self._snps(self, chrnum)
        res = conversion.ri2py(res)
        return res


class GeneSymbol(bioc.biobase.CharacterOrMIAME):
    __rname__ = 'gensym'
    pass

class GwSnpScreenResult(bioc.biobase.AssayData):
    __rname__ = "gwSnpScreenResult"


class CwSnpScreenResult(GwSnpScreenResult):
    __rname__ = "cwSnpScreenResult"

class MultiGwSnpScreenResult(GwSnpScreenResult):
    __rname__ = "multiGwSnpScreenResult"

class FilteredMultiGwSnpScreenResult(MultiGwSnpScreenResult):
    __rname__ = "filteredMultiGwSnpScreenResult"


_ggbase_dict = {
    'smlSet': SmlSet,
    'genesym': GeneSymbol,
    'gwSnpScreenResult': GwSnpScreenResult,
    'cwSnpScreenResult': CwSnpScreenResult,
   ' multiGwSnpScreenResult': MultiGwSnpScreenResult,
    'filteredMultiGwSnpScreenResult': FilteredMultiGwSnpScreenResult
    }

original_conversion = robjects.conversion.ri2py
def ggbase_conversion(robj):

    pyobj = original_conversion(robj)
    if isinstance(pyobj, robjects.RS4):
        rclass = [x for x in pyobj.rclass][0]
        try:
            cls = _ggbase_dict[rclass]
            pyobj = cls(pyobj)
        except KeyError, ke:
            pass

    return pyobj

robjects.conversion.ri2py = ggbase_conversion
