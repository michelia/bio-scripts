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

__rname__ = 'GGtools'

import rpy2.robjects.methods
import rpy2.robjects as robjects
import rpy2.robjects.packages
from rpy2.robjects.packages import importr

conversion = robjects.conversion
getmethod = robjects.baseenv.get("getMethod")

import bioc.biostrings, bioc.iranges, bioc.ggbase
__rpackage__ = importr(__rname__)


StrVector = robjects.StrVector

# disappeared without notice from the bioconductor package GGtools
# class HbTestResults(rpy2.robjects.methods.RS4):
#     __rname__ = 'hbTestResults'
#     __metaclass__ = rpy2.robjects.methods.RS4_Type
#     __accessors__ = (('pvals', 'GGtools', 'pvals', True, None),
#                      ('locs', 'GGtools', None, True, None),
#                      ('hscores', 'GGtools', None, False, None))

class MultffManager(rpy2.robjects.methods.RS4):
    __rname__ = 'multffManager'
    __metaclass__ = rpy2.robjects.methods.RS4_Type

class MaxChisq(rpy2.robjects.methods.RS4):
    __rname__ = 'maxchisq'
    __metaclass__ = rpy2.robjects.methods.RS4_Type

_ggtools_dict = {
#    'hbTestResults': HbTestResults,
    'multffManager': MultffManager,
    'maxchisq': MaxChisq
    }

original_conversion = robjects.conversion.ri2py
def ggtools_conversion(robj):

    pyobj = original_conversion(robj)
    if isinstance(pyobj, robjects.RS4):
        rclass = [x for x in pyobj.rclass][0]
        try:
            cls = _ggtools_dict[rclass]
            pyobj = cls(pyobj)
        except KeyError, ke:
            pass

    return pyobj

robjects.conversion.ri2py = ggtools_conversion
 
