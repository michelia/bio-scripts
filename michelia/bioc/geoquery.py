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
The module reflects some of the content of the R/Bioconductor
package *GEOquery*. It defines Python-level classes
for the R/S4 classes, and gives otherwise access to
R-level commands the usual :mod:`rpy2:robjects` way.

In complement to this module, it is recommended to also use
'rpy2.robjects.packages.importr' to expose the full content of *GEOquery*

>>> import bioc.geoquery

Now all methods exposed by the R package are accessible from 'bio.geoquery.__rpackage__',
with Python classes defined in bioc.geoquery mapping the R classes.

"""

import rpy2.robjects as robjects
import rpy2.robjects.packages

conversion = robjects.conversion

getmethod = robjects.baseenv.get("getMethod")

__rname__ = 'GEOquery'
robjects.packages.quiet_require(__rname__)


geoquery_env = robjects.baseenv['as.environment']('package:GEOquery')
__rpackage__ = rpy2.robjects.packages.SignatureTranslatedPackage(geoquery_env, __rname__)

StrVector = robjects.StrVector

class GEOData(robjects.methods.RS4):
    """ R 'virtual' class (that is an abstract class) """
    __metaclass__ = robjects.methods.RS4_Type

    __accessors__ = (
        ('Table', 'GEOquery', 'table', True, "table"),
        ('dataTable', 'GEOquery', 'datatable', True, "data table"),
        ('Meta', 'GEOquery', 'meta', True, "meta-data (sample information)"),
        ('Columns', 'GEOquery', 'columns', True, "column descriptions"),
        ('Accession', 'GEOquery', 'accession', True, "GEO accession number")
        )

class GEODataTable(robjects.methods.RS4):
    __metaclass__ = robjects.methods.RS4_Type

    __accessors__ = (
        ('Table', 'GEOquery', 'table', True, 'table'),
        ('Columns', 'GEOquery', 'columns', True, 'column descriptions')
        )

class GPL(GEOData):
    """ GEO 'platforms' entity """
    pass

class GSM(GEOData):
    """ GEO 'samples' entity """
    pass

class GDS(GEOData):
    """ GEO 'datasets' entity """
    pass

_geoquery_dict = {
    'GEOData': GEOData,
    'GEODataTable': GEODataTable,
    'GSM': GSM,
    'GPL': GPL,
    'GDS': GDS
    }

original_conversion = conversion.ri2py
def geoquery_conversion(robj):
    pyobj = original_conversion(robj)
    if isinstance(pyobj, robjects.RS4):
        rclass = [x for x in pyobj.rclass][0]
        try:
            cls = _geoquery_dict[rclass]
            pyobj = cls(pyobj)
        except KeyError, ke:
            pass

    return pyobj

conversion.ri2py = geoquery_conversion
