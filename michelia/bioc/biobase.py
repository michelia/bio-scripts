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
A module to model the Biobase library in Bioconductor

Copyright 2009-2010 - Laurent Gautier

"""

import rpy2.robjects.methods
import rpy2.robjects as robjects
import rpy2.robjects.packages
import rpy2.rlike.container as rlc

getmethod = robjects.baseenv.get("getMethod")

__rname__ = 'Biobase'
robjects.packages.quiet_require(__rname__)

biobase_env = robjects.baseenv['as.environment']('package:Biobase')
__rpackage__ = robjects.packages.SignatureTranslatedPackage(biobase_env, __rname__)

new = robjects.conversion.ri2py(robjects.methods.methods_env['new'])

StrVector = robjects.StrVector


def _setExtractDelegators(self):
    self.rx = robjects.vectors.ExtractDelegator(self)
    self.rx2 = robjects.vectors.DoubleExtractDelegator(self)        


class Versioned(rpy2.robjects.methods.RS4):
    __metaclass__ = rpy2.robjects.methods.RS4_Type
    
    __accessors__ = (('classVersion', 'Biobase', '__version__', True,
                      'maps Biobase::classVersion'),
                     )

class VersionedBiobase(Versioned):
    pass

class CharacterOrMIAME(rpy2.robjects.methods.RS4):
    """ Abstract class """
    pass

class AssayData(rpy2.robjects.methods.RS4):
    """ Abstract class. That class in a ClassUnionRepresentation
     in R, that a is way to create a parent class for existing classes.
     This is currently not modelled in Python. """
    __rname__ = 'AssayData'
    __metaclass__ = rpy2.robjects.methods.RS4_Type
    
    __accessors__ = (('featureNames', 'Biobase', 'featurenames', True,
                      'maps Biobase::featureNames'),
                     ('sampleNames', 'Biobase', 'samplenames', True,
                      'maps Biobase::samplenames'),
                     ('storageMode', 'Biobase', 'storagemode', True,
                      'maps Biobase::storageMode')
                     )



class AnnotatedDataFrame(Versioned):
    ''' An annotated data.frame as defined in the R package Biobase'''

    _combine = getmethod("combine", 
                         signature = StrVector(["AnnotatedDataFrame", 
                                                "AnnotatedDataFrame"]),
                         where="package:Biobase")
    _get_pData = getmethod("pData", 
                           signature = "AnnotatedDataFrame",
                           where="package:Biobase")
    _set_pData = getmethod("pData<-", 
                           signature = StrVector(["AnnotatedDataFrame",
                                                  "data.frame"]),
                           where="package:Biobase")
    _get_featureNames = getmethod("featureNames", 
                                  signature = "AnnotatedDataFrame",
                                  where="package:Biobase")
    _set_featureNames = getmethod("featureNames<-", 
                                  signature = StrVector(["AnnotatedDataFrame", ]),
                                  where="package:Biobase")


    @classmethod
    def new(cls,
            data = robjects.DataFrame(rlc.OrdDict()),
            varmetadata = robjects.DataFrame(rlc.OrdDict()),
            dimlabels = robjects.StrVector(("rowNames", "columnNames"))
            ):

        res = new("AnnotatedDataFrame", data = data,
                  varMetadata = varmetadata,
                  dimLabels = dimlabels)
        _setExtractDelegators(res)
        return res

    def combine(self, annotated_dataframe):
        """ Combine the instance with an other data.frame """
        res = self._combine(self, annotated_dataframe)
        return res

    def get_pdata(self):
        res = self._get_pData(self)
        return res

    def set_pdata(self, value):
        res = self._set_pData(self, value)
        self.__sexp__ = res.__sexp__

    pdata = property(get_pdata, set_pdata,
                     None,
                     "Property for both R's 'pData' and 'pData<-'")

    def get_featurenames(self):
        res = self._get_featureNames(self)
        return res

    def set_featurenames(self, value):
        res = self._set_featureNames(self, value)
        self.__sexp__ = res.__sexp__

    featurenames = property(get_featurenames, 
                            set_featurenames,
                            None,
                            "Property for both R's 'featureNames' and 'featureNames<-'")


class ESet(VersionedBiobase):
    ''' An eSet as defined in the R package Biobase.
        This class in defined as "virtual" in the R/S4 scheme,
        which can translates as "abstract" in more common OO terminologies. '''

    _assayData = getmethod("assayData", 
                           signature = StrVector(["eSet", ]), 
                           where="package:Biobase")
    _phenoData_get = getmethod("phenoData", 
                               signature = StrVector(["eSet", ]), 
                               where="package:Biobase")
    _phenoData_set = getmethod("phenoData<-", 
                               signature = StrVector(["eSet", 
                                                      "AnnotatedDataFrame"]), 
                               where="package:Biobase")
    _featureData_get = getmethod("featureData", 
                                 signature = StrVector(["eSet", ]), 
                                 where="package:Biobase")
    _featureData_set = getmethod("featureData<-", 
                                 signature = StrVector(["eSet", 
                                                        "AnnotatedDataFrame"]), 
                                 where="package:Biobase")
    _experimentData = getmethod("experimentData", 
                                signature = StrVector(["eSet", ]), 
                                where="package:Biobase")
    _annotation = getmethod("annotation", 
                            signature = StrVector(["eSet", ]), 
                            where="package:Biobase")

    def assaydata(self):
	res = self._assayData(self)
        return res

    def get_phenodata(self):
        res = self._phenoData_get(self)
        return res
    def set_phenodata(self, value):
        sexp = self._phenoData_set(self, value)
        self.__sexp__ = sexp
    property(get_phenodata, set_phenodata,
             None,
             "Property for both R's 'phenoData' and 'phenoData<-'")             
    
    def get_featuredata(self):
	res = self._featureData(self)
        return res
    def set_featuredata(self, value):
	sexp = self._featureData_set(self, value)
        self.__sexp__ = sexp
    property(get_featuredata, set_featuredata,
             None,
             "Property for both R's 'featureData' and 'featureData<-'")
        
    def experimentdata(self):
	res = self._experimentData(self)
        return res

    def annotation(self):
	res = self._annotation(self)
        return res

class ExpressionSet(ESet):

    _get_exprs = getmethod("exprs", 
                           signature = StrVector(["ExpressionSet", ]), 
                           where="package:Biobase")
    _set_exprs = getmethod("exprs<-", 
                           signature = StrVector(["ExpressionSet", "matrix"]), 
                           where="package:Biobase")
    _esApply = getmethod("esApply", 
                         signature = StrVector(["ExpressionSet", ]), 
                         where="package:Biobase")

    def __init__(self, *args, **kwargs):
        super(ExpressionSet, self).__init__(*args, **kwargs)
        self.rx = None
        self.rx2 = None
        
    @classmethod
    def new(cls,
            phenodata = new("AnnotatedDataFrame"),
            featuredata = new("AnnotatedDataFrame"),
            experimentdata = new("MIAME"),
            annotation = robjects.StrVector([]),
            exprs = new("matrix")):

        res = new("ExpressionSet", phenoData = phenodata,
                  featureData = featuredata,
                  experimentData = experimentdata,
                  annotation = annotation,
                  exprs = exprs)
        _setExtractDelegators(res)
        return res
            
    def get_exprs(self):
	res = self._get_exprs(self)
        return res 

    def set_exprs(self, value):
        res = self._set_exprs(self, value)
        self.__sexp__ = res.__sexp__

    exprs = property(get_exprs, set_exprs,
                     None,
                     "Property for both R's 'exprs' and 'exprs<-'")
    
    def esapply(self, margin, fun, *args):
        ''' "apply" a function "fun" along the "margin". '''
        #FIXME: kwargs not (yet) handled
        res = self._esapply(self, margin, fun, *args)
        return res

class SnpSet(ESet):
    pass

class NChannelSet(ESet):
    __metaclass__ = rpy2.robjects.methods.RS4_Type
    
    __accessors__ = (('channelNames', 'Biobase', 'channelnames', 
                      True, 'maps Biobase::channelNames'),
                     )

    __channel = getmethod("channel", 
                          signature = StrVector(("NChannelSet", 
                                                 "character")),
                          where = "Biobase")
    
    __selectchannels = getmethod("selectChannels", 
                                 signature = StrVector(("NChannelSet", 
                                                        "character")),
                                 where = "Biobase")
    

    def channel(self, name, **kwargs):
        res = self.__channel(name, **kwargs)
        res = conversion.ri2py(res)
        return res

    def selectchannels(self, names, **kwargs):
        res = self.__selectchannels(names, **kwargs)
        res = conversion.ri2py(res)
        return res

class MultiSet(ESet):
    pass

class MIAME(CharacterOrMIAME):
    __metaclass__ = rpy2.robjects.methods.RS4_Type
    
    __accessors__ = (('abstract', 'Biobase', 'abstract', True,
                      'maps Biobase::abstract'),
                     ('expinfo', 'Biobase', 'expinfo', True,
                      'maps Biobase::expinfo'),
                     ('hybridizations', 'Biobase', 'hybridizations', True,
                      'maps Biobase::hybridizations'),
                     ('normControls', 'Biobase', 'normcontrols', True,
                      'maps Biobase::normControls'),
                     ('notes', 'Biobase', 'notes', True,
                      'maps Biobase::notes'),
                     ('otherInfo', 'Biobase', 'otherinfo', True,
                      'maps Biobase::otherInfo'),
                     ('preproc', 'Biobase', 'preproc', True,
                      'maps Biobase::preproc'),
                     ('samples', 'Biobase', 'samples', True,
                      'maps Biobase::samples'),
                     )


_biobase_dict = {
    'AnnotatedDataFrame': AnnotatedDataFrame,
    'AssayData': AssayData,
    'chararacterORMIAME': CharacterOrMIAME,
    'eSet': ESet,
    'ExpressionSet': ExpressionSet,
    'MIAME': MIAME,
    'SnpSet': SnpSet,
    'Versioned': Versioned,
    'VersionedBiobase': VersionedBiobase
    }

original_conversion = robjects.conversion.ri2py
def biobase_conversion(robj):

    pyobj = original_conversion(robj)
    if isinstance(pyobj, robjects.RS4):
        rclass = [x for x in pyobj.rclass][0]
        try:
            cls = _biobase_dict[rclass]
            pyobj = cls(pyobj)
        except KeyError, ke:
            pass

    return pyobj

robjects.conversion.ri2py = biobase_conversion



