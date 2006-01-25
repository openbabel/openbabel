import sys
import dl
sys.setdlopenflags(sys.getdlopenflags() | dl.RTLD_GLOBAL)

# This file was created automatically by SWIG 1.3.27.
# Don't modify this file, modify the SWIG interface instead.

import _openbabel

def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name) or (name == "thisown"):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


def _swig_setattr_nondynamic_method(set):
    def set_attr(self,name,value):
        if hasattr(self,name) or (name in ("this", "thisown")):
            set(self,name,value)
        else:
            raise AttributeError("You cannot add attributes to %s" % self)
    return set_attr


class OBGlobalDataBase(object):
    def __repr__(self):
        return "<%s.%s; proxy of C++ OpenBabel::OBGlobalDataBase instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        newobj = _openbabel.new_OBGlobalDataBase(*args)
        self.this = newobj.this
        self.thisown = 1
        del newobj.thisown
    def __del__(self, destroy=_openbabel.delete_OBGlobalDataBase):
        try:
            if self.thisown: destroy(self)
        except: pass

    def Init(*args): return _openbabel.OBGlobalDataBase_Init(*args)
    def GetSize(*args): return _openbabel.OBGlobalDataBase_GetSize(*args)
    def SetReadDirectory(*args): return _openbabel.OBGlobalDataBase_SetReadDirectory(*args)
    def SetEnvironmentVariable(*args): return _openbabel.OBGlobalDataBase_SetEnvironmentVariable(*args)
    def ParseLine(*args): return _openbabel.OBGlobalDataBase_ParseLine(*args)

class OBGlobalDataBasePtr(OBGlobalDataBase):
    def __init__(self, this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = OBGlobalDataBase
_openbabel.OBGlobalDataBase_swigregister(OBGlobalDataBasePtr)

class OBElement(object):
    def __repr__(self):
        return "<%s.%s; proxy of C++ OpenBabel::OBElement instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        newobj = _openbabel.new_OBElement(*args)
        self.this = newobj.this
        self.thisown = 1
        del newobj.thisown
    def GetAtomicNum(*args): return _openbabel.OBElement_GetAtomicNum(*args)
    def GetSymbol(*args): return _openbabel.OBElement_GetSymbol(*args)
    def GetCovalentRad(*args): return _openbabel.OBElement_GetCovalentRad(*args)
    def GetVdwRad(*args): return _openbabel.OBElement_GetVdwRad(*args)
    def GetMass(*args): return _openbabel.OBElement_GetMass(*args)
    def GetMaxBonds(*args): return _openbabel.OBElement_GetMaxBonds(*args)
    def GetElectroNeg(*args): return _openbabel.OBElement_GetElectroNeg(*args)
    def GetIonization(*args): return _openbabel.OBElement_GetIonization(*args)
    def GetElectronAffinity(*args): return _openbabel.OBElement_GetElectronAffinity(*args)
    def GetName(*args): return _openbabel.OBElement_GetName(*args)
    def GetRed(*args): return _openbabel.OBElement_GetRed(*args)
    def GetGreen(*args): return _openbabel.OBElement_GetGreen(*args)
    def GetBlue(*args): return _openbabel.OBElement_GetBlue(*args)
    def __del__(self, destroy=_openbabel.delete_OBElement):
        try:
            if self.thisown: destroy(self)
        except: pass


class OBElementPtr(OBElement):
    def __init__(self, this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = OBElement
_openbabel.OBElement_swigregister(OBElementPtr)

class OBElementTable(OBGlobalDataBase):
    def __repr__(self):
        return "<%s.%s; proxy of C++ OpenBabel::OBElementTable instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        newobj = _openbabel.new_OBElementTable(*args)
        self.this = newobj.this
        self.thisown = 1
        del newobj.thisown
    def __del__(self, destroy=_openbabel.delete_OBElementTable):
        try:
            if self.thisown: destroy(self)
        except: pass

    def ParseLine(*args): return _openbabel.OBElementTable_ParseLine(*args)
    def GetNumberOfElements(*args): return _openbabel.OBElementTable_GetNumberOfElements(*args)
    def GetSize(*args): return _openbabel.OBElementTable_GetSize(*args)
    def GetAtomicNum(*args): return _openbabel.OBElementTable_GetAtomicNum(*args)
    def GetSymbol(*args): return _openbabel.OBElementTable_GetSymbol(*args)
    def GetVdwRad(*args): return _openbabel.OBElementTable_GetVdwRad(*args)
    def GetCovalentRad(*args): return _openbabel.OBElementTable_GetCovalentRad(*args)
    def GetMass(*args): return _openbabel.OBElementTable_GetMass(*args)
    def CorrectedBondRad(*args): return _openbabel.OBElementTable_CorrectedBondRad(*args)
    def CorrectedVdwRad(*args): return _openbabel.OBElementTable_CorrectedVdwRad(*args)
    def GetMaxBonds(*args): return _openbabel.OBElementTable_GetMaxBonds(*args)
    def GetElectroNeg(*args): return _openbabel.OBElementTable_GetElectroNeg(*args)
    def GetIonization(*args): return _openbabel.OBElementTable_GetIonization(*args)
    def GetElectronAffinity(*args): return _openbabel.OBElementTable_GetElectronAffinity(*args)
    def GetRGB(*args): return _openbabel.OBElementTable_GetRGB(*args)
    def GetName(*args): return _openbabel.OBElementTable_GetName(*args)

class OBElementTablePtr(OBElementTable):
    def __init__(self, this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = OBElementTable
_openbabel.OBElementTable_swigregister(OBElementTablePtr)

class OBIsotopeTable(OBGlobalDataBase):
    def __repr__(self):
        return "<%s.%s; proxy of C++ OpenBabel::OBIsotopeTable instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        newobj = _openbabel.new_OBIsotopeTable(*args)
        self.this = newobj.this
        self.thisown = 1
        del newobj.thisown
    def __del__(self, destroy=_openbabel.delete_OBIsotopeTable):
        try:
            if self.thisown: destroy(self)
        except: pass

    def GetSize(*args): return _openbabel.OBIsotopeTable_GetSize(*args)
    def ParseLine(*args): return _openbabel.OBIsotopeTable_ParseLine(*args)
    def GetExactMass(*args): return _openbabel.OBIsotopeTable_GetExactMass(*args)

class OBIsotopeTablePtr(OBIsotopeTable):
    def __init__(self, this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = OBIsotopeTable
_openbabel.OBIsotopeTable_swigregister(OBIsotopeTablePtr)

class OBTypeTable(OBGlobalDataBase):
    def __repr__(self):
        return "<%s.%s; proxy of C++ OpenBabel::OBTypeTable instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        newobj = _openbabel.new_OBTypeTable(*args)
        self.this = newobj.this
        self.thisown = 1
        del newobj.thisown
    def __del__(self, destroy=_openbabel.delete_OBTypeTable):
        try:
            if self.thisown: destroy(self)
        except: pass

    def ParseLine(*args): return _openbabel.OBTypeTable_ParseLine(*args)
    def GetSize(*args): return _openbabel.OBTypeTable_GetSize(*args)
    def SetFromType(*args): return _openbabel.OBTypeTable_SetFromType(*args)
    def SetToType(*args): return _openbabel.OBTypeTable_SetToType(*args)
    def Translate(*args): return _openbabel.OBTypeTable_Translate(*args)
    def GetFromType(*args): return _openbabel.OBTypeTable_GetFromType(*args)
    def GetToType(*args): return _openbabel.OBTypeTable_GetToType(*args)

class OBTypeTablePtr(OBTypeTable):
    def __init__(self, this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = OBTypeTable
_openbabel.OBTypeTable_swigregister(OBTypeTablePtr)

class OBResidueData(OBGlobalDataBase):
    def __repr__(self):
        return "<%s.%s; proxy of C++ OpenBabel::OBResidueData instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        newobj = _openbabel.new_OBResidueData(*args)
        self.this = newobj.this
        self.thisown = 1
        del newobj.thisown
    def ParseLine(*args): return _openbabel.OBResidueData_ParseLine(*args)
    def GetSize(*args): return _openbabel.OBResidueData_GetSize(*args)
    def SetResName(*args): return _openbabel.OBResidueData_SetResName(*args)
    def LookupBO(*args): return _openbabel.OBResidueData_LookupBO(*args)
    def LookupType(*args): return _openbabel.OBResidueData_LookupType(*args)
    def AssignBonds(*args): return _openbabel.OBResidueData_AssignBonds(*args)
    def __del__(self, destroy=_openbabel.delete_OBResidueData):
        try:
            if self.thisown: destroy(self)
        except: pass


class OBResidueDataPtr(OBResidueData):
    def __init__(self, this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = OBResidueData
_openbabel.OBResidueData_swigregister(OBResidueDataPtr)

FILE_SEP_CHAR = _openbabel.FILE_SEP_CHAR
class OBFormat(object):
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ OpenBabel::OBFormat instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def ReadMolecule(*args): return _openbabel.OBFormat_ReadMolecule(*args)
    def ReadChemObject(*args): return _openbabel.OBFormat_ReadChemObject(*args)
    def WriteMolecule(*args): return _openbabel.OBFormat_WriteMolecule(*args)
    def WriteChemObject(*args): return _openbabel.OBFormat_WriteChemObject(*args)
    def Description(*args): return _openbabel.OBFormat_Description(*args)
    def TargetClassDescription(*args): return _openbabel.OBFormat_TargetClassDescription(*args)
    def GetType(*args): return _openbabel.OBFormat_GetType(*args)
    def SpecificationURL(*args): return _openbabel.OBFormat_SpecificationURL(*args)
    def GetMIMEType(*args): return _openbabel.OBFormat_GetMIMEType(*args)
    def Flags(*args): return _openbabel.OBFormat_Flags(*args)
    def SkipObjects(*args): return _openbabel.OBFormat_SkipObjects(*args)
    def MakeNewInstance(*args): return _openbabel.OBFormat_MakeNewInstance(*args)
    def __del__(self, destroy=_openbabel.delete_OBFormat):
        try:
            if self.thisown: destroy(self)
        except: pass


class OBFormatPtr(OBFormat):
    def __init__(self, this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = OBFormat
_openbabel.OBFormat_swigregister(OBFormatPtr)

class CharPtrLess(object):
    def __repr__(self):
        return "<%s.%s; proxy of C++ OpenBabel::CharPtrLess instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __call__(*args): return _openbabel.CharPtrLess___call__(*args)
    def __init__(self, *args):
        newobj = _openbabel.new_CharPtrLess(*args)
        self.this = newobj.this
        self.thisown = 1
        del newobj.thisown
    def __del__(self, destroy=_openbabel.delete_CharPtrLess):
        try:
            if self.thisown: destroy(self)
        except: pass


class CharPtrLessPtr(CharPtrLess):
    def __init__(self, this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = CharPtrLess
_openbabel.CharPtrLess_swigregister(CharPtrLessPtr)

class OBConversion(object):
    def __repr__(self):
        return "<%s.%s; proxy of C++ OpenBabel::OBConversion instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        newobj = _openbabel.new_OBConversion(*args)
        self.this = newobj.this
        self.thisown = 1
        del newobj.thisown
    def __del__(self, destroy=_openbabel.delete_OBConversion):
        try:
            if self.thisown: destroy(self)
        except: pass

    RegisterFormat = staticmethod(_openbabel.OBConversion_RegisterFormat)
    FindFormat = staticmethod(_openbabel.OBConversion_FindFormat)
    FormatFromExt = staticmethod(_openbabel.OBConversion_FormatFromExt)
    FormatFromMIME = staticmethod(_openbabel.OBConversion_FormatFromMIME)
    GetNextFormat = staticmethod(_openbabel.OBConversion_GetNextFormat)
    Description = staticmethod(_openbabel.OBConversion_Description)
    def GetInStream(*args): return _openbabel.OBConversion_GetInStream(*args)
    def GetOutStream(*args): return _openbabel.OBConversion_GetOutStream(*args)
    def SetInStream(*args): return _openbabel.OBConversion_SetInStream(*args)
    def SetOutStream(*args): return _openbabel.OBConversion_SetOutStream(*args)
    def SetInAndOutFormats(*args): return _openbabel.OBConversion_SetInAndOutFormats(*args)
    def SetInFormat(*args): return _openbabel.OBConversion_SetInFormat(*args)
    def SetOutFormat(*args): return _openbabel.OBConversion_SetOutFormat(*args)
    def GetInFormat(*args): return _openbabel.OBConversion_GetInFormat(*args)
    def GetOutFormat(*args): return _openbabel.OBConversion_GetOutFormat(*args)
    def GetInFilename(*args): return _openbabel.OBConversion_GetInFilename(*args)
    def GetInPos(*args): return _openbabel.OBConversion_GetInPos(*args)
    def GetInLen(*args): return _openbabel.OBConversion_GetInLen(*args)
    def GetTitle(*args): return _openbabel.OBConversion_GetTitle(*args)
    def GetAuxConv(*args): return _openbabel.OBConversion_GetAuxConv(*args)
    def SetAuxConv(*args): return _openbabel.OBConversion_SetAuxConv(*args)
    INOPTIONS = _openbabel.OBConversion_INOPTIONS
    OUTOPTIONS = _openbabel.OBConversion_OUTOPTIONS
    GENOPTIONS = _openbabel.OBConversion_GENOPTIONS
    def IsOption(*args): return _openbabel.OBConversion_IsOption(*args)
    def GetOptions(*args): return _openbabel.OBConversion_GetOptions(*args)
    def AddOption(*args): return _openbabel.OBConversion_AddOption(*args)
    def RemoveOption(*args): return _openbabel.OBConversion_RemoveOption(*args)
    def SetOptions(*args): return _openbabel.OBConversion_SetOptions(*args)
    RegisterOptionParam = staticmethod(_openbabel.OBConversion_RegisterOptionParam)
    GetOptionParams = staticmethod(_openbabel.OBConversion_GetOptionParams)
    def Convert(*args): return _openbabel.OBConversion_Convert(*args)
    def FullConvert(*args): return _openbabel.OBConversion_FullConvert(*args)
    def AddChemObject(*args): return _openbabel.OBConversion_AddChemObject(*args)
    def GetChemObject(*args): return _openbabel.OBConversion_GetChemObject(*args)
    def IsLast(*args): return _openbabel.OBConversion_IsLast(*args)
    def IsFirstInput(*args): return _openbabel.OBConversion_IsFirstInput(*args)
    def GetOutputIndex(*args): return _openbabel.OBConversion_GetOutputIndex(*args)
    def SetOutputIndex(*args): return _openbabel.OBConversion_SetOutputIndex(*args)
    def SetMoreFilesToCome(*args): return _openbabel.OBConversion_SetMoreFilesToCome(*args)
    def SetOneObjectOnly(*args): return _openbabel.OBConversion_SetOneObjectOnly(*args)
    GetDefaultFormat = staticmethod(_openbabel.OBConversion_GetDefaultFormat)
    def Write(*args): return _openbabel.OBConversion_Write(*args)
    def WriteString(*args): return _openbabel.OBConversion_WriteString(*args)
    def WriteFile(*args): return _openbabel.OBConversion_WriteFile(*args)
    def Read(*args): return _openbabel.OBConversion_Read(*args)
    def ReadString(*args): return _openbabel.OBConversion_ReadString(*args)
    def ReadFile(*args): return _openbabel.OBConversion_ReadFile(*args)
    BatchFileName = staticmethod(_openbabel.OBConversion_BatchFileName)
    IncrementedFileName = staticmethod(_openbabel.OBConversion_IncrementedFileName)

class OBConversionPtr(OBConversion):
    def __init__(self, this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = OBConversion
_openbabel.OBConversion_swigregister(OBConversionPtr)

OBConversion_RegisterFormat = _openbabel.OBConversion_RegisterFormat

OBConversion_FindFormat = _openbabel.OBConversion_FindFormat

OBConversion_FormatFromExt = _openbabel.OBConversion_FormatFromExt

OBConversion_FormatFromMIME = _openbabel.OBConversion_FormatFromMIME

OBConversion_GetNextFormat = _openbabel.OBConversion_GetNextFormat

OBConversion_Description = _openbabel.OBConversion_Description

OBConversion_RegisterOptionParam = _openbabel.OBConversion_RegisterOptionParam

OBConversion_GetOptionParams = _openbabel.OBConversion_GetOptionParams

OBConversion_GetDefaultFormat = _openbabel.OBConversion_GetDefaultFormat

OBConversion_BatchFileName = _openbabel.OBConversion_BatchFileName

OBConversion_IncrementedFileName = _openbabel.OBConversion_IncrementedFileName

NOTREADABLE = _openbabel.NOTREADABLE
READONEONLY = _openbabel.READONEONLY
READBINARY = _openbabel.READBINARY
ZEROATOMSOK = _openbabel.ZEROATOMSOK
NOTWRITABLE = _openbabel.NOTWRITABLE
WRITEONEONLY = _openbabel.WRITEONEONLY
WRITEBINARY = _openbabel.WRITEBINARY
DEFAULTFORMAT = _openbabel.DEFAULTFORMAT
class OBResidue(object):
    def __repr__(self):
        return "<%s.%s; proxy of C++ OpenBabel::OBResidue instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        newobj = _openbabel.new_OBResidue(*args)
        self.this = newobj.this
        self.thisown = 1
        del newobj.thisown
    def __del__(self, destroy=_openbabel.delete_OBResidue):
        try:
            if self.thisown: destroy(self)
        except: pass

    def AddAtom(*args): return _openbabel.OBResidue_AddAtom(*args)
    def InsertAtom(*args): return _openbabel.OBResidue_InsertAtom(*args)
    def RemoveAtom(*args): return _openbabel.OBResidue_RemoveAtom(*args)
    def Clear(*args): return _openbabel.OBResidue_Clear(*args)
    def SetName(*args): return _openbabel.OBResidue_SetName(*args)
    def SetNum(*args): return _openbabel.OBResidue_SetNum(*args)
    def SetChain(*args): return _openbabel.OBResidue_SetChain(*args)
    def SetChainNum(*args): return _openbabel.OBResidue_SetChainNum(*args)
    def SetIdx(*args): return _openbabel.OBResidue_SetIdx(*args)
    def SetAtomID(*args): return _openbabel.OBResidue_SetAtomID(*args)
    def SetHetAtom(*args): return _openbabel.OBResidue_SetHetAtom(*args)
    def SetSerialNum(*args): return _openbabel.OBResidue_SetSerialNum(*args)
    def GetName(*args): return _openbabel.OBResidue_GetName(*args)
    def GetNum(*args): return _openbabel.OBResidue_GetNum(*args)
    def GetNumAtoms(*args): return _openbabel.OBResidue_GetNumAtoms(*args)
    def GetChain(*args): return _openbabel.OBResidue_GetChain(*args)
    def GetChainNum(*args): return _openbabel.OBResidue_GetChainNum(*args)
    def GetIdx(*args): return _openbabel.OBResidue_GetIdx(*args)
    def GetResKey(*args): return _openbabel.OBResidue_GetResKey(*args)
    def GetAtoms(*args): return _openbabel.OBResidue_GetAtoms(*args)
    def GetBonds(*args): return _openbabel.OBResidue_GetBonds(*args)
    def GetAtomID(*args): return _openbabel.OBResidue_GetAtomID(*args)
    def GetSerialNum(*args): return _openbabel.OBResidue_GetSerialNum(*args)
    def GetAminoAcidProperty(*args): return _openbabel.OBResidue_GetAminoAcidProperty(*args)
    def GetAtomProperty(*args): return _openbabel.OBResidue_GetAtomProperty(*args)
    def GetResidueProperty(*args): return _openbabel.OBResidue_GetResidueProperty(*args)
    def IsHetAtom(*args): return _openbabel.OBResidue_IsHetAtom(*args)
    def IsResidueType(*args): return _openbabel.OBResidue_IsResidueType(*args)
    def BeginAtom(*args): return _openbabel.OBResidue_BeginAtom(*args)
    def NextAtom(*args): return _openbabel.OBResidue_NextAtom(*args)

class OBResiduePtr(OBResidue):
    def __init__(self, this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = OBResidue
_openbabel.OBResidue_swigregister(OBResiduePtr)

OB_4RING_ATOM = _openbabel.OB_4RING_ATOM
OB_3RING_ATOM = _openbabel.OB_3RING_ATOM
OB_AROMATIC_ATOM = _openbabel.OB_AROMATIC_ATOM
OB_RING_ATOM = _openbabel.OB_RING_ATOM
OB_CSTEREO_ATOM = _openbabel.OB_CSTEREO_ATOM
OB_ACSTEREO_ATOM = _openbabel.OB_ACSTEREO_ATOM
OB_DONOR_ATOM = _openbabel.OB_DONOR_ATOM
OB_ACCEPTOR_ATOM = _openbabel.OB_ACCEPTOR_ATOM
OB_CHIRAL_ATOM = _openbabel.OB_CHIRAL_ATOM
OB_POS_CHIRAL_ATOM = _openbabel.OB_POS_CHIRAL_ATOM
OB_NEG_CHIRAL_ATOM = _openbabel.OB_NEG_CHIRAL_ATOM
class OBAtom(object):
    def __repr__(self):
        return "<%s.%s; proxy of C++ OpenBabel::OBAtom instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        newobj = _openbabel.new_OBAtom(*args)
        self.this = newobj.this
        self.thisown = 1
        del newobj.thisown
    def __del__(self, destroy=_openbabel.delete_OBAtom):
        try:
            if self.thisown: destroy(self)
        except: pass

    def Clear(*args): return _openbabel.OBAtom_Clear(*args)
    def SetIdx(*args): return _openbabel.OBAtom_SetIdx(*args)
    def SetHyb(*args): return _openbabel.OBAtom_SetHyb(*args)
    def SetAtomicNum(*args): return _openbabel.OBAtom_SetAtomicNum(*args)
    def SetIsotope(*args): return _openbabel.OBAtom_SetIsotope(*args)
    def SetImplicitValence(*args): return _openbabel.OBAtom_SetImplicitValence(*args)
    def IncrementImplicitValence(*args): return _openbabel.OBAtom_IncrementImplicitValence(*args)
    def DecrementImplicitValence(*args): return _openbabel.OBAtom_DecrementImplicitValence(*args)
    def SetFormalCharge(*args): return _openbabel.OBAtom_SetFormalCharge(*args)
    def SetSpinMultiplicity(*args): return _openbabel.OBAtom_SetSpinMultiplicity(*args)
    def SetType(*args): return _openbabel.OBAtom_SetType(*args)
    def SetPartialCharge(*args): return _openbabel.OBAtom_SetPartialCharge(*args)
    def SetCoordPtr(*args): return _openbabel.OBAtom_SetCoordPtr(*args)
    def SetVector(*args): return _openbabel.OBAtom_SetVector(*args)
    def SetResidue(*args): return _openbabel.OBAtom_SetResidue(*args)
    def SetAromatic(*args): return _openbabel.OBAtom_SetAromatic(*args)
    def UnsetAromatic(*args): return _openbabel.OBAtom_UnsetAromatic(*args)
    def SetClockwiseStereo(*args): return _openbabel.OBAtom_SetClockwiseStereo(*args)
    def SetAntiClockwiseStereo(*args): return _openbabel.OBAtom_SetAntiClockwiseStereo(*args)
    def SetPositiveStereo(*args): return _openbabel.OBAtom_SetPositiveStereo(*args)
    def SetNegativeStereo(*args): return _openbabel.OBAtom_SetNegativeStereo(*args)
    def UnsetStereo(*args): return _openbabel.OBAtom_UnsetStereo(*args)
    def SetInRing(*args): return _openbabel.OBAtom_SetInRing(*args)
    def SetChiral(*args): return _openbabel.OBAtom_SetChiral(*args)
    def ClearCoordPtr(*args): return _openbabel.OBAtom_ClearCoordPtr(*args)
    def GetFormalCharge(*args): return _openbabel.OBAtom_GetFormalCharge(*args)
    def GetAtomicNum(*args): return _openbabel.OBAtom_GetAtomicNum(*args)
    def GetIsotope(*args): return _openbabel.OBAtom_GetIsotope(*args)
    def GetSpinMultiplicity(*args): return _openbabel.OBAtom_GetSpinMultiplicity(*args)
    def GetAtomicMass(*args): return _openbabel.OBAtom_GetAtomicMass(*args)
    def GetExactMass(*args): return _openbabel.OBAtom_GetExactMass(*args)
    def GetIdx(*args): return _openbabel.OBAtom_GetIdx(*args)
    def GetCoordinateIdx(*args): return _openbabel.OBAtom_GetCoordinateIdx(*args)
    def GetCIdx(*args): return _openbabel.OBAtom_GetCIdx(*args)
    def GetValence(*args): return _openbabel.OBAtom_GetValence(*args)
    def GetHyb(*args): return _openbabel.OBAtom_GetHyb(*args)
    def GetImplicitValence(*args): return _openbabel.OBAtom_GetImplicitValence(*args)
    def GetHvyValence(*args): return _openbabel.OBAtom_GetHvyValence(*args)
    def GetHeteroValence(*args): return _openbabel.OBAtom_GetHeteroValence(*args)
    def GetType(*args): return _openbabel.OBAtom_GetType(*args)
    def GetX(*args): return _openbabel.OBAtom_GetX(*args)
    def GetY(*args): return _openbabel.OBAtom_GetY(*args)
    def GetZ(*args): return _openbabel.OBAtom_GetZ(*args)
    def x(*args): return _openbabel.OBAtom_x(*args)
    def y(*args): return _openbabel.OBAtom_y(*args)
    def z(*args): return _openbabel.OBAtom_z(*args)
    def GetCoordinate(*args): return _openbabel.OBAtom_GetCoordinate(*args)
    def GetVector(*args): return _openbabel.OBAtom_GetVector(*args)
    def GetPartialCharge(*args): return _openbabel.OBAtom_GetPartialCharge(*args)
    def GetResidue(*args): return _openbabel.OBAtom_GetResidue(*args)
    def GetNewBondVector(*args): return _openbabel.OBAtom_GetNewBondVector(*args)
    def GetBond(*args): return _openbabel.OBAtom_GetBond(*args)
    def GetNextAtom(*args): return _openbabel.OBAtom_GetNextAtom(*args)
    def BeginBonds(*args): return _openbabel.OBAtom_BeginBonds(*args)
    def EndBonds(*args): return _openbabel.OBAtom_EndBonds(*args)
    def BeginBond(*args): return _openbabel.OBAtom_BeginBond(*args)
    def NextBond(*args): return _openbabel.OBAtom_NextBond(*args)
    def BeginNbrAtom(*args): return _openbabel.OBAtom_BeginNbrAtom(*args)
    def NextNbrAtom(*args): return _openbabel.OBAtom_NextNbrAtom(*args)
    def GetDistance(*args): return _openbabel.OBAtom_GetDistance(*args)
    def GetAngle(*args): return _openbabel.OBAtom_GetAngle(*args)
    def NewResidue(*args): return _openbabel.OBAtom_NewResidue(*args)
    def DeleteResidue(*args): return _openbabel.OBAtom_DeleteResidue(*args)
    def AddBond(*args): return _openbabel.OBAtom_AddBond(*args)
    def InsertBond(*args): return _openbabel.OBAtom_InsertBond(*args)
    def DeleteBond(*args): return _openbabel.OBAtom_DeleteBond(*args)
    def ClearBond(*args): return _openbabel.OBAtom_ClearBond(*args)
    def CountFreeOxygens(*args): return _openbabel.OBAtom_CountFreeOxygens(*args)
    def ImplicitHydrogenCount(*args): return _openbabel.OBAtom_ImplicitHydrogenCount(*args)
    def ExplicitHydrogenCount(*args): return _openbabel.OBAtom_ExplicitHydrogenCount(*args)
    def MemberOfRingCount(*args): return _openbabel.OBAtom_MemberOfRingCount(*args)
    def MemberOfRingSize(*args): return _openbabel.OBAtom_MemberOfRingSize(*args)
    def SmallestBondAngle(*args): return _openbabel.OBAtom_SmallestBondAngle(*args)
    def AverageBondAngle(*args): return _openbabel.OBAtom_AverageBondAngle(*args)
    def BOSum(*args): return _openbabel.OBAtom_BOSum(*args)
    def KBOSum(*args): return _openbabel.OBAtom_KBOSum(*args)
    def HtoMethyl(*args): return _openbabel.OBAtom_HtoMethyl(*args)
    def SetHybAndGeom(*args): return _openbabel.OBAtom_SetHybAndGeom(*args)
    def HasResidue(*args): return _openbabel.OBAtom_HasResidue(*args)
    def IsHydrogen(*args): return _openbabel.OBAtom_IsHydrogen(*args)
    def IsCarbon(*args): return _openbabel.OBAtom_IsCarbon(*args)
    def IsNitrogen(*args): return _openbabel.OBAtom_IsNitrogen(*args)
    def IsOxygen(*args): return _openbabel.OBAtom_IsOxygen(*args)
    def IsSulfur(*args): return _openbabel.OBAtom_IsSulfur(*args)
    def IsPhosphorus(*args): return _openbabel.OBAtom_IsPhosphorus(*args)
    def IsAromatic(*args): return _openbabel.OBAtom_IsAromatic(*args)
    def IsInRing(*args): return _openbabel.OBAtom_IsInRing(*args)
    def IsInRingSize(*args): return _openbabel.OBAtom_IsInRingSize(*args)
    def IsHeteroatom(*args): return _openbabel.OBAtom_IsHeteroatom(*args)
    def IsNotCorH(*args): return _openbabel.OBAtom_IsNotCorH(*args)
    def IsConnected(*args): return _openbabel.OBAtom_IsConnected(*args)
    def IsOneThree(*args): return _openbabel.OBAtom_IsOneThree(*args)
    def IsOneFour(*args): return _openbabel.OBAtom_IsOneFour(*args)
    def IsCarboxylOxygen(*args): return _openbabel.OBAtom_IsCarboxylOxygen(*args)
    def IsPhosphateOxygen(*args): return _openbabel.OBAtom_IsPhosphateOxygen(*args)
    def IsSulfateOxygen(*args): return _openbabel.OBAtom_IsSulfateOxygen(*args)
    def IsNitroOxygen(*args): return _openbabel.OBAtom_IsNitroOxygen(*args)
    def IsAmideNitrogen(*args): return _openbabel.OBAtom_IsAmideNitrogen(*args)
    def IsPolarHydrogen(*args): return _openbabel.OBAtom_IsPolarHydrogen(*args)
    def IsNonPolarHydrogen(*args): return _openbabel.OBAtom_IsNonPolarHydrogen(*args)
    def IsAromaticNOxide(*args): return _openbabel.OBAtom_IsAromaticNOxide(*args)
    def IsChiral(*args): return _openbabel.OBAtom_IsChiral(*args)
    def IsAxial(*args): return _openbabel.OBAtom_IsAxial(*args)
    def IsClockwise(*args): return _openbabel.OBAtom_IsClockwise(*args)
    def IsAntiClockwise(*args): return _openbabel.OBAtom_IsAntiClockwise(*args)
    def IsPositiveStereo(*args): return _openbabel.OBAtom_IsPositiveStereo(*args)
    def IsNegativeStereo(*args): return _openbabel.OBAtom_IsNegativeStereo(*args)
    def HasChiralitySpecified(*args): return _openbabel.OBAtom_HasChiralitySpecified(*args)
    def HasChiralVolume(*args): return _openbabel.OBAtom_HasChiralVolume(*args)
    def IsHbondAcceptor(*args): return _openbabel.OBAtom_IsHbondAcceptor(*args)
    def IsHbondDonor(*args): return _openbabel.OBAtom_IsHbondDonor(*args)
    def IsHbondDonorH(*args): return _openbabel.OBAtom_IsHbondDonorH(*args)
    def HasAlphaBetaUnsat(*args): return _openbabel.OBAtom_HasAlphaBetaUnsat(*args)
    def HasBondOfOrder(*args): return _openbabel.OBAtom_HasBondOfOrder(*args)
    def CountBondsOfOrder(*args): return _openbabel.OBAtom_CountBondsOfOrder(*args)
    def HasNonSingleBond(*args): return _openbabel.OBAtom_HasNonSingleBond(*args)
    def HasSingleBond(*args): return _openbabel.OBAtom_HasSingleBond(*args)
    def HasDoubleBond(*args): return _openbabel.OBAtom_HasDoubleBond(*args)
    def HasAromaticBond(*args): return _openbabel.OBAtom_HasAromaticBond(*args)
    def MatchesSMARTS(*args): return _openbabel.OBAtom_MatchesSMARTS(*args)

class OBAtomPtr(OBAtom):
    def __init__(self, this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = OBAtom
_openbabel.OBAtom_swigregister(OBAtomPtr)

OB_AROMATIC_BOND = _openbabel.OB_AROMATIC_BOND
OB_WEDGE_BOND = _openbabel.OB_WEDGE_BOND
OB_HASH_BOND = _openbabel.OB_HASH_BOND
OB_RING_BOND = _openbabel.OB_RING_BOND
OB_TORUP_BOND = _openbabel.OB_TORUP_BOND
OB_TORDOWN_BOND = _openbabel.OB_TORDOWN_BOND
OB_KSINGLE_BOND = _openbabel.OB_KSINGLE_BOND
OB_KDOUBLE_BOND = _openbabel.OB_KDOUBLE_BOND
OB_KTRIPLE_BOND = _openbabel.OB_KTRIPLE_BOND
OB_CLOSURE_BOND = _openbabel.OB_CLOSURE_BOND
class OBBond(object):
    def __repr__(self):
        return "<%s.%s; proxy of C++ OpenBabel::OBBond instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        newobj = _openbabel.new_OBBond(*args)
        self.this = newobj.this
        self.thisown = 1
        del newobj.thisown
    def __del__(self, destroy=_openbabel.delete_OBBond):
        try:
            if self.thisown: destroy(self)
        except: pass

    def SetIdx(*args): return _openbabel.OBBond_SetIdx(*args)
    def SetBO(*args): return _openbabel.OBBond_SetBO(*args)
    def SetBegin(*args): return _openbabel.OBBond_SetBegin(*args)
    def SetEnd(*args): return _openbabel.OBBond_SetEnd(*args)
    def SetLength(*args): return _openbabel.OBBond_SetLength(*args)
    def Set(*args): return _openbabel.OBBond_Set(*args)
    def SetKSingle(*args): return _openbabel.OBBond_SetKSingle(*args)
    def SetKDouble(*args): return _openbabel.OBBond_SetKDouble(*args)
    def SetKTriple(*args): return _openbabel.OBBond_SetKTriple(*args)
    def SetAromatic(*args): return _openbabel.OBBond_SetAromatic(*args)
    def SetHash(*args): return _openbabel.OBBond_SetHash(*args)
    def SetWedge(*args): return _openbabel.OBBond_SetWedge(*args)
    def SetUp(*args): return _openbabel.OBBond_SetUp(*args)
    def SetDown(*args): return _openbabel.OBBond_SetDown(*args)
    def SetInRing(*args): return _openbabel.OBBond_SetInRing(*args)
    def SetClosure(*args): return _openbabel.OBBond_SetClosure(*args)
    def UnsetAromatic(*args): return _openbabel.OBBond_UnsetAromatic(*args)
    def UnsetKekule(*args): return _openbabel.OBBond_UnsetKekule(*args)
    def GetBO(*args): return _openbabel.OBBond_GetBO(*args)
    def GetBondOrder(*args): return _openbabel.OBBond_GetBondOrder(*args)
    def GetFlags(*args): return _openbabel.OBBond_GetFlags(*args)
    def GetBeginAtomIdx(*args): return _openbabel.OBBond_GetBeginAtomIdx(*args)
    def GetEndAtomIdx(*args): return _openbabel.OBBond_GetEndAtomIdx(*args)
    def GetBeginAtom(*args): return _openbabel.OBBond_GetBeginAtom(*args)
    def GetEndAtom(*args): return _openbabel.OBBond_GetEndAtom(*args)
    def GetNbrAtom(*args): return _openbabel.OBBond_GetNbrAtom(*args)
    def GetEquibLength(*args): return _openbabel.OBBond_GetEquibLength(*args)
    def GetLength(*args): return _openbabel.OBBond_GetLength(*args)
    def GetNbrAtomIdx(*args): return _openbabel.OBBond_GetNbrAtomIdx(*args)
    def IsAromatic(*args): return _openbabel.OBBond_IsAromatic(*args)
    def IsInRing(*args): return _openbabel.OBBond_IsInRing(*args)
    def IsRotor(*args): return _openbabel.OBBond_IsRotor(*args)
    def IsAmide(*args): return _openbabel.OBBond_IsAmide(*args)
    def IsPrimaryAmide(*args): return _openbabel.OBBond_IsPrimaryAmide(*args)
    def IsSecondaryAmide(*args): return _openbabel.OBBond_IsSecondaryAmide(*args)
    def IsEster(*args): return _openbabel.OBBond_IsEster(*args)
    def IsCarbonyl(*args): return _openbabel.OBBond_IsCarbonyl(*args)
    def IsSingle(*args): return _openbabel.OBBond_IsSingle(*args)
    def IsDouble(*args): return _openbabel.OBBond_IsDouble(*args)
    def IsTriple(*args): return _openbabel.OBBond_IsTriple(*args)
    def IsKSingle(*args): return _openbabel.OBBond_IsKSingle(*args)
    def IsKDouble(*args): return _openbabel.OBBond_IsKDouble(*args)
    def IsKTriple(*args): return _openbabel.OBBond_IsKTriple(*args)
    def IsClosure(*args): return _openbabel.OBBond_IsClosure(*args)
    def IsUp(*args): return _openbabel.OBBond_IsUp(*args)
    def IsDown(*args): return _openbabel.OBBond_IsDown(*args)
    def IsWedge(*args): return _openbabel.OBBond_IsWedge(*args)
    def IsHash(*args): return _openbabel.OBBond_IsHash(*args)
    def IsDoubleBondGeometry(*args): return _openbabel.OBBond_IsDoubleBondGeometry(*args)

class OBBondPtr(OBBond):
    def __init__(self, this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = OBBond
_openbabel.OBBond_swigregister(OBBondPtr)

OB_SSSR_MOL = _openbabel.OB_SSSR_MOL
OB_RINGFLAGS_MOL = _openbabel.OB_RINGFLAGS_MOL
OB_AROMATIC_MOL = _openbabel.OB_AROMATIC_MOL
OB_ATOMTYPES_MOL = _openbabel.OB_ATOMTYPES_MOL
OB_CHIRALITY_MOL = _openbabel.OB_CHIRALITY_MOL
OB_PCHARGE_MOL = _openbabel.OB_PCHARGE_MOL
OB_HYBRID_MOL = _openbabel.OB_HYBRID_MOL
OB_IMPVAL_MOL = _openbabel.OB_IMPVAL_MOL
OB_KEKULE_MOL = _openbabel.OB_KEKULE_MOL
OB_CLOSURE_MOL = _openbabel.OB_CLOSURE_MOL
OB_H_ADDED_MOL = _openbabel.OB_H_ADDED_MOL
OB_PH_CORRECTED_MOL = _openbabel.OB_PH_CORRECTED_MOL
OB_AROM_CORRECTED_MOL = _openbabel.OB_AROM_CORRECTED_MOL
OB_CHAINS_MOL = _openbabel.OB_CHAINS_MOL
OB_TCHARGE_MOL = _openbabel.OB_TCHARGE_MOL
OB_TSPIN_MOL = _openbabel.OB_TSPIN_MOL
OB_CURRENT_CONFORMER = _openbabel.OB_CURRENT_CONFORMER
class OBMol(object):
    def __repr__(self):
        return "<%s.%s; proxy of C++ OpenBabel::OBMol instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        newobj = _openbabel.new_OBMol(*args)
        self.this = newobj.this
        self.thisown = 1
        del newobj.thisown
    def __del__(self, destroy=_openbabel.delete_OBMol):
        try:
            if self.thisown: destroy(self)
        except: pass

    def __iadd__(*args): return _openbabel.OBMol___iadd__(*args)
    def ReserveAtoms(*args): return _openbabel.OBMol_ReserveAtoms(*args)
    def CreateAtom(*args): return _openbabel.OBMol_CreateAtom(*args)
    def CreateBond(*args): return _openbabel.OBMol_CreateBond(*args)
    def DestroyAtom(*args): return _openbabel.OBMol_DestroyAtom(*args)
    def DestroyBond(*args): return _openbabel.OBMol_DestroyBond(*args)
    def AddAtom(*args): return _openbabel.OBMol_AddAtom(*args)
    def AddBond(*args): return _openbabel.OBMol_AddBond(*args)
    def AddResidue(*args): return _openbabel.OBMol_AddResidue(*args)
    def InsertAtom(*args): return _openbabel.OBMol_InsertAtom(*args)
    def DeleteAtom(*args): return _openbabel.OBMol_DeleteAtom(*args)
    def DeleteBond(*args): return _openbabel.OBMol_DeleteBond(*args)
    def DeleteResidue(*args): return _openbabel.OBMol_DeleteResidue(*args)
    def NewAtom(*args): return _openbabel.OBMol_NewAtom(*args)
    def NewResidue(*args): return _openbabel.OBMol_NewResidue(*args)
    def BeginModify(*args): return _openbabel.OBMol_BeginModify(*args)
    def EndModify(*args): return _openbabel.OBMol_EndModify(*args)
    def GetMod(*args): return _openbabel.OBMol_GetMod(*args)
    def IncrementMod(*args): return _openbabel.OBMol_IncrementMod(*args)
    def DecrementMod(*args): return _openbabel.OBMol_DecrementMod(*args)
    def GetFlags(*args): return _openbabel.OBMol_GetFlags(*args)
    def GetTitle(*args): return _openbabel.OBMol_GetTitle(*args)
    def NumAtoms(*args): return _openbabel.OBMol_NumAtoms(*args)
    def NumBonds(*args): return _openbabel.OBMol_NumBonds(*args)
    def NumHvyAtoms(*args): return _openbabel.OBMol_NumHvyAtoms(*args)
    def NumResidues(*args): return _openbabel.OBMol_NumResidues(*args)
    def NumRotors(*args): return _openbabel.OBMol_NumRotors(*args)
    def GetAtom(*args): return _openbabel.OBMol_GetAtom(*args)
    def GetFirstAtom(*args): return _openbabel.OBMol_GetFirstAtom(*args)
    def GetBond(*args): return _openbabel.OBMol_GetBond(*args)
    def GetResidue(*args): return _openbabel.OBMol_GetResidue(*args)
    def GetInternalCoord(*args): return _openbabel.OBMol_GetInternalCoord(*args)
    def GetTorsion(*args): return _openbabel.OBMol_GetTorsion(*args)
    def GetFormula(*args): return _openbabel.OBMol_GetFormula(*args)
    def GetEnergy(*args): return _openbabel.OBMol_GetEnergy(*args)
    def GetMolWt(*args): return _openbabel.OBMol_GetMolWt(*args)
    def GetExactMass(*args): return _openbabel.OBMol_GetExactMass(*args)
    def GetTotalCharge(*args): return _openbabel.OBMol_GetTotalCharge(*args)
    def GetTotalSpinMultiplicity(*args): return _openbabel.OBMol_GetTotalSpinMultiplicity(*args)
    def GetDimension(*args): return _openbabel.OBMol_GetDimension(*args)
    def GetCoordinates(*args): return _openbabel.OBMol_GetCoordinates(*args)
    def GetSSSR(*args): return _openbabel.OBMol_GetSSSR(*args)
    def AutomaticFormalCharge(*args): return _openbabel.OBMol_AutomaticFormalCharge(*args)
    def AutomaticPartialCharge(*args): return _openbabel.OBMol_AutomaticPartialCharge(*args)
    def SetTitle(*args): return _openbabel.OBMol_SetTitle(*args)
    def SetFormula(*args): return _openbabel.OBMol_SetFormula(*args)
    def SetEnergy(*args): return _openbabel.OBMol_SetEnergy(*args)
    def SetDimension(*args): return _openbabel.OBMol_SetDimension(*args)
    def SetTotalCharge(*args): return _openbabel.OBMol_SetTotalCharge(*args)
    def SetTotalSpinMultiplicity(*args): return _openbabel.OBMol_SetTotalSpinMultiplicity(*args)
    def SetInternalCoord(*args): return _openbabel.OBMol_SetInternalCoord(*args)
    def SetAutomaticFormalCharge(*args): return _openbabel.OBMol_SetAutomaticFormalCharge(*args)
    def SetAutomaticPartialCharge(*args): return _openbabel.OBMol_SetAutomaticPartialCharge(*args)
    def SetAromaticPerceived(*args): return _openbabel.OBMol_SetAromaticPerceived(*args)
    def SetSSSRPerceived(*args): return _openbabel.OBMol_SetSSSRPerceived(*args)
    def SetRingAtomsAndBondsPerceived(*args): return _openbabel.OBMol_SetRingAtomsAndBondsPerceived(*args)
    def SetAtomTypesPerceived(*args): return _openbabel.OBMol_SetAtomTypesPerceived(*args)
    def SetChainsPerceived(*args): return _openbabel.OBMol_SetChainsPerceived(*args)
    def SetChiralityPerceived(*args): return _openbabel.OBMol_SetChiralityPerceived(*args)
    def SetPartialChargesPerceived(*args): return _openbabel.OBMol_SetPartialChargesPerceived(*args)
    def SetHybridizationPerceived(*args): return _openbabel.OBMol_SetHybridizationPerceived(*args)
    def SetImplicitValencePerceived(*args): return _openbabel.OBMol_SetImplicitValencePerceived(*args)
    def SetKekulePerceived(*args): return _openbabel.OBMol_SetKekulePerceived(*args)
    def SetClosureBondsPerceived(*args): return _openbabel.OBMol_SetClosureBondsPerceived(*args)
    def SetHydrogensAdded(*args): return _openbabel.OBMol_SetHydrogensAdded(*args)
    def SetCorrectedForPH(*args): return _openbabel.OBMol_SetCorrectedForPH(*args)
    def SetAromaticCorrected(*args): return _openbabel.OBMol_SetAromaticCorrected(*args)
    def SetSpinMultiplicityAssigned(*args): return _openbabel.OBMol_SetSpinMultiplicityAssigned(*args)
    def SetFlags(*args): return _openbabel.OBMol_SetFlags(*args)
    def UnsetAromaticPerceived(*args): return _openbabel.OBMol_UnsetAromaticPerceived(*args)
    def UnsetPartialChargesPerceived(*args): return _openbabel.OBMol_UnsetPartialChargesPerceived(*args)
    def UnsetImplicitValencePerceived(*args): return _openbabel.OBMol_UnsetImplicitValencePerceived(*args)
    def UnsetFlag(*args): return _openbabel.OBMol_UnsetFlag(*args)
    def DoTransformations(*args): return _openbabel.OBMol_DoTransformations(*args)
    ClassDescription = staticmethod(_openbabel.OBMol_ClassDescription)
    def Clear(*args): return _openbabel.OBMol_Clear(*args)
    def RenumberAtoms(*args): return _openbabel.OBMol_RenumberAtoms(*args)
    def ToInertialFrame(*args): return _openbabel.OBMol_ToInertialFrame(*args)
    def Translate(*args): return _openbabel.OBMol_Translate(*args)
    def Rotate(*args): return _openbabel.OBMol_Rotate(*args)
    def Kekulize(*args): return _openbabel.OBMol_Kekulize(*args)
    def PerceiveKekuleBonds(*args): return _openbabel.OBMol_PerceiveKekuleBonds(*args)
    def NewPerceiveKekuleBonds(*args): return _openbabel.OBMol_NewPerceiveKekuleBonds(*args)
    def DeleteHydrogen(*args): return _openbabel.OBMol_DeleteHydrogen(*args)
    def DeleteHydrogens(*args): return _openbabel.OBMol_DeleteHydrogens(*args)
    def DeleteNonPolarHydrogens(*args): return _openbabel.OBMol_DeleteNonPolarHydrogens(*args)
    def AddHydrogens(*args): return _openbabel.OBMol_AddHydrogens(*args)
    def AddPolarHydrogens(*args): return _openbabel.OBMol_AddPolarHydrogens(*args)
    def StripSalts(*args): return _openbabel.OBMol_StripSalts(*args)
    def ConvertDativeBonds(*args): return _openbabel.OBMol_ConvertDativeBonds(*args)
    def CorrectForPH(*args): return _openbabel.OBMol_CorrectForPH(*args)
    def AssignSpinMultiplicity(*args): return _openbabel.OBMol_AssignSpinMultiplicity(*args)
    def Center(*args): return _openbabel.OBMol_Center(*args)
    def SetTorsion(*args): return _openbabel.OBMol_SetTorsion(*args)
    def FindSSSR(*args): return _openbabel.OBMol_FindSSSR(*args)
    def FindRingAtomsAndBonds(*args): return _openbabel.OBMol_FindRingAtomsAndBonds(*args)
    def FindChiralCenters(*args): return _openbabel.OBMol_FindChiralCenters(*args)
    def FindChildren(*args): return _openbabel.OBMol_FindChildren(*args)
    def FindLargestFragment(*args): return _openbabel.OBMol_FindLargestFragment(*args)
    def ContigFragList(*args): return _openbabel.OBMol_ContigFragList(*args)
    def Align(*args): return _openbabel.OBMol_Align(*args)
    def ConnectTheDots(*args): return _openbabel.OBMol_ConnectTheDots(*args)
    def PerceiveBondOrders(*args): return _openbabel.OBMol_PerceiveBondOrders(*args)
    def FindTorsions(*args): return _openbabel.OBMol_FindTorsions(*args)
    def GetGTDVector(*args): return _openbabel.OBMol_GetGTDVector(*args)
    def GetGIVector(*args): return _openbabel.OBMol_GetGIVector(*args)
    def GetGIDVector(*args): return _openbabel.OBMol_GetGIDVector(*args)
    def Has2D(*args): return _openbabel.OBMol_Has2D(*args)
    def Has3D(*args): return _openbabel.OBMol_Has3D(*args)
    def HasNonZeroCoords(*args): return _openbabel.OBMol_HasNonZeroCoords(*args)
    def HasAromaticPerceived(*args): return _openbabel.OBMol_HasAromaticPerceived(*args)
    def HasSSSRPerceived(*args): return _openbabel.OBMol_HasSSSRPerceived(*args)
    def HasRingAtomsAndBondsPerceived(*args): return _openbabel.OBMol_HasRingAtomsAndBondsPerceived(*args)
    def HasAtomTypesPerceived(*args): return _openbabel.OBMol_HasAtomTypesPerceived(*args)
    def HasChiralityPerceived(*args): return _openbabel.OBMol_HasChiralityPerceived(*args)
    def HasPartialChargesPerceived(*args): return _openbabel.OBMol_HasPartialChargesPerceived(*args)
    def HasHybridizationPerceived(*args): return _openbabel.OBMol_HasHybridizationPerceived(*args)
    def HasImplicitValencePerceived(*args): return _openbabel.OBMol_HasImplicitValencePerceived(*args)
    def HasKekulePerceived(*args): return _openbabel.OBMol_HasKekulePerceived(*args)
    def HasClosureBondsPerceived(*args): return _openbabel.OBMol_HasClosureBondsPerceived(*args)
    def HasChainsPerceived(*args): return _openbabel.OBMol_HasChainsPerceived(*args)
    def HasHydrogensAdded(*args): return _openbabel.OBMol_HasHydrogensAdded(*args)
    def HasAromaticCorrected(*args): return _openbabel.OBMol_HasAromaticCorrected(*args)
    def IsCorrectedForPH(*args): return _openbabel.OBMol_IsCorrectedForPH(*args)
    def HasSpinMultiplicityAssigned(*args): return _openbabel.OBMol_HasSpinMultiplicityAssigned(*args)
    def IsChiral(*args): return _openbabel.OBMol_IsChiral(*args)
    def Empty(*args): return _openbabel.OBMol_Empty(*args)
    def NumConformers(*args): return _openbabel.OBMol_NumConformers(*args)
    def SetConformers(*args): return _openbabel.OBMol_SetConformers(*args)
    def AddConformer(*args): return _openbabel.OBMol_AddConformer(*args)
    def SetConformer(*args): return _openbabel.OBMol_SetConformer(*args)
    def CopyConformer(*args): return _openbabel.OBMol_CopyConformer(*args)
    def DeleteConformer(*args): return _openbabel.OBMol_DeleteConformer(*args)
    def GetConformer(*args): return _openbabel.OBMol_GetConformer(*args)
    def BeginConformer(*args): return _openbabel.OBMol_BeginConformer(*args)
    def NextConformer(*args): return _openbabel.OBMol_NextConformer(*args)
    def GetConformers(*args): return _openbabel.OBMol_GetConformers(*args)
    def BeginAtom(*args): return _openbabel.OBMol_BeginAtom(*args)
    def NextAtom(*args): return _openbabel.OBMol_NextAtom(*args)
    def BeginBond(*args): return _openbabel.OBMol_BeginBond(*args)
    def NextBond(*args): return _openbabel.OBMol_NextBond(*args)
    def BeginResidue(*args): return _openbabel.OBMol_BeginResidue(*args)
    def NextResidue(*args): return _openbabel.OBMol_NextResidue(*args)
    def BeginInternalCoord(*args): return _openbabel.OBMol_BeginInternalCoord(*args)
    def NextInternalCoord(*args): return _openbabel.OBMol_NextInternalCoord(*args)

class OBMolPtr(OBMol):
    def __init__(self, this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = OBMol
_openbabel.OBMol_swigregister(OBMolPtr)

OBMol_ClassDescription = _openbabel.OBMol_ClassDescription

class OBInternalCoord(object):
    def __repr__(self):
        return "<%s.%s; proxy of C++ OpenBabel::OBInternalCoord instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    _a = property(_openbabel.OBInternalCoord__a_get, _openbabel.OBInternalCoord__a_set)
    _b = property(_openbabel.OBInternalCoord__b_get, _openbabel.OBInternalCoord__b_set)
    _c = property(_openbabel.OBInternalCoord__c_get, _openbabel.OBInternalCoord__c_set)
    _dst = property(_openbabel.OBInternalCoord__dst_get, _openbabel.OBInternalCoord__dst_set)
    _ang = property(_openbabel.OBInternalCoord__ang_get, _openbabel.OBInternalCoord__ang_set)
    _tor = property(_openbabel.OBInternalCoord__tor_get, _openbabel.OBInternalCoord__tor_set)
    def __init__(self, *args):
        newobj = _openbabel.new_OBInternalCoord(*args)
        self.this = newobj.this
        self.thisown = 1
        del newobj.thisown
    def __del__(self, destroy=_openbabel.delete_OBInternalCoord):
        try:
            if self.thisown: destroy(self)
        except: pass


class OBInternalCoordPtr(OBInternalCoord):
    def __init__(self, this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = OBInternalCoord
_openbabel.OBInternalCoord_swigregister(OBInternalCoordPtr)


CartesianToInternal = _openbabel.CartesianToInternal

InternalToCartesian = _openbabel.InternalToCartesian

NewExtension = _openbabel.NewExtension
BUFF_SIZE = _openbabel.BUFF_SIZE

get_rmat = _openbabel.get_rmat

ob_make_rmat = _openbabel.ob_make_rmat

qtrfit = _openbabel.qtrfit

superimpose = _openbabel.superimpose


tokenize = _openbabel.tokenize

ThrowError = _openbabel.ThrowError
cvar = _openbabel.cvar

