# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

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


class OBResidue(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, OBResidue, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, OBResidue, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ OpenBabel::OBResidue instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, OBResidue, 'this', _openbabel.new_OBResidue(*args))
        _swig_setattr(self, OBResidue, 'thisown', 1)
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
    def HasData(*args): return _openbabel.OBResidue_HasData(*args)
    def DeleteData(*args): return _openbabel.OBResidue_DeleteData(*args)
    def SetData(*args): return _openbabel.OBResidue_SetData(*args)
    def DataSize(*args): return _openbabel.OBResidue_DataSize(*args)
    def GetData(*args): return _openbabel.OBResidue_GetData(*args)
    def BeginData(*args): return _openbabel.OBResidue_BeginData(*args)
    def EndData(*args): return _openbabel.OBResidue_EndData(*args)

class OBResiduePtr(OBResidue):
    def __init__(self, this):
        _swig_setattr(self, OBResidue, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, OBResidue, 'thisown', 0)
        _swig_setattr(self, OBResidue,self.__class__,OBResidue)
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
class OBAtom(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, OBAtom, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, OBAtom, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ OpenBabel::OBAtom instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, OBAtom, 'this', _openbabel.new_OBAtom(*args))
        _swig_setattr(self, OBAtom, 'thisown', 1)
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
    def SetVector(*args): return _openbabel.OBAtom_SetVector(*args)
    def SetResidue(*args): return _openbabel.OBAtom_SetResidue(*args)
    def SetCoordPtr(*args): return _openbabel.OBAtom_SetCoordPtr(*args)
    def SetAromatic(*args): return _openbabel.OBAtom_SetAromatic(*args)
    def UnsetAromatic(*args): return _openbabel.OBAtom_UnsetAromatic(*args)
    def SetClockwiseStereo(*args): return _openbabel.OBAtom_SetClockwiseStereo(*args)
    def SetAntiClockwiseStereo(*args): return _openbabel.OBAtom_SetAntiClockwiseStereo(*args)
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
    def x(*args): return _openbabel.OBAtom_x(*args)
    def GetY(*args): return _openbabel.OBAtom_GetY(*args)
    def y(*args): return _openbabel.OBAtom_y(*args)
    def GetZ(*args): return _openbabel.OBAtom_GetZ(*args)
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
    def HasChiralitySpecified(*args): return _openbabel.OBAtom_HasChiralitySpecified(*args)
    def HasAlphaBetaUnsat(*args): return _openbabel.OBAtom_HasAlphaBetaUnsat(*args)
    def HasBondOfOrder(*args): return _openbabel.OBAtom_HasBondOfOrder(*args)
    def CountBondsOfOrder(*args): return _openbabel.OBAtom_CountBondsOfOrder(*args)
    def HasNonSingleBond(*args): return _openbabel.OBAtom_HasNonSingleBond(*args)
    def HasSingleBond(*args): return _openbabel.OBAtom_HasSingleBond(*args)
    def HasDoubleBond(*args): return _openbabel.OBAtom_HasDoubleBond(*args)
    def HasAromaticBond(*args): return _openbabel.OBAtom_HasAromaticBond(*args)
    def MatchesSMARTS(*args): return _openbabel.OBAtom_MatchesSMARTS(*args)
    def HasData(*args): return _openbabel.OBAtom_HasData(*args)
    def DeleteData(*args): return _openbabel.OBAtom_DeleteData(*args)
    def SetData(*args): return _openbabel.OBAtom_SetData(*args)
    def DataSize(*args): return _openbabel.OBAtom_DataSize(*args)
    def GetData(*args): return _openbabel.OBAtom_GetData(*args)
    def BeginData(*args): return _openbabel.OBAtom_BeginData(*args)
    def EndData(*args): return _openbabel.OBAtom_EndData(*args)

class OBAtomPtr(OBAtom):
    def __init__(self, this):
        _swig_setattr(self, OBAtom, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, OBAtom, 'thisown', 0)
        _swig_setattr(self, OBAtom,self.__class__,OBAtom)
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
class OBBond(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, OBBond, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, OBBond, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ OpenBabel::OBBond instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, OBBond, 'this', _openbabel.new_OBBond(*args))
        _swig_setattr(self, OBBond, 'thisown', 1)
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
    def HasData(*args): return _openbabel.OBBond_HasData(*args)
    def DeleteData(*args): return _openbabel.OBBond_DeleteData(*args)
    def SetData(*args): return _openbabel.OBBond_SetData(*args)
    def DataSize(*args): return _openbabel.OBBond_DataSize(*args)
    def GetData(*args): return _openbabel.OBBond_GetData(*args)
    def BeginData(*args): return _openbabel.OBBond_BeginData(*args)
    def EndData(*args): return _openbabel.OBBond_EndData(*args)

class OBBondPtr(OBBond):
    def __init__(self, this):
        _swig_setattr(self, OBBond, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, OBBond, 'thisown', 0)
        _swig_setattr(self, OBBond,self.__class__,OBBond)
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
class OBMol(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, OBMol, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, OBMol, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ OpenBabel::OBMol instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, OBMol, 'this', _openbabel.new_OBMol(*args))
        _swig_setattr(self, OBMol, 'thisown', 1)
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
    def Compress(*args): return _openbabel.OBMol_Compress(*args)
    def UnCompress(*args): return _openbabel.OBMol_UnCompress(*args)
    def BeginAccess(*args): return _openbabel.OBMol_BeginAccess(*args)
    def EndAccess(*args): return _openbabel.OBMol_EndAccess(*args)
    def HasData(*args): return _openbabel.OBMol_HasData(*args)
    def DeleteData(*args): return _openbabel.OBMol_DeleteData(*args)
    def SetData(*args): return _openbabel.OBMol_SetData(*args)
    def DataSize(*args): return _openbabel.OBMol_DataSize(*args)
    def GetData(*args): return _openbabel.OBMol_GetData(*args)
    def BeginData(*args): return _openbabel.OBMol_BeginData(*args)
    def EndData(*args): return _openbabel.OBMol_EndData(*args)
    def GetFlags(*args): return _openbabel.OBMol_GetFlags(*args)
    def GetTitle(*args): return _openbabel.OBMol_GetTitle(*args)
    def GetInputType(*args): return _openbabel.OBMol_GetInputType(*args)
    def GetOutputType(*args): return _openbabel.OBMol_GetOutputType(*args)
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
    def GetEnergy(*args): return _openbabel.OBMol_GetEnergy(*args)
    def GetMolWt(*args): return _openbabel.OBMol_GetMolWt(*args)
    def GetExactMass(*args): return _openbabel.OBMol_GetExactMass(*args)
    def GetTotalCharge(*args): return _openbabel.OBMol_GetTotalCharge(*args)
    def GetTotalSpinMultiplicity(*args): return _openbabel.OBMol_GetTotalSpinMultiplicity(*args)
    def GetCoordinates(*args): return _openbabel.OBMol_GetCoordinates(*args)
    def GetSSSR(*args): return _openbabel.OBMol_GetSSSR(*args)
    def IsCompressed(*args): return _openbabel.OBMol_IsCompressed(*args)
    def AutomaticFormalCharge(*args): return _openbabel.OBMol_AutomaticFormalCharge(*args)
    def AutomaticPartialCharge(*args): return _openbabel.OBMol_AutomaticPartialCharge(*args)
    def SetTitle(*args): return _openbabel.OBMol_SetTitle(*args)
    def SetEnergy(*args): return _openbabel.OBMol_SetEnergy(*args)
    def SetTotalCharge(*args): return _openbabel.OBMol_SetTotalCharge(*args)
    def SetTotalSpinMultiplicity(*args): return _openbabel.OBMol_SetTotalSpinMultiplicity(*args)
    def SetInputType(*args): return _openbabel.OBMol_SetInputType(*args)
    def SetOutputType(*args): return _openbabel.OBMol_SetOutputType(*args)
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
    def UnsetAromaticPerceived(*args): return _openbabel.OBMol_UnsetAromaticPerceived(*args)
    def UnsetPartialChargesPerceived(*args): return _openbabel.OBMol_UnsetPartialChargesPerceived(*args)
    def UnsetImplicitValencePerceived(*args): return _openbabel.OBMol_UnsetImplicitValencePerceived(*args)
    def UnsetFlag(*args): return _openbabel.OBMol_UnsetFlag(*args)
    def SetFlags(*args): return _openbabel.OBMol_SetFlags(*args)
    def Clear(*args): return _openbabel.OBMol_Clear(*args)
    def RenumberAtoms(*args): return _openbabel.OBMol_RenumberAtoms(*args)
    def ToInertialFrame(*args): return _openbabel.OBMol_ToInertialFrame(*args)
    def Translate(*args): return _openbabel.OBMol_Translate(*args)
    def Rotate(*args): return _openbabel.OBMol_Rotate(*args)
    def Kekulize(*args): return _openbabel.OBMol_Kekulize(*args)
    def PerceiveKekuleBonds(*args): return _openbabel.OBMol_PerceiveKekuleBonds(*args)
    def NewPerceiveKekuleBonds(*args): return _openbabel.OBMol_NewPerceiveKekuleBonds(*args)
    def start_kekulize(*args): return _openbabel.OBMol_start_kekulize(*args)
    def expand_kekulize(*args): return _openbabel.OBMol_expand_kekulize(*args)
    def getorden(*args): return _openbabel.OBMol_getorden(*args)
    def expandcycle(*args): return _openbabel.OBMol_expandcycle(*args)
    def DeleteHydrogen(*args): return _openbabel.OBMol_DeleteHydrogen(*args)
    def DeleteHydrogens(*args): return _openbabel.OBMol_DeleteHydrogens(*args)
    def DeleteNonPolarHydrogens(*args): return _openbabel.OBMol_DeleteNonPolarHydrogens(*args)
    def AddHydrogens(*args): return _openbabel.OBMol_AddHydrogens(*args)
    def AddPolarHydrogens(*args): return _openbabel.OBMol_AddPolarHydrogens(*args)
    def StripSalts(*args): return _openbabel.OBMol_StripSalts(*args)
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
    def BeginAtom(*args): return _openbabel.OBMol_BeginAtom(*args)
    def NextAtom(*args): return _openbabel.OBMol_NextAtom(*args)
    def BeginBond(*args): return _openbabel.OBMol_BeginBond(*args)
    def NextBond(*args): return _openbabel.OBMol_NextBond(*args)
    def BeginResidue(*args): return _openbabel.OBMol_BeginResidue(*args)
    def NextResidue(*args): return _openbabel.OBMol_NextResidue(*args)
    def BeginInternalCoord(*args): return _openbabel.OBMol_BeginInternalCoord(*args)
    def NextInternalCoord(*args): return _openbabel.OBMol_NextInternalCoord(*args)
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
    def AssignResidueBonds(*args): return _openbabel.OBMol_AssignResidueBonds(*args)

class OBMolPtr(OBMol):
    def __init__(self, this):
        _swig_setattr(self, OBMol, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, OBMol, 'thisown', 0)
        _swig_setattr(self, OBMol,self.__class__,OBMol)
_openbabel.OBMol_swigregister(OBMolPtr)

class OBInternalCoord(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, OBInternalCoord, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, OBInternalCoord, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ OpenBabel::OBInternalCoord instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["_a"] = _openbabel.OBInternalCoord__a_set
    __swig_getmethods__["_a"] = _openbabel.OBInternalCoord__a_get
    if _newclass:_a = property(_openbabel.OBInternalCoord__a_get, _openbabel.OBInternalCoord__a_set)
    __swig_setmethods__["_b"] = _openbabel.OBInternalCoord__b_set
    __swig_getmethods__["_b"] = _openbabel.OBInternalCoord__b_get
    if _newclass:_b = property(_openbabel.OBInternalCoord__b_get, _openbabel.OBInternalCoord__b_set)
    __swig_setmethods__["_c"] = _openbabel.OBInternalCoord__c_set
    __swig_getmethods__["_c"] = _openbabel.OBInternalCoord__c_get
    if _newclass:_c = property(_openbabel.OBInternalCoord__c_get, _openbabel.OBInternalCoord__c_set)
    __swig_setmethods__["_dst"] = _openbabel.OBInternalCoord__dst_set
    __swig_getmethods__["_dst"] = _openbabel.OBInternalCoord__dst_get
    if _newclass:_dst = property(_openbabel.OBInternalCoord__dst_get, _openbabel.OBInternalCoord__dst_set)
    __swig_setmethods__["_ang"] = _openbabel.OBInternalCoord__ang_set
    __swig_getmethods__["_ang"] = _openbabel.OBInternalCoord__ang_get
    if _newclass:_ang = property(_openbabel.OBInternalCoord__ang_get, _openbabel.OBInternalCoord__ang_set)
    __swig_setmethods__["_tor"] = _openbabel.OBInternalCoord__tor_set
    __swig_getmethods__["_tor"] = _openbabel.OBInternalCoord__tor_get
    if _newclass:_tor = property(_openbabel.OBInternalCoord__tor_get, _openbabel.OBInternalCoord__tor_set)
    def __init__(self, *args):
        _swig_setattr(self, OBInternalCoord, 'this', _openbabel.new_OBInternalCoord(*args))
        _swig_setattr(self, OBInternalCoord, 'thisown', 1)
    def __del__(self, destroy=_openbabel.delete_OBInternalCoord):
        try:
            if self.thisown: destroy(self)
        except: pass


class OBInternalCoordPtr(OBInternalCoord):
    def __init__(self, this):
        _swig_setattr(self, OBInternalCoord, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, OBInternalCoord, 'thisown', 0)
        _swig_setattr(self, OBInternalCoord,self.__class__,OBInternalCoord)
_openbabel.OBInternalCoord_swigregister(OBInternalCoordPtr)


CartesianToInternal = _openbabel.CartesianToInternal

InternalToCartesian = _openbabel.InternalToCartesian

NewExtension = _openbabel.NewExtension

SetInputType = _openbabel.SetInputType

SetOutputType = _openbabel.SetOutputType
BUFF_SIZE = _openbabel.BUFF_SIZE

rint = _openbabel.rint

snprintf = _openbabel.snprintf

strncasecmp = _openbabel.strncasecmp

get_rmat = _openbabel.get_rmat

ob_make_rmat = _openbabel.ob_make_rmat

qtrfit = _openbabel.qtrfit

superimpose = _openbabel.superimpose

tokenize = _openbabel.tokenize

ThrowError = _openbabel.ThrowError
cvar = _openbabel.cvar

