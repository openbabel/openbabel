
%csconst(1);
//Rename GetType methods. GetType() is a member of the base C# object
//and should retain its default meaning in all derived classes
%rename(GetAtomType) OpenBabel::OBAtom::GetType();
%rename(GetBondType) OpenBabel::OBBond::GetType();
%rename(GetFormatType) OpenBabel::OBFormat::GetType();
%rename(GetRingType) OpenBabel::OBRing::GetType();
%rename(SetAtomType) OpenBabel::OBAtom::SetType(const char *);
%rename(SetBondType) OpenBabel::OBBond::SetType(const char *);
%rename(SetFormatType) OpenBabel::OBFormat::SetType(const char *);
%rename(SetRingType) OpenBabel::OBRing::SetType(const char *);

//renamed to avoid confusion with similarly named types in
//the System.Windows.Media.Media3d namespace and because all C# class names should start with a capital
//letter
%rename(OBVector3) *::vector3;
%rename(OBMatrix3x3) OpenBabel::matrix3x3;
%rename(OBTransform3d) OpenBabel::transform3d;

//renamed these because all public methods of a C# class should start with
//a capital letter. As swig for c# matures this may become uneccesary.
%rename(DistSq) OpenBabel::vector3::distSq(const vector3 &) const;
%rename(RandomUnitVector) OpenBabel::vector3::randomUnitVector(OBRandom *);
%rename(RandomUnitVector) OpenBabel::vector3::randomUnitVector();
%rename(Normalize) OpenBabel::vector3::normalize();
//changed this name slightly to match DistSq(vector3)
%rename(LengthSq) OpenBabel::vector3::length_2() const;
%rename(Length) OpenBabel::vector3::length() const;
%rename(CreateOrthoVector) OpenBabel::vector3::createOrthoVector(vector3 &) const;

//marks certain classes as partial for people who want to add methods and properties
//when using OBDotNet to build their own apis (ChemSharp will provide an example of this)
//this change can be made to other classes if needed in the future.
%typemap(csclassmodifiers) OpenBabel::OBAtom "public partial class"
%typemap(csclassmodifiers) OpenBabel::OBMol "public partial class"
%typemap(csclassmodifiers) OpenBabel::OBBond "public partial class"
%typemap(csclassmodifiers) OpenBabel::matrix3x3 "public partial class"
%typemap(csclassmodifiers) OpenBabel::OBRing "public partial class"
%typemap(csclassmodifiers) OpenBabel::vector3 "public partial class"
%typemap(csclassmodifiers) OpenBabel::OBElementTable "public partial class"

//adds the ignored operators back into the OBVector3 class
//with companion methods for use by .Net languages that don't support operator
//overloading
%typemap(cscode) OpenBabel::vector3
%{
	public static OBVector3 Add(OBVector3 vecA, OBVector3 vecB)
	{
		return new OBVector3(vecA.x()+vecB.x(), vecA.y()+vecB.y(),vecA.z()+vecB.z());
	}
	
	public static OBVector3 operator +(OBVector3 vecA, OBVector3 vecB)
	{
		return new OBVector3(vecA.x()+vecB.x(), vecA.y()+vecB.y(),vecA.z()+vecB.z());
	}  
	
	public static OBVector3 Mul(double d, OBVector3 vec)
	{
		return new OBVector3(d*vec.x(),d*vec.y(),d*vec.z());
	}
	
	public static OBVector3 operator *(double d, OBVector3 vec)
	{
		return new OBVector3(d*vec.x(),d*vec.y(),d*vec.z());
	}

	public static OBVector3 operator *(OBVector3 vec, double d)
	{
		return new OBVector3(d*vec.x(),d*vec.y(),d*vec.z());
	}
		
	public static OBVector3 Sub(OBVector3 vecA, OBVector3 vecB)
	{
		return new OBVector3(vecA.x()-vecB.x(), vecA.y()-vecB.y(),vecA.z()-vecB.z());
	}
	
	public static OBVector3 operator -(OBVector3 vecA, OBVector3 vecB)
	{
		return new OBVector3(vecA.x()-vecB.x(), vecA.y()-vecB.y(),vecA.z()-vecB.z());
	}
	
	public void Negate()
	{
		SetX(-x());
		SetY(-y());
		SetZ(-z());	
	}
	
	public static OBVector3 operator -(OBVector3 vec)
	{
		return new OBVector3(-vec.x(),-vec.y(),-vec.z());
	}

	public double this[int index]
	{
		get
		{
		if(index == 0)
			return x();
		else if(index == 1)
			return y();
		else if(index == 2)
			return z();
		else
			throw new IndexOutOfRangeException("Largest allowable index is 2");		
		}
	}
	
	public static OBVector3 VX
	{
		get{return new OBVector3(1,0,0);}
	}
	
	public static OBVector3 VY
	{
		get{return new OBVector3(0,1,0);}
	}
	
	public static OBVector3 VZ
	{
		get{return new OBVector3(0,0,1);}
	}
	
	public static OBVector3 VZero
	{
		get{return new OBVector3(0,0,0);}
	}
%}

%include "carrays.i"

%define WRAP_ARRAY(TYPE, NAME)
%array_class(TYPE,NAME)
%typemap(cstype) TYPE * "NAME"

%typemap(csout, excode=SWIGEXCODE) TYPE * 
{
	IntPtr cPtr = $imcall;$excode
	$csclassname ret = null;
	if (cPtr != IntPtr.Zero)
	{
        ret = new $csclassname(cPtr,true);
	}
	return NAME.frompointer(ret);
}
%typemap(csin, pre="     $csclassname temp$csinput = $csinput.cast();")
	TYPE * "$csclassname.getCPtr(temp$csinput)"
	
%typemap(csvarin, excode=SWIGEXCODE2) TYPE * %{
set
{
    $csclassname tempvalue = value.cast();
    $imcall;$excode
} %}

%typemap(csvarout, excode=SWIGEXCODE2) TYPE * %{
get
{
    IntPtr cPtr = $imcall;$excode;
    $csclassname tmp = new $csclassname(cPtr,true);$excode
    return NAME.frompointer(tmp);
} %}
%enddef

WRAP_ARRAY(double,double_array);


%include "std_pair.i"


%ignore OpenBabel::OBRingSearch;
%ignore OpenBabel::OBSSMatch;
%ignore OpenBabel::DoubleType;
%ignore OpenBabel::CharPtrLess;
%ignore OpenBabel::obLogBuf;
%ignore OpenBabel::OBMol::SetTitle(std::string &);
%ignore OpenBabel::OBFingerPrint::Tanimoto(const std::vector<unsigned int>&, const unsigned int*);
%ignore OpenBabel::OBRing::SetType(std::string &);
%ignore OpenBabel::OBAngle::SetAtoms(triple<OBAtom*,OBAtom*,OBAtom*> &);
%ignore OpenBabel::OBTorsion::AddTorsion(quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> &);
%ignore OpenBabel::OBGridData::GetAxes( double x[3], double y[3], double z[3] ) const;
%ignore OpenBabel::DoubleMultiply( unsigned int,unsigned int,DoubleType*);
%ignore OpenBabel::DoubleAdd( DoubleType*,unsigned int);
%ignore OpenBabel::DoubleModulus( DoubleType*,unsigned int);
%ignore OpenBabel::OBAtom::SetCoordPtr(double **);
%ignore OpenBabel::OBAtom::ClearCoordPtr();
%ignore OpenBabel::OBAtom::GetCIdx();;
%ignore OpenBabel::OBTypeTable::Translate(std::string &, const std::string &);
%ignore OpenBabel::vector3::Set(const double *c);
%ignore OpenBabel::vector3::Get();
%ignore *::operator->() const;

//ignored because of difficulty with wrapping type info
//type info could be marshalled to a C# struct if it is
//required later.
%ignore OpenBabel::OBFormat::GetType();

//ignored to remove failed mapping of maps and because C# classes
//aren't going to be used by the babel program 
%ignore OpenBabel::OBOp;
%ignore OpenBabel::OBBase::DoTransformations(const std::map<std::string, std::string>*);


//ignored because IO doesn't really belong in this class and streams don't work yet 
//*** check if this is ignored without the & operator in the sig
%ignore OpenBabel::OBSmartsPattern::WriteMapList(std::ostream);

//ignore some deprecated methods
%ignore OpenBabel::OBAtom::GetNextAtom();
%ignore OpenBabel::OBMol::CreateAtom();
%ignore OpenBabel::OBMol::CreateBond();
%ignore OpenBabel::OBMol::CreateResidue();
%ignore OpenBabel::OBMol::GetFirstAtom();
%ignore OpenBabel::OBMol::HasAromaticCorrected();
%ignore OpenBabel::ThrowError(char*);
%ignore OpenBabel::ThrowError(std::string &);
%ignore OpenBabel::OBConversion::GetNextFormat(Formatpos&, const char*&,OBFormat*&);
%ignore OpenBabel::OBAtom::GetCIdx();
%ignore OpenBabel::OBAtom::GetNextAtom();
%ignore OpenBabel::OBAtom::KBOSum();
%ignore OpenBabel::OBBond::Visit;
%ignore OpenBabel::OBBond::SetBO(int);
%ignore OpenBabel::OBBond::GetBO();
%ignore OpenBabel::OBBond::SetKSingle();
%ignore OpenBabel::OBBond::SetKDouble();
%ignore OpenBabel::OBBond::SetKTriple();
%ignore OpenBabel::OBBond::IsKSingle();
%ignore OpenBabel::OBBond::IsKDouble();
%ignore OpenBabel::OBBond::IsKTriple();
%ignore OpenBabel::OBResidueData::LookupBO(const std::string &);
%ignore OpenBabel::OBResidueData::AssignBonds(OBMol &,OBBitVec &);
%ignore OpenBabel::OBFloatGrid::GetMin(double *);
%ignore OpenBabel::OBFloatGrid::GetMax(double *);
%ignore OpenBabel::OBFloatGrid::GetSpacing(double &);
%ignore OpenBabel::OBFloatGrid::GetDim(int *);
%ignore OpenBabel::OBFloatGrid::SetLimits(const double[3], const double[3], const double[3],const double[3]);
%ignore OpenBabel::OBFloatGrid::GetVals();
%ignore OpenBabel::OBFloatGrid::SetVals(double *);
%ignore OpenBabel::OBRing::PathSize();
%ignore OpenBabel::OBTypeTable::Translate(char *, const char *);
%ignore OpenBabel::OBAtomTyper::CorrectAromaticNitrogens(OBMol &);
%ignore OpenBabel::OBConversion::GetNextFormat(Formatpos&, const char*&,OBFormat*&);
%ignore OpenBabel::OBBitVec::Empty();
%ignore OpenBabel::OBBitVec::BitIsOn(int);
%ignore OpenBabel::OBForceField::UpdateCoordinates(OBMol &);
%ignore OpenBabel::OBForceField::UpdateConformers(OBMol &);
%ignore OpenBabel::OBProxGrid;
%ignore OpenBabel::OBScoreGrid;
%ignore OpenBabel::score_t;

//moved these to static methods of the OBVector3 class which should
//be more intuitive for C# programmers
%ignore OpenBabel::VX;
%ignore OpenBabel::VY;
%ignore OpenBabel::VZ;
%ignore OpenBabel::VZero;

//ignoring the clone method of OBGenericData and its subclasses
//must implementing class just return null and if we really want to 
//clone we can do it on the managed side and implement a proper ICloneable/Clone
%ignore *::Clone(OBBase*) const;

//not needed in C#, programmers can use methods of the string class
//to do the same thing
%ignore OpenBabel::NewExtension(string &,char *);

//deprecated but not ignored
//OBElementTable::GetAtomicNum(const char *)
//OBResidueData::AssignBonds(OBMol &,OBBitVec &);
//
// these methods are not ignored because the alternative
// is an operator not yet wrapped.
//matrix3x3::Set(int row,int column, double v);
//matrix3x3::Get(int row,int column);

//ignore these global arrays, they are
//replaced by collections defined in
//the C# module class
%ignore OpenBabel::ElemDesc;
%ignore OpenBabel::Residue;

//Ignore these global constants to streamline the
//imclass and avoid using unmanaged function calls 
//to retrieve global variables with constant values. 
//Use moduleclass.ElemDesc.Count to get the same value. 
//If these are desired later add them to the module class
//as const values using pragmas.
%ignore OpenBabel::ElemNo;
%ignore OpenBabel::ResNo;

//these methods will be ignored until we figure out
//how to mod the wrappers in  std_map.i.
%ignore OpenBabel::OBPlugin::Begin(const char*);
%ignore OpenBabel::OBPlugin::End(const char*);

//***OBPlugin::GetMap() is not being ignored
%ignore OpenBabel::OBPlugin::GetMap() const;
%ignore OpenBabel::OBConversion::GetOptions(Option_type);

//ignore this constructor and method, multi-dimensional arrays are
//not well supported within .Net interop
%ignore OpenBabel::matrix3x3(double[3][3]);
%ignore OpenBabel::transform3d(double[3][3], double[3]);

//ignore these methods that return an iterator
//swig has issues wrapping the iterators and the
//anonymous IEnumerables added to some typemaps replace
//these methods where necessary.
%ignore OpenBabel::OBMol::BeginAtom(OBAtomIterator &);
%ignore OpenBabel::OBMol::NextAtom(OBAtomIterator &);
%ignore OpenBabel::OBMol::BeginBond(OBBondIterator &);
%ignore OpenBabel::OBMol::NextBond(OBBondIterator &);
%ignore OpenBabel::OBMol::BeginResidue(OBResidueIterator &);
%ignore OpenBabel::OBMol::NextResidue(OBResidueIterator &);
%ignore OpenBabel::OBMol::BeginInternalCoord(std::vector<OBInternalCoord*>::iterator &);
%ignore OpenBabel::OBMol::NextInternalCoord(std::vector<OBInternalCoord*>::iterator &);
%ignore OpenBabel::OBMol::BeginAtoms();
%ignore OpenBabel::OBMol::EndAtoms();
%ignore OpenBabel::OBMol::BeginBonds();
%ignore OpenBabel::OBMol::EndBonds();
%ignore OpenBabel::OBMol::BeginResidues();
%ignore OpenBabel::OBMol::EndResidues();

//***conformer methods not yet replaced
%ignore OpenBabel::OBMol::BeginConformer(std::vector<double*>::iterator&);
%ignore OpenBabel::OBMol::NextConformer(std::vector<double*>::iterator&);

%ignore OpenBabel::OBAtom::BeginAtoms();
%ignore OpenBabel::OBAtom::EndAtoms();
%ignore OpenBabel::OBAtom::BeginBonds();
%ignore OpenBabel::OBAtom::EndBonds();
%ignore OpenBabel::OBAtom::BeginBond(OBBondIterator &);
%ignore OpenBabel::OBAtom::NextBond(OBBondIterator &);
%ignore OpenBabel::OBAtom::BeginNbrAtom(OBAtomIterator &);
%ignore OpenBabel::OBAtom::NextNbrAtom(OBAtomIterator &);

//***should OBSmartsPattern be Enumerable?
%ignore OpenBabel::OBSmartsPattern::BeginMList();
%ignore OpenBabel::OBSmartsPattern::EndMList();

//***these are not yet replaced, should OBSetData be Enumberable?
%ignore OpenBabel::OBSetData::GetBegin();
%ignore OpenBabel::OBSetData::GetEnd();

//*** not yet replaced with enumerables
%ignore OpenBabel::OBRingData::BeginRings();
%ignore OpenBabel::OBRingData::EndRings();
%ignore OpenBabel::OBRingData::BeginRing(std::vector<OBRing*>::iterator &);
%ignore OpenBabel::OBRingData::NextRing(std::vector<OBRing*>::iterator &);



//these methods provide C# IEnumerables to allow easy iteration and
//enable the use of generic extension methods.
%typemap(csimports) OpenBabel::OBMol "
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
"
%typemap(cscode) OpenBabel::OBMol
%{
	public IEnumerable<OBAtom> Atoms()
	{
		OBMolAtomIter iter = new OBMolAtomIter(this);
		
		while(iter.MoveNext())
		{
			yield return iter.Current;
		}
	}
	
	public IEnumerable<OBBond> Bonds()
	{
		OBMolBondIter iter = new OBMolBondIter(this);
		
		while(iter.MoveNext())
		{
			yield return iter.Current;
		}
	}
	
	public IEnumerable<OBResidue> Residues()
	{
		OBResidueIter iter = new OBResidueIter(this);
		
		while(iter.MoveNext())
		{
			yield return iter.Current;
		}
	}
	
	public IEnumerable<OBRing> Rings()
	{
		return GetSSSR().Cast<OBRing>();
	}
	
%}

%typemap(csimports) OpenBabel::OBAtom "
using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
"
%typemap(cscode) OpenBabel::OBAtom
%{
	public IEnumerable<OBBond> Bonds()
	{
		OBAtomBondIter iter = new OBAtomBondIter(this);
		
		while(iter.MoveNext())
		{
			yield return iter.Current;
		}
	}
	
	public IEnumerable<OBAtom> Neighbors()
	{
		OBAtomAtomIter iter = new OBAtomAtomIter(this);
		
		while(iter.MoveNext())
		{
			yield return iter.Current;
		}
	}
%}

%module openbabelcsharp
%{
// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/obutil.h>
#include <openbabel/rand.h>
#include <openbabel/math/vector3.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/math/transform3d.h>
#include <openbabel/generic.h>

#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/residue.h>
#include <openbabel/internalcoord.h>

#include <openbabel/ring.h>
#include <openbabel/obconversion.h>
#include <openbabel/oberror.h>
#include <openbabel/plugin.h>
#include <openbabel/fingerprint.h>
#include <openbabel/descriptor.h>
#include <openbabel/format.h>

#include <openbabel/forcefield.h>
#include <openbabel/op.h>

#include <openbabel/bitvec.h>
#include <openbabel/data.h>
#include <openbabel/parsmart.h>
#include <openbabel/alias.h>
#include <openbabel/atomclass.h>

#include <openbabel/kinetics.h>
#include <openbabel/rotamer.h>

%}

%module VecMath
%{
#include <openbabel/math/vector3.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/math/transform3d.h>
%}



%ignore *::operator=;
%ignore *::operator-=;
%ignore *::operator+=;
%ignore *::operator*=;
%ignore *::operator/=;
%ignore *::operator <<;
%ignore *::operator==;
%ignore *::operator-;
//%ignore *::operator*;
%ignore *::operator !=;

%include "std_string.i"
%include "std_map.i"
%include "std_vector.i"

%template (VectorInt)             std::vector<int>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(VectorInt, std::vector<int>);
// Note that the following line will fail if the space between 
// the two greater-than signs is removed!
%template (VectorVecInt)                 std::vector<std::vector<int> >;

//**why does this not wrap correctly???
//%template (vvUInt)                std::vector<std::vector<unsigned int> >;
%template (VectorUShort)	      std::vector<unsigned short>;
%template (VectorUInt)     std::vector<unsigned int>;
%template (VectorDouble)          std::vector<double>;
%template (VectorString)          std::vector<std::string>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBVector3, OpenBabel::vector3);
%template (VectorOBVector3)              std::vector<OpenBabel::vector3>;

//**check this after build
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBExternalBond,OpenBabel::OBExternalBond);
%template (VectorOBExternalBond)	std::vector<OpenBabel::OBExternalBond>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBMol, OpenBabel::OBMol);
%template (VectorMol)     std::vector<OpenBabel::OBMol>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBBond, OpenBabel::OBBond);
%template (VectorBond)    std::vector<OpenBabel::OBBond>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBResidue, OpenBabel::OBResidue);
%template (VectorResidue) std::vector<OpenBabel::OBResidue>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBRing, OpenBabel::OBRing);
%template (VectorRing)    std::vector<OpenBabel::OBRing>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBTorsion, OpenBabel::OBTorsion);
%template (VectorTorsion)     std::vector<OpenBabel::OBTorsion>;

// Note that vectors of pointers need slightly different syntax
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBRing, OpenBabel::OBRing*);
%template (VectorpRing)   std::vector<OpenBabel::OBRing*>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBGenericData, OpenBabel::OBGenericData*);
%template (VectorpData)    std::vector<OpenBabel::OBGenericData*>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBInternalCoord, OpenBabel::OBInternalCoord*);
%template (VectorpInternalCoord)		std::vector<OpenBabel::OBInternalCoord*>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBAtom, OpenBabel::OBAtom*);
%template (VectorpAtom)		std::vector<OpenBabel::OBAtom*>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBBond, OpenBabel::OBBond*);
%template (VectorpBond)		std::vector<OpenBabel::OBBond*>;


%import <openbabel/babelconfig.h>

%include <openbabel/data.h>
%include <openbabel/rand.h>
%include <openbabel/obutil.h>
%include <openbabel/math/vector3.h>
%include <openbabel/math/matrix3x3.h>
%include <openbabel/math/transform3d.h>
%include <openbabel/math/spacegroup.h>
%include <openbabel/base.h>

%include <openbabel/generic.h>
%define CAST_GENERICDATA_TO(subclass)
%extend OpenBabel::OBGenericData {
    OpenBabel::OB ## subclass *To ## subclass() {
	    return (OpenBabel::OB ## subclass *) $self;
    }
};
%enddef
CAST_GENERICDATA_TO(AngleData)
CAST_GENERICDATA_TO(AtomClassData)
CAST_GENERICDATA_TO(ChiralData)
CAST_GENERICDATA_TO(CommentData)
CAST_GENERICDATA_TO(ConformerData)
CAST_GENERICDATA_TO(ExternalBondData)
CAST_GENERICDATA_TO(GridData)
CAST_GENERICDATA_TO(MatrixData)
CAST_GENERICDATA_TO(NasaThermoData)
CAST_GENERICDATA_TO(PairData)
// CAST_GENERICDATA_TO(PairTemplate)
CAST_GENERICDATA_TO(RateData)
CAST_GENERICDATA_TO(RotamerList)
CAST_GENERICDATA_TO(RotationData)
CAST_GENERICDATA_TO(SerialNums)
CAST_GENERICDATA_TO(SetData)
CAST_GENERICDATA_TO(SymmetryData)
CAST_GENERICDATA_TO(TorsionData)
CAST_GENERICDATA_TO(UnitCell)
CAST_GENERICDATA_TO(VectorData)
CAST_GENERICDATA_TO(VibrationData)
CAST_GENERICDATA_TO(VirtualBond)

%include <openbabel/griddata.h> // Needs to come after generic.h

%include <openbabel/chains.h>
//# %import <openbabel/bitvec.h>
%include <openbabel/typer.h>

%include <openbabel/plugin.h>

%include <openbabel/oberror.h>
%include <openbabel/format.h>
%include <openbabel/obconversion.h>
%include <openbabel/residue.h>
%include <openbabel/internalcoord.h>
%include <openbabel/atom.h>
%include <openbabel/bond.h>
%ignore OpenBabel::OBMol::SetData;
%include <openbabel/mol.h>
%include <openbabel/ring.h>
%include <openbabel/parsmart.h>
%include <openbabel/alias.h>
%include <openbabel/atomclass.h>

%include <openbabel/fingerprint.h>
%include <openbabel/descriptor.h>
%include <openbabel/forcefield.h>

%include <openbabel/op.h>

%include <openbabel/bitvec.h>

# The following %ignores avoid warning messages due to shadowed classes.
# This does not imply a loss of functionality as (in this case)
# the shadowed class is identical (from the point of view of SWIG) to
# the shadowing class.
# This is because C++ references (&) are transformed by SWIG back into
# pointers, so that OBAtomIter(OBMol &) would be treated the same as
# OBAtomIter(OBMol *).

%ignore OBAtomAtomIter(OBAtom &);
%ignore OBAtomBondIter(OBAtom &);
%ignore OBMolAngleIter(OBMol &);
%ignore OBMolAtomIter(OBMol &);
%ignore OBMolAtomBFSIter(OBMol &);
%ignore OBMolAtomDFSIter(OBMol &);
%ignore OBMolAtomBFSIter(OBMol &, int);
%ignore OBMolAtomDFSIter(OBMol &, int);
%ignore OBMolBondIter(OBMol &);
%ignore OBMolPairIter(OBMol &);
%ignore OBMolRingIter(OBMol &);
%ignore OBMolTorsionIter(OBMol &);
%ignore OBResidueIter(OBMol &);
%ignore OBResidueAtomIter(OBResidue &);

//macro for wrapping iterators
%define WRAPITERATOR(NAME,CTYPE,RETYPE)
%ignore CTYPE::NAME();
%csmethodmodifiers CTYPE::operator* "protected";
%csmethodmodifiers CTYPE::operator++ "protected";
%rename(obAdvance) CTYPE::operator++;
%csmethodmodifiers CTYPE::operator bool() const "protected";
%rename(obHasNext) CTYPE::operator bool;
//***redundant?
%typemap(cstype) CTYPE "$csclassname"
%typemap(csinterfaces) CTYPE "IEnumerator<RETYPE>"
%typemap(csimports) CTYPE "
using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
"
%typemap(cscode) CTYPE
%{
	protected bool iterating;

	public virtual RETYPE Current
	{
		get{return iterating ? __ref__() : null;}
	}
	
	object System.Collections.IEnumerator.Current
	{
		get{return Current;}
	}
	
	public virtual bool MoveNext()
	{
         if (obHasNext())
         {
             if (!iterating)
             {
                 iterating = true;
                 return true;
             }
             obAdvance();
             bool ret = obHasNext();

             if (!ret)
                 iterating = false;
                
             return ret;
         }
            
         return false;
	}
	
	public void Reset()
	{
		throw new InvalidOperationException("Reset is not a supported operation");
	}
%} 
%enddef


WRAPITERATOR(OBMolAtomIter,OpenBabel::OBMolAtomIter,OBAtom);
WRAPITERATOR(OBMolAtomDFSIter,OpenBabel::OBMolAtomDFSIter,OBAtom);
WRAPITERATOR(OBMolAtomBFSIter,OpenBabel::OBMolAtomBFSIter,OBAtom);
WRAPITERATOR(OBMolBondIter,OpenBabel::OBMolBondIter,OBBond);
WRAPITERATOR(OBMolAngleIter,OpenBabel::OBMolAngleIter,VectorUInt);
WRAPITERATOR(OBAtomAtomIter,OpenBabel::OBAtomAtomIter,OBAtom)
WRAPITERATOR(OBAtomBondIter,OpenBabel::OBAtomBondIter,OBBond);
WRAPITERATOR(OBMolRingIter,OpenBabel::OBRingIter,OBRing);
WRAPITERATOR(OBMolTorsionIter,OpenBabel::OBTorsionIter,VectorUInt);
WRAPITERATOR(OBResidueIter,OpenBabel::OBResidueIter,OBResidue);
WRAPITERATOR(OBResidueAtomIter,OpenBabel::OBResidueAtomIter,OBAtom);
WRAPITERATOR(OBMolPairIter,OpenBabel::OBMolPairIter,VectorUInt)

//***add CurrentDepth to BFS iterator


%include <openbabel/obiter.h>
