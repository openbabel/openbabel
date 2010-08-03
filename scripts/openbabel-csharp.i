
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

%rename(AtomRefType) atomreftype;
%rename(ErrorQualifier) errorQualifier;
%rename(OBMessageLevel) obMessageLevel;
%rename(Rate_Type) rate_type;
%rename(Reaction_Type) reaction_type;

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
//when using OBDotNet to build their own apis
//this change can be made to other classes if needed in the future.
%typemap(csclassmodifiers) OpenBabel::OBAtom "public partial class"
%typemap(csclassmodifiers) OpenBabel::OBMol "public partial class"
%typemap(csclassmodifiers) OpenBabel::OBBond "public partial class"
%typemap(csclassmodifiers) OpenBabel::OBRing "public partial class"



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
  
//simplified public Downcast method
//this is defined up here because something
//lower down in the file interferes with it
//
%typemap(cscode) OpenBabel::OBGenericData
%{

  public virtual DType Downcast<DType>() where DType : OBGenericData
  {
      string derivedType = typeof(DType).Name;
      derivedType = derivedType.StartsWith("OB") ? derivedType.Substring(2) :derivedType;
      string castMethodName = "";
      
      castMethodName = string.Format("{0}{1}","To",derivedType);
      
      System.Reflection.MethodInfo castMethod = typeof(OBGenericData).GetMethod(castMethodName);
      
      if(castMethod == null)
        throw new InvalidCastException("No explicit downcast is defined for " + derivedType);
      
      return (DType)castMethod.Invoke(this,null);
   }
%}

//disable Downcast method in derived types
%define DISABLE_DOWNCAST(CSCLASS)
%typemap(cscode) OpenBabel::CSCLASS
%{
  public override DType Downcast<DType>()
  {
      throw new NotImplementedException("Downcast<DTYPE> is not implemented in " + "CSCLASS");
   }
%}
%enddef
DISABLE_DOWNCAST(AliasData);
DISABLE_DOWNCAST(OBAngleData);
DISABLE_DOWNCAST(OBAtomClassData);
DISABLE_DOWNCAST(OBChiralData);
DISABLE_DOWNCAST(OBCommentData);
DISABLE_DOWNCAST(OBConformerData);
DISABLE_DOWNCAST(OBExternalBondData);
DISABLE_DOWNCAST(OBGridData);
DISABLE_DOWNCAST(OBMatrixData);
DISABLE_DOWNCAST(OBNasaThermoData);
DISABLE_DOWNCAST(OBPairData);
DISABLE_DOWNCAST(OBPairTemplate);
DISABLE_DOWNCAST(OBRateData);
DISABLE_DOWNCAST(OBRotamerList);
DISABLE_DOWNCAST(OBRotationData);
DISABLE_DOWNCAST(OBSerialNums);
DISABLE_DOWNCAST(OBSetData);
DISABLE_DOWNCAST(OBSymmetryData);
DISABLE_DOWNCAST(OBTorsionData);
DISABLE_DOWNCAST(OBUnitCell);
DISABLE_DOWNCAST(OBVectorData);
DISABLE_DOWNCAST(OBVibrationData);
DISABLE_DOWNCAST(OBVirtualBond);

//extra method to provide access to a rotor
//iterator
%ignore OpenBabel::OBRotorList::BeginRotors();
%ignore OpenBabel::OBRotorList::EndRotors();
%ignore OpenBabel::OBRotorList::BeginRotor(OBRotorIterator &);
%ignore OpenBabel::OBRotorList::NextRotor(OBRotorIterator &);
%extend OpenBabel::OBRotorList
{
	std::vector<OpenBabel::OBRotor*> GetRotors()
	{
		std::vector<OpenBabel::OBRotor*> rotors($self->BeginRotors(),$self->EndRotors());
		return rotors;
	}
};


%include "carrays.i"
%define WRAP_ARRAY(TYPE, NAME)
//can't seal array classes because of the dispose method
//%typemap(csclassmodifiers) NAME "public sealed class"
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
WRAP_ARRAY(double,CDoubleArray);
WRAP_ARRAY(int, CIntArray);
WRAP_ARRAY(unsigned char, CByteArray);
WRAP_ARRAY(unsigned int, CUIntArray);

//**Trying to set up array marshalling in pinvoke class.
//Finishing this is a priority for v0.3.
//
//%include "arrays_csharp.i"
//CSHARP_ARRAYS(double, double)
//OpenBabel::OBMol::AddConformer(double* f);
//%apply double INPUT[] { double* c }
//OpenBabel::OBMol::SetCoordinates(double* c);
//%apply double INPUT[] { double* arg0,double* arg1 }
//OpenBabel::get_rmat(double* arg0, double* arg1, int arg3);

%ignore OpenBabel::OBRingSearch;
%ignore OpenBabel::OBSSMatch;
%ignore OpenBabel::DoubleType;
%ignore OpenBabel::CharPtrLess;
%ignore OpenBabel::obLogBuf;
%ignore OpenBabel::OBMol::SetTitle(std::string &);
//I can't get this overload to be ignored.
%ignore OpenBabel::OBFingerPrint::Tanimoto(const std::vector<unsigned int>&, const unsigned int* );
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
 
//remove ignore to access OBOp for using Gen3D but still ignore unsupported
//methods
%ignore OpenBabel::OBOp::DoOps(OBBase*, OpMap*);
%ignore OpenBabel::OBOp::Do(OBBase*, OpMap*, const char*);
%ignore OpenBabel::OBOp::Do(OBBase*,OpMap*);
%ignore OpenBabel::OBBase::DoTransformations(const std::map<std::string, std::string>*);


//ignored because streams don't work yet 
//*** check if this is ignored without the & operator in the sig
%ignore OpenBabel::OBSmartsPattern::WriteMapList(std::ostream);
//***not being ignored
%ignore OpenBabel::OBConversion::GetInPos() const;
%ignore OpenBabel::FastSearchIndexer::Add(OBBase*, std::streampos);

//ignore until std::pair<int,int> is wrapped successfuly
%ignore OpenBabel::OBAromaticTyper;
%ignore OpenBabel::aromtyper;


//No longer ignoring the deprecated methods.
//I've decided there is no real point to doing so.

//ignore until std::pair is wrapped correctly
//or add overloads using the extend directive that take
//std::vector<std::vector<int>> 
%ignore OpenBabel::OBSmartsPattern::RestrictedMatch(OBMol &, std::vector<std::pair<int,int> > &, bool);
%ignore OpenBabel::OBSmartsPattern::RestrictedMatch(OBMol &, std::vector<std::pair<int,int> > &);
%ignore OpenBabel::SmartsLexReplace(std::string &,std::vector<std::pair<std::string,std::string> > &);

//%ignore OpenBabel::OBRotor::SetDihedralAtoms(int[4]);


//ignoring the clone method of OBGenericData and its subclasses
//most implementing class just return null and if we really want to 
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

//wrappers for ElemDesc and Residue global arrays with readonly
//collections this is a temporary work around until I write an immutable
//list type that can be declared with an object initializer
%pragma(csharp) moduleimports=
%{
	using System;
	using System.Collections.Generic;
	using System.Collections.ObjectModel;
	using System.Runtime.InteropServices;
%}

%pragma(csharp) modulecode=
%{
	private static readonly IList<string> elemDesc = new List<string>()
	{
		" N  "," CA "," C  "," O  ",@" C\ "," OT ",
		" S  "," P  "," O1P"," O2P", " O5*"," C5*"," C4*",
		" O4*"," C3*"," O3*"," C2*"," O2*"," C1*"," CA2",
		" SG "," N1 "," N2 "," N3 "," N4 "," N6 "," O2 "," O4 "," O6 ",
	};
	
	public static ReadOnlyCollection<string> ElemDesc
	{
		get{return new ReadOnlyCollection<string>(elemDesc);}
	} 

	private static readonly IList<string> residue = new List<string>()
	{
    "ALA", /* 8.4% */     "GLY", /* 8.3% */
    "LEU", /* 8.0% */     "SER", /* 7.5% */
    "VAL", /* 7.1% */     "THR", /* 6.4% */
    "LYS", /* 5.8% */     "ASP", /* 5.5% */
    "ILE", /* 5.2% */     "ASN", /* 4.9% */
    "GLU", /* 4.9% */     "PRO", /* 4.4% */
    "ARG", /* 3.8% */     "PHE", /* 3.7% */
    "GLN", /* 3.5% */     "TYR", /* 3.5% */
    "HIS", /* 2.3% */     "CYS", /* 2.0% */
    "MET", /* 1.8% */     "TRP", /* 1.4% */
    "ASX", "GLX", "PCA", "HYP",
    "  A", "  C", "  G", "  T",
    "  U", " +U", "  I", "1MA",
    "5MC", "OMC", "1MG", "2MG",
    "M2G", "7MG", "OMG", " YG",
    "H2U", "5MU", "PSU",
    "UNK", "ACE", "FOR", "HOH",
    "DOD", "SO4", "PO4", "NAD",
    "COA", "NAP", "NDP"
	};
	
	public static ReadOnlyCollection<string> Residue
	{
		get{return new ReadOnlyCollection<string>(residue);}
	}	
%}

//make the module class static
%pragma(csharp) moduleclassmodifiers="public static class"

//Ignore these global constants to streamline the
//imclass and avoid using unmanaged function calls 
//to retrieve global variables with constant values. 
//Use moduleclass.ElemDesc.Count to get the same value. 
//If these are desired later add them to the module class
//as const values using pragmas
%ignore OpenBabel::ElemNo;
%ignore OpenBabel::ResNo;

//these methods will be ignored until we figure out
//how to mod the wrappers in  std_map.i.
%ignore OpenBabel::OBPlugin::Begin(const char*);
%ignore OpenBabel::OBPlugin::End(const char*);

//***OBPlugin::GetMap() is not being ignored
%ignore OpenBabel::OBPlugin::GetMap() const;
%ignore OpenBabel::OBConversion::GetOptions(Option_type);
%ignore OpenBabel::OBConversion::GetNextFormat(Formatpos&, const char*&,OBFormat*&);

//ignore this constructor and method, multi-dimensional arrays are
//not well supported within .Net interop
%ignore matrix3x3(double[3][3]);
%ignore transform3d(double[3][3], double[3]);
%ignore OpenBabel::OBMol::Rotate(const double[3][3]);
%ignore OpenBabel::ob_make_rmat(double[3][3],double[9]);
%ignore OpenBabel::qtrfit (double *,double *,int size,double[3][3]);
%ignore OpenBabel::OBAngleData::FillAngleArray(int**, unsigned int&);
%ignore OpenBabel::SpaceGroup::Transform(const vector3 &) const;
%ignore OpenBabel::SpaceGroup::NextTransform(transform3dIterator &i) const;
%ignore OpenBabel::SpaceGroup::BeginTransform(transform3dIterator &i) const;

//ignore these methods using iterators
//swig has issues wrapping the iterators and the
//anonymous IEnumerables added to some typemaps replace
//these methods where necessary.
%ignore OpenBabel::OBBase::BeginData();
%ignore OpenBabel::OBBase::EndData();
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
%ignore OpenBabel::OBAtom::BeginAtoms();
%ignore OpenBabel::OBAtom::EndAtoms();
%ignore OpenBabel::OBAtom::BeginBonds();
%ignore OpenBabel::OBAtom::EndBonds();
%ignore OpenBabel::OBAtom::BeginBond(OBBondIterator &);
%ignore OpenBabel::OBAtom::NextBond(OBBondIterator &);
%ignore OpenBabel::OBAtom::BeginNbrAtom(OBAtomIterator &);
%ignore OpenBabel::OBAtom::NextNbrAtom(OBAtomIterator &);
%ignore OpenBabel::OBAtom::BeginNbrAtom(OBBondIterator &);
%ignore OpenBabel::OBAtom::NextNbrAtom(OBBondIterator &);
%ignore OpenBabel::OBAtom::InsertBond(OBBondIterator &,OBBond*);
%ignore OpenBabel::OBResidue::BeginAtoms();
%ignore OpenBabel::OBResidue::EndAtoms();
%ignore OpenBabel::OBResidue::BeginAtom(std::vector<OBAtom*>::iterator &);
%ignore OpenBabel::OBResidue::NextAtom(std::vector<OBAtom*>::iterator &);

//use C#  mol.GetConformers()added below
%ignore OpenBabel::OBMol::BeginConformer(std::vector<double*>::iterator&);
%ignore OpenBabel::OBMol::NextConformer(std::vector<double*>::iterator&);
%ignore OpenBabel::OBMol::GetConformers();

//Use OBSmartsPattern.GetMapList().GetEnumerator() 
//to obtain an equivalent iterator.
%ignore OpenBabel::OBSmartsPattern::BeginMList();
%ignore OpenBabel::OBSmartsPattern::EndMList();

//Use OBSetData.GetData().GetEnumerator()
//to obtain an equivalent iterator.
%ignore OpenBabel::OBSetData::GetBegin();
%ignore OpenBabel::OBSetData::GetEnd();

//Use OBRingData.GetAllData().GetEnumerator() 
//to obtain an equivalent iterator.
%ignore OpenBabel::OBRingData::BeginRings();
%ignore OpenBabel::OBRingData::EndRings();
%ignore OpenBabel::OBRingData::BeginRing(std::vector<OBRing*>::iterator &);
%ignore OpenBabel::OBRingData::NextRing(std::vector<OBRing*>::iterator &);
%ignore OpenBabel::Swab(int);

//these methods provide C# IEnumerables to allow easy iteration and
//enable LINQ queries.
%typemap(csimports) OpenBabel::OBMol "
using System;
using System.Collections.Generic;
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
	
	public IEnumerable<OBMol> Fragments()
	{
		OBMolAtomDFSIter iter = new OBMolAtomDFSIter(this);
		OBMol ret = new OBMol();
		while(GetNextFragment(iter,ret))
		{
			yield return ret;
		}
	}
	
	//Temporary workaround until the std::vector<double*>
	//proxy is functional. The proxy class will implement
	//IEnumerable<CDoubleArray> any code using this method
	//will remain compatible with future versions.
	public IEnumerable<CDoubleArray> GetConformers()
	{
		int n = NumConformers();
		
		if(n == 0)
			yield return null;
		
		for(int i = 0; i < n; i++)
			yield return GetConformer(i);
	}
	
	
	//note - this method could potentially trigger 
	//GC issues by deleting CDoubleArray instances referenced
	//elsewhere in the code. 
	public void SetConformers(IEnumerable<CDoubleArray> confs)
	{	
		int n = NumConformers();
		
		if(n != 0)
			for(int i=0;i<n;i++)
				DeleteConformer(i);
				
		foreach(CDoubleArray cda in confs)
			AddConformer(cda);
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

//Changed module class name to comply with
//C# naming conventions
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
#include <openbabel/builder.h>
#include <openbabel/op.h>

#include <openbabel/bitvec.h>
#include <openbabel/data.h>
#include <openbabel/parsmart.h>
#include <openbabel/alias.h>
#include <openbabel/atomclass.h>

#include <openbabel/kinetics.h>
//OBReaction can't be mapped properly
//until shared_ptr is mapped
//#include <openbabel/reaction.h>
#include <openbabel/rotor.h>
#include <openbabel/rotamer.h>

%}

//ignore some currently unsupported operators
%ignore *::operator=;
%ignore *::operator-=;
%ignore *::operator+=;
%ignore *::operator*=;
%ignore *::operator/=;
%ignore *::operator <<;
%ignore *::operator==;
%ignore *::operator-;
%ignore *::operator*;
%ignore *::operator !=;
%ignore *::operator ++;
%ignore *::operator bool;
%ignore *::operator [];
%ignore *::operator /;
%ignore *::operator +;

%include "std_string.i"
%include "std_vector.i"

%typemap(cscode) std::vector<T>
%{
//ForEach method to simplify working around the foreach/delegate interaction "bug"
#if !SWIG_DOTNET_1
  public void ForEach(Action<CSTYPE> action)
  {
    if(action == null)
		throw new ArgumentNullException("action");
	
	for(int i = 0 ;i<Count; i++)
		action(getitem(i));
  }
#endif
%}

%template (VectorInt)             std::vector<int>;
// Note that the following line will fail if the space between 
// the two greater-than signs is removed!
%template (VectorVecInt)     std::vector<std::vector<int> >;
%template (VectorUInt)            std::vector<unsigned int>;
%template (VectorVecUInt)         std::vector<std::vector<unsigned int> >;
%template (VectorUShort)	      std::vector<unsigned short>;
%template (VectorDouble)          std::vector<double>;
%template (VectorVecDouble)		  std::vector<std::vector<double> >;
%template (VectorString)          std::vector<std::string>;
%template (VectorOBVector3)              std::vector<OpenBabel::vector3>;
%template (VectorVecOBVector3)    std::vector<std::vector<OpenBabel::vector3> >;

%template (VectorOBExternalBond)	std::vector<OpenBabel::OBExternalBond>;
%template (VectorMol)     std::vector<OpenBabel::OBMol>;
%template (VectorBond)    std::vector<OpenBabel::OBBond>;
%template (VectorResidue) std::vector<OpenBabel::OBResidue>;
%template (VectorRing)    std::vector<OpenBabel::OBRing>;
%template (VectorTorsion)     std::vector<OpenBabel::OBTorsion>;


// Note that vectors of pointers need slightly different syntax
%template (VectorpRing)   std::vector<OpenBabel::OBRing*>;
%template (VectorpData)    std::vector<OpenBabel::OBGenericData*>;
%template (VectorpInternalCoord)		std::vector<OpenBabel::OBInternalCoord*>;
%template (VectorpAtom)		std::vector<OpenBabel::OBAtom*>;
%template (VectorpBond)		std::vector<OpenBabel::OBBond*>;
%template (VectorpRotor)		std::vector<OpenBabel::OBRotor*>;

//the typemap for wrapping std::vector<double*> is going to need
//some customization to work with the CDoubleArray type.
//%template (VectorpDouble)		  std::vector<double*>;

%import <openbabel/babelconfig.h>

%warnfilter(516) OpenBabel::OBElementTable; // Ignoring std::string methods in favour of char* ones
%include <openbabel/data.h>
%include <openbabel/rand.h>
%include <openbabel/obutil.h>
%warnfilter(516) OpenBabel::vector3; // Using the const x(), y() and z() in favour of the non-const
%include <openbabel/math/vector3.h>
%warnfilter(503) OpenBabel::matrix3x3; // Not wrapping any of the overloaded operators
%include <openbabel/math/matrix3x3.h>
%include <openbabel/math/transform3d.h>
%warnfilter(516) OpenBabel::SpaceGroup; // Ignoring std::string methods in favour of char* ones
%include <openbabel/math/spacegroup.h>

# CloneData should be used instead of the following method
%ignore OpenBabel::OBBase::SetData;
%warnfilter(516) OpenBabel::OBBase; // Ignoring std::string methods in favour of char* ones
%include <openbabel/base.h>

%warnfilter(516) OpenBabel::OBPairData; // Ignoring std::string methods in favour of char* ones
%warnfilter(516) OpenBabel::OBSetData;
%warnfilter(516) OpenBabel::OBCommentData;

//replacement for method return unsupported std::pair
%ignore OpenBabel::OBTorsion::GetBC;
%extend OpenBabel::OBTorsion
{
	std::vector<OpenBabel::OBAtom*> GetBC()
	{
		std::vector<OpenBabel::OBAtom*> bcAtoms(2);
        
		bcAtoms[0] = $self->GetBC().first;
		bcAtoms[1] = $self->GetBC().second;
		
		return bcAtoms;
	}
};
%include <openbabel/generic.h>



//extending typemap to work around
//lack of support for void*
%ignore OpenBabel::OBRotor::GetRotAtoms;
%extend OpenBabel::OBRotor
{
	int* GetRotAtoms()
	{
		return (int*)$self->GetRotAtoms();	
	}
	
};

//FastSearch won't be fully functional until we figure out
//how to handle streams.
//We'd also have to do our own implementation of the multimap wrapper
%ignore OpenBabel::FastSearch::FindSimilar(OBBase*, multimap<double, unsigned int>&,int);
%ignore OpenBabel::FastSearch::FindSimilar(OBBase*, multimap<double, unsigned int>&,double);

//individual cast methods are private as the child classes
//do not need to inherit them. A C# Downcast<T> method
//is defined above.
//%csmethodifier must proceed %extend declaration
%define CAST_GENERICDATA_TO(subclass)
%csmethodmodifiers OpenBabel::OBGenericData::To ## subclass() "private";
%extend OpenBabel::OBGenericData {
    OpenBabel::OB ## subclass *To ## subclass() {
	    return (OpenBabel::OB ## subclass *) $self;
    }
};
%enddef

%typemap(cstype) OpenBabel::OBRotamerList* "OBRotamerList";

//why is AliasData not supported?
CAST_GENERICDATA_TO(AngleData);
CAST_GENERICDATA_TO(AtomClassData);
CAST_GENERICDATA_TO(ChiralData);
CAST_GENERICDATA_TO(CommentData);
CAST_GENERICDATA_TO(ConformerData);
CAST_GENERICDATA_TO(ExternalBondData);
CAST_GENERICDATA_TO(GridData);
CAST_GENERICDATA_TO(MatrixData);
CAST_GENERICDATA_TO(NasaThermoData);
CAST_GENERICDATA_TO(PairData);
// CAST_GENERICDATA_TO(PairTemplate);
CAST_GENERICDATA_TO(RateData);
CAST_GENERICDATA_TO(RotamerList);
CAST_GENERICDATA_TO(RotationData);
//OBSerialNum class will not be functional until std::map wrappers
//are functional.
//CAST_GENERICDATA_TO(SerialNums);
CAST_GENERICDATA_TO(SetData);
CAST_GENERICDATA_TO(SymmetryData);
CAST_GENERICDATA_TO(TorsionData);
CAST_GENERICDATA_TO(UnitCell);
CAST_GENERICDATA_TO(VectorData);
CAST_GENERICDATA_TO(VibrationData);
CAST_GENERICDATA_TO(VirtualBond);


%include <openbabel/griddata.h> // Needs to come after generic.h

%include <openbabel/chains.h>
%include <openbabel/typer.h>

// To avoid warning in plugin.h about "Nothing known about std::binary_function"
namespace std { 
        template <T1, T2, T3>
        class binary_function {}; 
}
%template(dummy) std::binary_function <const char *, const char *, bool>;
%include <openbabel/plugin.h>

// To avoid warning in oberror.h about "Nothing known about std::stringbuf"
namespace std { class stringbuf {}; }
%warnfilter(503) OpenBabel::OBError; // Not wrapping any of the overloaded operators
%include <openbabel/oberror.h>
%include <openbabel/format.h>
%include <openbabel/obconversion.h>
%include <openbabel/residue.h>
%include <openbabel/internalcoord.h>
%warnfilter(516) OpenBabel::OBAtom; // Using non-const version of GetVector
%include <openbabel/atom.h>
%warnfilter(516) OpenBabel::OBBond; // Using non-const versions of GetBeginAtom, GetEndAtom
%include <openbabel/bond.h>
%ignore OpenBabel::OBMol::SetData;
%include <openbabel/mol.h>
%include <openbabel/ring.h>
%warnfilter(516) OpenBabel::OBSmartsPattern; // Using non-const versions of GetSMARTS
%include <openbabel/parsmart.h>
%warnfilter(516) OpenBabel::AliasData; // Ignoring std::string methods in favour of char* ones
%include <openbabel/alias.h>
%include <openbabel/atomclass.h>
%ignore OpenBabel::FptIndex;
%include <openbabel/fingerprint.h>
%include <openbabel/descriptor.h>

# Ignore shadowed methods
%ignore OpenBabel::OBForceField::VectorSubtract(const double *const, const double *const, double *);
%ignore OpenBabel::OBForceField::VectorMultiply(const double *const, const double, double *);
%warnfilter(516) OpenBabel::OBForceField; // Ignoring std::string methods in favour of char* ones
%include <openbabel/forcefield.h>

%include <openbabel/builder.h>
%include <openbabel/op.h>

%warnfilter(503) OpenBabel::OBBitVec; // Not wrapping any of the overloaded operators
%include <openbabel/bitvec.h>

%include <openbabel/rotor.h>
%ignore OpenBabel::Swab;
%include <openbabel/rotamer.h>

//wrapping boost::shared_ptr
//this does not work yet
//%include <boost/shared_ptr.hpp> 
//%include "boost_shared_ptr.i"
//%include "shared_ptr.i"
//SWIG_SHARED_PTR_TYPEMAPS(OBMolSharedPtr, ,OpenBabel::OBMol)
//%include <openbabel/reaction.h>

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
%ignore OBMolBondBFSIter(OBMol &);
%ignore OBMolBondBFSIter(OBMol &, int);
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
WRAPITERATOR(OBMolBondBFSIter,OpenBabel::OBMolBondBFSIter,OBBond);
WRAPITERATOR(OBMolAngleIter,OpenBabel::OBMolAngleIter,VectorUInt);
WRAPITERATOR(OBAtomAtomIter,OpenBabel::OBAtomAtomIter,OBAtom)
WRAPITERATOR(OBAtomBondIter,OpenBabel::OBAtomBondIter,OBBond);
WRAPITERATOR(OBMolRingIter,OpenBabel::OBRingIter,OBRing);
WRAPITERATOR(OBMolTorsionIter,OpenBabel::OBTorsionIter,VectorUInt);
WRAPITERATOR(OBResidueIter,OpenBabel::OBResidueIter,OBResidue);
WRAPITERATOR(OBResidueAtomIter,OpenBabel::OBResidueAtomIter,OBAtom);
WRAPITERATOR(OBMolPairIter,OpenBabel::OBMolPairIter,VectorUInt)

%include <openbabel/obiter.h>
