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
    IntPtr cPtr = $imcall;
    $csclassname tmp = new $csclassname(cPtr,true);$excode
    return NAME.frompointer(tmp);
} %}
%enddef

WRAP_ARRAY(double,double_array)

%module openbabelcsharp
%{
// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/obutil.h>
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

#include <openbabel/kinetics.h>
#include <openbabel/rotamer.h>

#include <openbabel/chains.h>
#include <openbabel/obiter.h>
%}

%module VecMath
%{
#include <openbabel/math/vector3.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/math/transform3d.h>
%}

%ignore *::operator=;
%ignore *::operator++;
%ignore *::operator-=;
%ignore *::operator+=;
%ignore *::operator bool;
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

%template (vectorInt)             std::vector<int>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(vectorInt, std::vector<int>);
// Note that the following line will fail if the space between 
// the two greater-than signs is removed!
%template (vvInt)                 std::vector<std::vector<int> >;
%template (vectorUnsignedInt)     std::vector<unsigned int>;
%template (vectorDouble)          std::vector<double>;
%template (vectorString)          std::vector<std::string>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBVector3, OpenBabel::vector3);
%template (vectorOBVector3)              std::vector<OpenBabel::vector3>;

SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBMol, OpenBabel::OBMol);
%template (vectorMol)     std::vector<OpenBabel::OBMol>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBBond, OpenBabel::OBBond);
%template (vectorBond)    std::vector<OpenBabel::OBBond>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBResidue, OpenBabel::OBResidue);
%template (vectorResidue) std::vector<OpenBabel::OBResidue>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBRing, OpenBabel::OBRing);
%template (vectorRing)    std::vector<OpenBabel::OBRing>;

// Note that vectors of pointers need slightly different syntax
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBRing, OpenBabel::OBRing*);
%template (vectorpRing)   std::vector<OpenBabel::OBRing*>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBGenericData, OpenBabel::OBGenericData*);
%template (vectorpData)    std::vector<OpenBabel::OBGenericData*>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBInternalCoord, OpenBabel::OBInternalCoord*);
%template (vectorpInternalCoord)		std::vector<OpenBabel::OBInternalCoord*>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBAtom, OpenBabel::OBAtom*);
%template (vectorpAtom)		std::vector<OpenBabel::OBAtom*>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBBond, OpenBabel::OBBond*);
%template (vectorpBond)		std::vector<OpenBabel::OBBond*>;


%import <openbabel/babelconfig.h>

%include <openbabel/data.h>
%include <openbabel/obutil.h>
%include <openbabel/math/vector3.h>
%include <openbabel/math/matrix3x3.h>
%include <openbabel/math/transform3d.h>
%include <openbabel/math/spacegroup.h>
%include <openbabel/bitvec.h>

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

%include <openbabel/fingerprint.h>
%include <openbabel/descriptor.h>
%include <openbabel/forcefield.h>

%include <openbabel/op.h>

// The following %ignores avoid warning messages due to shadowed classes.
// This does not imply a loss of functionality as (in this case)
// the shadowed class is identical (from the point of view of SWIG) to
// the shadowing class.
// This is because C++ references (&) are transformed by SWIG back into
// pointers, so that OBAtomIter(OBMol &) would be treated the same as
// OBAtomIter(OBMol *).

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



%include <openbabel/obiter.h>
