%module openbabel_R
%{
#include <R.h>
// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif


#include <openbabel/obutil.h>
#include <openbabel/math/vector3.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/math/transform3d.h>
#include <openbabel/math/spacegroup.h>

#include <openbabel/generic.h>
#include <openbabel/griddata.h>

#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/reaction.h>
#include <openbabel/residue.h>
#include <openbabel/internalcoord.h>

#include <openbabel/ring.h>
#include <openbabel/obconversion.h>
#include <openbabel/oberror.h>
#include <openbabel/plugin.h>
#include <openbabel/fingerprint.h>
#include <openbabel/descriptor.h>
#include <openbabel/format.h>

#include <openbabel/builder.h>
#include <openbabel/op.h>

#include <openbabel/bitvec.h>
#include <openbabel/data.h>
#include <openbabel/parsmart.h>
#include <openbabel/alias.h>

#include <openbabel/kinetics.h>
#include <openbabel/rotor.h>
#include <openbabel/rotamer.h>
#include <openbabel/spectrophore.h>

#include <openbabel/chargemodel.h>
#include <openbabel/graphsym.h>
#include <openbabel/isomorphism.h>
#include <openbabel/query.h>
#include <openbabel/canon.h>

#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/squareplanar.h>
#include <openbabel/stereo/bindings.h>

#include <openbabel/chains.h>
#include <openbabel/obiter.h>
%}


#ifdef HAVE_EIGEN
%{
#include <openbabel/forcefield.h>
#include <openbabel/conformersearch.h>
#include <openbabel/math/align.h>
%}
#endif

%include "std_map.i"
%include "std_vector.i"
%include "std_string.i"
%include "std_pair.i"
%include "cpointer.i"

%pointer_class(std::string,stringp)
%typemap("rtype")  const std::string & "character"; //add typemap for const references
%typemap("rtype")  std::string *, std::string &; // remove erroneous typemap in swig libraries

%typemap(scoercein) const std::string &
      %{  $input = as($input, "character");     %}
%typemap(scoercein) std::string *, std::string &; 



%inline %{

    #include <iostream>
    std::istream* istreamFromString(std::string str){
       return new std::istringstream(str);
    }
    std::istream* istreamFromFile(const char *filename) {
       return new std::ifstream(filename);
    }
    std::ostream* ostreamToString(){
       return new std::ostringstream();
    }
    std::ostream* ostreamToFile(const char *filename) {
       return new std::ofstream(filename);
    }
    const std::string stringFromOstream(const std::ostringstream* os){
     // return (os->str()).c_str();
      return os->str();
    }

%}


namespace std {

%define VVTEMPLATE_WRAP(name, T) 
%feature("ignore") vector< vector<T> >::append;
%feature("ignore") vector< vector<T> >::assign;
%feature("ignore") vector< vector<T> >::back;
%feature("ignore") vector< vector<T> >::begin;
%feature("ignore") vector< vector<T> >::capacity;
%feature("ignore") vector< vector<T> >::clear;
%feature("ignore") vector< vector<T> >::empty;
%feature("ignore") vector< vector<T> >::end;
%feature("ignore") vector< vector<T> >::erase;
%feature("ignore") vector< vector<T> >::front;
%feature("ignore") vector< vector<T> >::get_allocator;
%feature("ignore") vector< vector<T> >::insert;
%feature("ignore") vector< vector<T> >::pop;
%feature("ignore") vector< vector<T> >::pop_back;
%feature("ignore") vector< vector<T> >::push_back;
%feature("ignore") vector< vector<T> >::rbegin;
%feature("ignore") vector< vector<T> >::rend;
%feature("ignore") vector< vector<T> >::reserve;
%feature("ignore") vector< vector<T> >::resize;
%feature("ignore") vector< vector<T> >::swap;
%template(vectorv ## name) vector< vector<T> >;
%enddef

%define VECTORTEMPLATE_WRAP(vectorname, T) 
%feature("ignore") vector<T>::append;
%feature("ignore") vector<T>::assign;
%feature("ignore") vector<T>::back;
%feature("ignore") vector<T>::begin;
%feature("ignore") vector<T>::capacity;
%feature("ignore") vector<T>::clear;
%feature("ignore") vector<T>::empty;
%feature("ignore") vector<T>::end;
%feature("ignore") vector<T>::erase;
%feature("ignore") vector<T>::front;
%feature("ignore") vector<T>::get_allocator;
%feature("ignore") vector<T>::insert;
%feature("ignore") vector<T>::pop;
%feature("ignore") vector<T>::pop_back;
%feature("ignore") vector<T>::push_back;
%feature("ignore") vector<T>::rbegin;
%feature("ignore") vector<T>::rend;
%feature("ignore") vector<T>::reserve;
%feature("ignore") vector<T>::resize;
%feature("ignore") vector<T>::swap;
%template(vector ## vectorname) vector<T>;
%enddef

%define VECTORPAIRTEMPLATE_WRAP(vectorname, T1, T2) 
%feature("ignore") vector< pair<T1, T2> >::append;
%feature("ignore") vector< pair<T1, T2> >::assign;
%feature("ignore") vector< pair<T1, T2> >::back;
%feature("ignore") vector< pair<T1, T2> >::begin;
%feature("ignore") vector< pair<T1, T2> >::capacity;
%feature("ignore") vector< pair<T1, T2> >::clear;
%feature("ignore") vector< pair<T1, T2> >::empty;
%feature("ignore") vector< pair<T1, T2> >::end;
%feature("ignore") vector< pair<T1, T2> >::erase;
%feature("ignore") vector< pair<T1, T2> >::front;
%feature("ignore") vector< pair<T1, T2> >::get_allocator;
%feature("ignore") vector< pair<T1, T2> >::insert;
%feature("ignore") vector< pair<T1, T2> >::pop;
%feature("ignore") vector< pair<T1, T2> >::pop_back;
%feature("ignore") vector< pair<T1, T2> >::push_back;
%feature("ignore") vector< pair<T1, T2> >::rbegin;
%feature("ignore") vector< pair<T1, T2> >::rend;
%feature("ignore") vector< pair<T1, T2> >::reserve;
%feature("ignore") vector< pair<T1, T2> >::resize;
%feature("ignore") vector< pair<T1, T2> >::swap;
%template(vpair ## vectorname) vector< pair<T1, T2> >;
%enddef

VECTORTEMPLATE_WRAP(Int, int)
VECTORTEMPLATE_WRAP(UnsignedInt, unsigned int)
VVTEMPLATE_WRAP(Int, int)
VECTORTEMPLATE_WRAP(Double, double)
VECTORTEMPLATE_WRAP(String, std::string)
VECTORTEMPLATE_WRAP(Vector3, OpenBabel::vector3)
VVTEMPLATE_WRAP(Vector3, OpenBabel::vector3)
VECTORTEMPLATE_WRAP(OBMol, OpenBabel::OBMol)
VECTORTEMPLATE_WRAP(OBBond, OpenBabel::OBBond)
VECTORTEMPLATE_WRAP(OBResidue, OpenBabel::OBResidue)
VECTORTEMPLATE_WRAP(OBRing, OpenBabel::OBRing)
VECTORTEMPLATE_WRAP(pOBRing, OpenBabel::OBRing*)
VECTORTEMPLATE_WRAP(pOBGenericData, OpenBabel::OBGenericData*)
VECTORTEMPLATE_WRAP(pOBInternalCoord, OpenBabel::OBInternalCoord*)

%template(pairUIntUInt) pair<unsigned int, unsigned int>;
VECTORPAIRTEMPLATE_WRAP(UIntUInt, unsigned int, unsigned int);
%template(vvpairUIntUInt) vector< vector< pair<unsigned int, unsigned int> > >;
}

%define CAST_GENERICDATA_TO(subclass)
%inline %{
OpenBabel::OB ## subclass *to ## subclass(OpenBabel::OBGenericData *data) {
    return (OpenBabel::OB ## subclass *) data;
}
%}
%enddef
%inline %{ // cant use macro -- AliasData not OBAliasData
OpenBabel::AliasData *toAliasData(OpenBabel::OBGenericData *data) {
    return (OpenBabel::AliasData*) data;
}
%}
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

%ignore *::operator=;
%ignore *::operator==;
%ignore *::operator!=;
%ignore *::operator<<;
%ignore *::operator();
%ignore *::operator*;
%ignore *::operator*=;
%ignore *::operator+;
%ignore *::operator+=;
%ignore *::operator-;
%ignore *::operator-=;
%ignore *::operator++;
%ignore *::operator--;
%ignore *::operator/;
%ignore *::operator/=;
%ignore *::operator bool;
%ignore *::operator[];

%import <openbabel/babelconfig.h>

%warnfilter(516) OpenBabel::OBElementTable; // Ignoring std::string methods in favour of char* ones
%include <openbabel/data.h>
%include <openbabel/obutil.h>
%warnfilter(516) OpenBabel::vector3; // Using the const x(), y() and z() in favour of the non-const
%include <openbabel/math/vector3.h>
%warnfilter(503) OpenBabel::matrix3x3; // Not wrapping any of the overloaded operators
%include <openbabel/math/matrix3x3.h>
%include <openbabel/math/transform3d.h>
%warnfilter(516) OpenBabel::SpaceGroup; // Ignoring std::string methods in favour of char* ones
%include <openbabel/math/spacegroup.h>
%warnfilter(503) OpenBabel::OBBitVec; // Not wrapping any of the overloaded operators
%include <openbabel/bitvec.h>

// CloneData should be used instead of the following method
%ignore OpenBabel::OBBase::SetData;
%ignore OpenBabel::OBBase::GetData(char const *);
%ignore OpenBabel::OBBase::HasData(char const *);
%include <openbabel/base.h>

%warnfilter(516) OpenBabel::OBPairData; // Ignoring std::string methods in favour of char* ones
%warnfilter(516) OpenBabel::OBSetData;
%warnfilter(516) OpenBabel::OBCommentData;
%warnfilter(516) OpenBabel::OBUnitCell; // Ignoring const methods in favour of non-const ones
%include <openbabel/generic.h>
%include <openbabel/griddata.h>

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
%ignore OpenBabel::OBConversion::FindFormat(const char *);
%ignore OpenBabel::OBConversion::FormatFromExt(const char *);
%include <openbabel/obconversion.h>
%include <openbabel/residue.h>
%include <openbabel/internalcoord.h>

%typemap(javacode) OpenBabel::OBAtom
%{
    private int currentDepth = 0;
    public void SetCurrentDepth(int d) {currentDepth = d;}
    public int GetCurrentDepth() {return currentDepth;}
%}
%ignore OpenBabel::OBAtom::GetVector() const;
%ignore OpenBabel::OBAtom::IsHydrogen const;
%ignore OpenBabel::OBAtom::SetType(char const *);
%include <openbabel/atom.h>
%ignore OpenBabel::OBBond::GetBeginAtom() const;
%ignore OpenBabel::OBBond::GetEndAtom() const;
%include <openbabel/bond.h>
%include <openbabel/reaction.h>

// Remove C++ iterators
%define IGNORE_ITER(parent, iteree)
%ignore OpenBabel::parent::Begin ## iteree ## s;
%ignore OpenBabel::parent::End ## iteree ## s;
%ignore OpenBabel::parent::Begin ## iteree;
%ignore OpenBabel::parent::Next ## iteree;
%enddef
IGNORE_ITER(OBMol, Bond)
IGNORE_ITER(OBMol, Atom)
IGNORE_ITER(OBMol, Residue)
%include <openbabel/mol.h>
%include <openbabel/ring.h>
%warnfilter(516) OpenBabel::OBSmartsPattern; // Using non-const versions of GetSMARTS
%include <openbabel/parsmart.h>
%warnfilter(516) OpenBabel::AliasData; // Ignoring std::string methods in favour of char* ones
%include <openbabel/alias.h>
%ignore OpenBabel::FptIndex;
%include <openbabel/fingerprint.h>
%ignore OpenBabel::OBDescriptor::LessThan;
%include <openbabel/descriptor.h>
%include <openbabel/builder.h>
%include <openbabel/op.h>

%include <openbabel/chargemodel.h>
%include <openbabel/graphsym.h>
%include <openbabel/isomorphism.h>
%include <openbabel/query.h>
%include <openbabel/canon.h>
%include <openbabel/stereo/stereo.h>

// Ignore shadowed method
%ignore OpenBabel::OBRotor::GetRotAtoms() const;
%include <openbabel/rotor.h>
%ignore OpenBabel::Swab;
%include <openbabel/rotamer.h>
%include <openbabel/spectrophore.h>
#ifdef HAVE_EIGEN

// Ignore shadowed methods
%ignore OpenBabel::OBForceField::VectorSubtract(const double *const, const double *const, double *);
%ignore OpenBabel::OBForceField::VectorMultiply(const double *const, const double, double *);
%warnfilter(516) OpenBabel::OBForceField; // Ignoring std::string methods in favour of char* ones
%include <openbabel/forcefield.h>


%include <openbabel/conformersearch.h>
%include <openbabel/math/align.h>
#endif

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
%ignore OBMolBondBFSIter(OBMol &);
%ignore OBMolBondBFSIter(OBMol &, int);
%ignore OBMolPairIter(OBMol &);
%ignore OBMolRingIter(OBMol &);
%ignore OBMolTorsionIter(OBMol &);
%ignore OBResidueIter(OBMol &);
%ignore OBResidueAtomIter(OBResidue &);

%ignore SpaceGroup::RegisterSpaceGroup;

%define WRAPITERATOR(NAME, RETURNS, OPTIONAL)
 %rename(_next) OpenBabel::NAME::next;
 %rename(inc) OpenBabel::NAME::operator++;
 %rename(HasMore) OpenBabel::NAME::operator bool;
 %typemap(javainterfaces) OpenBabel::NAME "Iterable<RETURNS>, Iterator<RETURNS>";
 %typemap(javaimports) OpenBabel::NAME "  import org.openbabel.RETURNS;
  import java.util.Iterator;"
 %typemap(javacode) OpenBabel::NAME
 %{
	public Iterator<RETURNS> iterator() {
		return this;
	}

	public boolean hasNext() {
		return HasMore();
	}

	public RETURNS next() {
		RETURNS atom = __ref__();
		OPTIONAL
		inc();
		return atom;
	}
	
	public void remove() {
	}	
 %}
%enddef

WRAPITERATOR(OBMolAtomIter, OBAtom, );
WRAPITERATOR(OBMolAtomDFSIter, OBAtom, );
WRAPITERATOR(OBMolAtomBFSIter, OBAtom, atom.SetCurrentDepth(CurrentDepth()););
WRAPITERATOR(OBMolBondIter, OBBond, );
WRAPITERATOR(OBMolAngleIter, vectorUnsignedInt, );
WRAPITERATOR(OBAtomAtomIter, OBAtom, )
WRAPITERATOR(OBAtomBondIter, OBBond, );
WRAPITERATOR(OBMolRingIter, OBRing, );
WRAPITERATOR(OBMolTorsionIter, vectorUnsignedInt, );
WRAPITERATOR(OBResidueIter, OBResidue, );
WRAPITERATOR(OBResidueAtomIter, OBAtom, );
WRAPITERATOR(OBMolPairIter, vectorUnsignedInt, );

%include <openbabel/obiter.h>

%include "stereo.i"
