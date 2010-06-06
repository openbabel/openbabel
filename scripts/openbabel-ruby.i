%module OpenBabel

// These fields are renamed to valid constant names
%rename(U1MA) OpenBabel::OBResidueIndex::_1MA;
%rename(U1MG) OpenBabel::OBResidueIndex::_1MG;
%rename(U2MG) OpenBabel::OBResidueIndex::_2MG;
%rename(U7MG) OpenBabel::OBResidueIndex::_7MG;
%rename(U5MU) OpenBabel::OBResidueIndex::_5MU;
%rename(U5MC) OpenBabel::OBResidueIndex::_5MC;

%{
// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif


#include <openbabel/obutil.h>
#include <openbabel/rand.h>
#include <openbabel/math/vector3.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/generic.h>
#include <openbabel/griddata.h>

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
#include <openbabel/rotor.h>
#include <openbabel/rotamer.h>

%}

%include "std_map.i"
%include "std_vector.i"
%include "std_string.i"

namespace std {
%template (VectorInt)		  vector<int>;
%template (VectorUnsignedInt)     vector<unsigned int>;
%template (VVInt)		      vector< vector<int> >;
%template (VectorDouble) 	vector<double>;
%template (VectorString)		  vector<std::string>;
%template (VVector3)		  vector<OpenBabel::vector3>;

%template (VectorMol)		  vector<OpenBabel::OBMol>;
%template (VectorBond)		vector<OpenBabel::OBBond>;
%template (VectorResidue)	vector<OpenBabel::OBResidue>;
%template (VectorRing)		vector<OpenBabel::OBRing>;
%template (VectorpRing)		vector<OpenBabel::OBRing*>;
%template (VectorGenericData)    vector<OpenBabel::OBGenericData*>;
}


%inline %{
OpenBabel::OBPairData *toPairData(OpenBabel::OBGenericData *data) {
	return (OpenBabel::OBPairData *) data;
}

OpenBabel::OBUnitCell *toUnitCell(OpenBabel::OBGenericData *data) {
	return (OpenBabel::OBUnitCell *) data;
}
%}

// These methods are renamed to valid method names
%rename(inc)   *::operator++;
%rename(good)  *::operator bool;
%rename(deref) *::operator->;
%rename(add)  *::operator+=;
%rename(idx)  *::operator[];
%ignore *::DescribeBits;
%ignore *::operator=;
%ignore *::operator*=;
%ignore *::operator/=;
%ignore *::operator-=;
%ignore *::operator!=;

%import <openbabel/babelconfig.h>

%include <openbabel/data.h>
%include <openbabel/rand.h>
%include <openbabel/obutil.h>
%include <openbabel/math/vector3.h>
%warnfilter(503) OpenBabel::matrix3x3; // Not wrapping any of the overloaded operators
%include <openbabel/math/matrix3x3.h>

%import <openbabel/math/spacegroup.h>

# CloneData should be used instead of the following method
%ignore OpenBabel::OBBase::SetData;
%include <openbabel/base.h>
%include <openbabel/generic.h>
%include <openbabel/griddata.h>

%import <openbabel/chains.h>
//# %import <openbabel/bitvec.h>
%import <openbabel/typer.h>

// To avoid warning in oberror.h about "Nothing known about std::binary_function"
namespace std { 
        template <T1, T2, T3>
        class binary_function {}; 
}
%template(Dummy) std::binary_function <const char *, const char *, bool>;
%include <openbabel/plugin.h>

// To avoid warning in oberror.h about "Nothing known about std::stringbuf"
namespace std { class stringbuf {}; }
%warnfilter(503) OpenBabel::OBError; // Not wrapping any of the overloaded operators
%include <openbabel/oberror.h>
%include <openbabel/format.h>
%include <openbabel/obconversion.h>
%include <openbabel/residue.h>
%include <openbabel/internalcoord.h>
%include <openbabel/atom.h>
%include <openbabel/bond.h>
%include <openbabel/mol.h>
%include <openbabel/ring.h>
%include <openbabel/parsmart.h>
%include <openbabel/alias.h>
%include <openbabel/atomclass.h>
%ignore OpenBabel::FptIndex;
%include <openbabel/fingerprint.h>
%include <openbabel/descriptor.h>

# Ignore shadowed methods
%ignore OpenBabel::OBForceField::VectorSubtract(const double *const, const double *const, double *);
%ignore OpenBabel::OBForceField::VectorMultiply(const double *const, const double, double *);
%include <openbabel/forcefield.h>

%include <openbabel/op.h>

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


%ignore *::operator=;

%include <openbabel/obiter.h>
