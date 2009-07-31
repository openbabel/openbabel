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
#include <openbabel/griddata.h>
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

#include <openbabel/data.h>
#include <openbabel/parsmart.h>
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
%ignore *::DescribeBits;
%import <openbabel/babelconfig.h>

%include <openbabel/data.h>
%include <openbabel/rand.h>
%include <openbabel/obutil.h>
%include <openbabel/math/vector3.h>
%import <openbabel/math/matrix3x3.h>

%import <openbabel/math/spacegroup.h>
%include <openbabel/base.h>
%include <openbabel/generic.h>
%include <openbabel/griddata.h>

%import <openbabel/chains.h>
//# %import <openbabel/bitvec.h>
%import <openbabel/typer.h>

%include <openbabel/plugin.h>

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

%include <openbabel/fingerprint.h>
%include <openbabel/descriptor.h>
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
%ignore OBMolAtomBFSIter;
%ignore OBMolAtomDFSIter;
%ignore OBMolBondIter(OBMol &);
%ignore OBMolPairIter(OBMol &);
%ignore OBMolRingIter(OBMol &);
%ignore OBMolTorsionIter(OBMol &);
%ignore OBResidueIter(OBMol &);
%ignore OBResidueAtomIter(OBResidue &);


%ignore *::operator=;

%include <openbabel/obiter.h>
