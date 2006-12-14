%module openbabel

%{
// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <math.h>

#include <openbabel/obutil.h>
#include <openbabel/math/vector3.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/generic.h>

#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/ring.h>
#include <openbabel/obconversion.h>
#include <openbabel/oberror.h>

#include <openbabel/fingerprint.h>

#include <openbabel/data.h>
#include <openbabel/parsmart.h>

#include <openbabel/bitvec.h>
#include <openbabel/chains.h>
#include <openbabel/typer.h>
%}

%include "std_map.i"
%include "std_vector.i"
%include "std_string.i"

namespace std {
%template (vectorInt)		  vector<int>;
%template (vectorUnsignedInt)     vector<unsigned int>;
%template (vvInt)		      vector< vector<int> >;
%template (vectorDouble) 	vector<double>;
%template (vVector3)		  vector<OpenBabel::vector3>;

%template (vectorMol)		  vector<OpenBabel::OBMol>;
%template (vectorBond)		vector<OpenBabel::OBBond>;
%template (vectorResidue)	vector<OpenBabel::OBResidue>;
%template (vectorRing)		vector<OpenBabel::OBRing>;
%template (vectorData)    vector<OpenBabel::OBGenericData*>;
}

%import <openbabel/babelconfig.h>

%apply std::string &OUTPUT { std::string &to };
%include <openbabel/data.h>
%include <openbabel/obutil.h>
%include <openbabel/math/vector3.h>
%include <openbabel/math/matrix3x3.h>
%include <openbabel/generic.h>

%include <openbabel/base.h>

%include <openbabel/chains.h>
%include <openbabel/bitvec.h>
%include <openbabel/typer.h>

%include <openbabel/oberror.h>
%include <openbabel/obconversion.h>
%include <openbabel/residue.h>
%include <openbabel/atom.h>
%include <openbabel/bond.h>
%include <openbabel/mol.h>
%include <openbabel/ring.h>
%include <openbabel/parsmart.h>
%include <openbabel/fingerprint.h>

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
%ignore OBMolBondIter(OBMol &);
%ignore OBMolPairIter(OBMol &);
%ignore OBMolRingIter(OBMol &);
%ignore OBMolTorsionIter(OBMol &);
%ignore OBResidueIter(OBMol &);
%ignore OBResidueAtomIter(OBResidue &);

# These methods are renamed to valid Python method names, as otherwise
# they cannot be used from Python

%rename(inc) OpenBabel::OBAtomAtomIter::operator++;
%rename(inc) OpenBabel::OBAtomBondIter::operator++;
%rename(inc) OpenBabel::OBMolAngleIter::operator++;
%rename(inc) OpenBabel::OBMolAtomIter::operator++;
%rename(inc) OpenBabel::OBMolAtomBFSIter::operator++;
%rename(inc) OpenBabel::OBMolAtomDFSIter::operator++;
%rename(inc) OpenBabel::OBMolBondIter::operator++;
%rename(inc) OpenBabel::OBMolPairIter::operator++;
%rename(inc) OpenBabel::OBMolRingIter::operator++;
%rename(inc) OpenBabel::OBMolTorsionIter::operator++;
%rename(inc) OpenBabel::OBResidueIter::operator++;
%rename(inc) OpenBabel::OBResidueAtomIter::operator++;

%rename(good) OpenBabel::OBAtomAtomIter::operator bool;
%rename(good) OpenBabel::OBAtomBondIter::operator bool;
%rename(good) OpenBabel::OBMolAngleIter::operator bool;
%rename(good) OpenBabel::OBMolAtomIter::operator bool;
%rename(good) OpenBabel::OBMolAtomBFSIter::operator bool;
%rename(good) OpenBabel::OBMolAtomDFSIter::operator bool;
%rename(good) OpenBabel::OBMolBondIter::operator bool;
%rename(good) OpenBabel::OBMolPairIter::operator bool;
%rename(good) OpenBabel::OBMolRingIter::operator bool;
%rename(good) OpenBabel::OBMolTorsionIter::operator bool;
%rename(good) OpenBabel::OBResidueIter::operator bool;
%rename(good) OpenBabel::OBResidueAtomIter::operator bool;

%rename(deref) OpenBabel::OBAtomAtomIter::operator->;
%rename(deref) OpenBabel::OBAtomBondIter::operator->;
%rename(deref) OpenBabel::OBMolAngleIter::operator->;
%rename(deref) OpenBabel::OBMolAtomIter::operator->;
%rename(deref) OpenBabel::OBMolAtomBFSIter::operator->;
%rename(deref) OpenBabel::OBMolAtomDFSIter::operator->;
%rename(deref) OpenBabel::OBMolBondIter::operator->;
%rename(deref) OpenBabel::OBMolPairIter::operator->;
%rename(deref) OpenBabel::OBMolRingIter::operator->;
%rename(deref) OpenBabel::OBMolTorsionIter::operator->;
%rename(deref) OpenBabel::OBResidueIter::operator->;
%rename(deref) OpenBabel::OBResidueAtomIter::operator->;

%include <openbabel/obiter.h>
