%module openbabel

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
%ignore *::operator*;
%ignore *::operator !=;

%include "std_vector.i"
%include "std_string.i"


%template (vectorInt)		  std::vector<int>;
%template (vectorUnsignedInt)     std::vector<unsigned int>;
//SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(vector<int>, std::vector<int>);
//%template (vvInt)		      std::vector< std::vector<int> >;
%template (vectorDouble) 	std::vector<double>;
%template (vectorString)		  std::vector<std::string>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(vector3, OpenBabel::vector3);
%template (vVector3)		  std::vector<OpenBabel::vector3>;

SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBMol, OpenBabel::OBMol);
%template (vectorMol)		  std::vector<OpenBabel::OBMol>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBBond, OpenBabel::OBBond);
%template (vectorBond)		std::vector<OpenBabel::OBBond>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBResidue, OpenBabel::OBResidue);
%template (vectorResidue)	std::vector<OpenBabel::OBResidue>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBRing, OpenBabel::OBRing);
%template (vectorRing)		std::vector<OpenBabel::OBRing>;
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBRing*, OpenBabel::OBRing*);
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(OBGenericData*, OpenBabel::OBGenericData*);
// Uncommenting the following ...
//%template (vectorpRing)		std::vector<OpenBabel::OBRing*>;
//%template (vectorData)    std::vector<OpenBabel::OBGenericData*>;
// ... gives ...
// Error	3	Pointers and fixed size buffers may only be used in an unsafe context

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
%ignore OBMolBondIter(OBMol &);
%ignore OBMolPairIter(OBMol &);
%ignore OBMolRingIter(OBMol &);
%ignore OBMolTorsionIter(OBMol &);
%ignore OBResidueIter(OBMol &);
%ignore OBResidueAtomIter(OBResidue &);



%include <openbabel/obiter.h>
