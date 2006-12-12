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

#include <openbabel/fingerprint.h>

#include <openbabel/data.h>
#include <openbabel/parsmart.h>
%}

%include "std_list.i"
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
%include <openbabel/rand.h>
%include <openbabel/obutil.h>
%include <openbabel/math/vector3.h>
%import <openbabel/math/matrix3x3.h>
%include <openbabel/generic.h>

%include <openbabel/base.h>

%import <openbabel/chains.h>
%import <openbabel/bitvec.h>
%import <openbabel/typer.h>

%include <openbabel/oberror.h>
%include <openbabel/obconversion.h>
%include <openbabel/residue.h>
%include <openbabel/internalcoord.h>
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

%ignore OBMolAtomIter(OBMol &);
%ignore OBMolAtomDFSIter(OBMol &);
%ignore OBMolAtomBFSIter(OBMol &);
%ignore OBMolBondIter(OBMol &);
%ignore OBResidueIter(OBMol &);
%ignore OBAtomAtomIter(OBAtom &);
%ignore OBAtomBondIter(OBAtom &);
%ignore OBResidueAtomIter(OBResidue &);

%rename(inc) OpenBabel::OBMolAtomIter::operator++;
%rename(inc) OpenBabel::OBMolBondIter::operator++;
%rename(inc) OpenBabel::OBAtomAtomIter::operator++;
%rename(inc) OpenBabel::OBAtomBondIter::operator++;
%rename(good) OpenBabel::OBMolAtomIter::operator bool;
%rename(good) OpenBabel::OBMolBondIter::operator bool;
%rename(good) OpenBabel::OBAtomAtomIter::operator bool;
%rename(good) OpenBabel::OBAtomBondIter::operator bool;
%rename(good) OpenBabel::OBMolAtomIter::operator bool;
%rename(deref) OpenBabel::OBMolAtomIter::operator->;
%rename(deref) OpenBabel::OBMolBondIter::operator->;
%rename(deref) OpenBabel::OBAtomAtomIter::operator->;
%rename(deref) OpenBabel::OBAtomBondIter::operator->;

%rename(inc) OpenBabel::OBResidueIter::operator++;
%rename(inc) OpenBabel::OBResidueAtomIter::operator++;
%rename(good) OpenBabel::OBResidueIter::operator bool;
%rename(good) OpenBabel::OBResidueAtomIter::operator bool;
%rename(deref) OpenBabel::OBResidueAtomIter::operator->;
%rename(deref) OpenBabel::OBResidueBondIter::operator->;

%include <openbabel/obiter.h>

%include "carrays.i"
%array_class(double, doubleArray)
%pythoncode %{
def double_array(mylist):
    """Create a C array of doubles from a list."""
    c = doubleArray(len(mylist))
    for i,v in enumerate(mylist):
        c[i] = v
    return c
%}
        
