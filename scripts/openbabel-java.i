%module net.sourceforge.openbabel

%{
// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif


#include "obutil.h"
#include "math/vector3.h"
#include "math/matrix3x3.h"
#include "generic.h"

#include "base.h"
#include "mol.h"
#include "ring.h"
#include "obconversion.h"
#include "oberror.h"

#include "data.h"
#include "parsmart.h"
%}

%include "std_map.i"
%include "std_vector.i"
%include "std_string.i"

namespace std {
%template (vectorInt)		  vector<int>;
%template (vvInt)		      vector< vector<int> >;
%template (vectorDouble) 	vector<double>;
%template (vVector3)		  vector<OpenBabel::vector3>;

%template (vectorMol)		  vector<OpenBabel::OBMol>;
%template (vectorBond)		vector<OpenBabel::OBBond>;
%template (vectorResidue)	vector<OpenBabel::OBResidue>;
%template (vectorRing)		vector<OpenBabel::OBRing>;
%template (vectorData)    vector<OpenBabel::OBGenericData*>;
}

%import "babelconfig.h"

%apply std::string &OUTPUT { std::string &to };
%include "data.h"
%include "obutil.h"
%include "math/vector3.h"
%import "math/matrix3x3.h"
%include "generic.h"

%include "base.h"

%import "chains.h"
%import "bitvec.h"
%import "typer.h"

%include "oberror.h"
%include "obconversion.h"
%include "residue.h"
%include "atom.h"
%include "bond.h"
%include "mol.h"
%include "ring.h"
%include "parsmart.h"

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

%include "obiter.h"
