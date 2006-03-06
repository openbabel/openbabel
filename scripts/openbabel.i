// needed to work around bug in SWIG -- can't completely override module from command-line
#ifdef SWIGPERL
%module "Chemistry::OpenBabel"
#else
%module openbabel
#endif

%{
#include "obutil.h"
#include "math/vector3.h"

#include "mol.h"
#include "generic.h"
#include "ring.h"
#include "obconversion.h"

#include "data.h"
#include "parsmart.h"
%}

%include "std_list.i"
%include "std_map.i"
%include "std_vector.i"
%include "std_string.i"

namespace std {
%template (vectorInt)		vector<int>;
%template (vvInt)		vector< vector<int> >;
%template (vectorDouble) 	vector<double>;
%template (vVector3)		vector<OpenBabel::vector3>;

%template (vectorMol)		vector<OpenBabel::OBMol>;
%template (vectorBond)		vector<OpenBabel::OBBond>;
%template (vectorResidue)	vector<OpenBabel::OBResidue>;
%template (vectorRing)		vector<OpenBabel::OBRing>;
}

%import "babelconfig.h"

%include "data.h"
%include "obutil.h"
%include "math/vector3.h"

%import "base.h"
%import "chains.h"
// %import "math/vector3.h"
%import "bitvec.h"
%import "generic.h"
%import "typer.h"
%import "oberror.h"

%include "obconversion.h"
%include "mol.h"
%include "ring.h"
%include "parsmart.h"
