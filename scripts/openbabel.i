// needed to work around bug in SWIG -- can't completely override module from command-line
#ifdef SWIGPERL
%module "Chemistry::OpenBabel"
#else
%module openbabel
#endif

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

#include "data.h"
#include "parsmart.h"
%}

%include "std_list.i"
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

%include "data.h"
%include "obutil.h"
%include "math/vector3.h"
%import "math/matrix3x3.h"
%include "generic.h"

%include "base.h"

%import "chains.h"
%import "bitvec.h"
%import "typer.h"
%import "oberror.h"

%include "obconversion.h"
%include "mol.h"
%include "ring.h"
%include "parsmart.h"
