// needed to work around bug in SWIG -- can't completely override module from command-line
#ifdef SWIGPERL
%module "Chemistry::OpenBabel"
#else
%module openbabel
#endif

%{
#include "mol.h"
#include "obconversion.h"
#include "data.h"
%}

%pythoncode %{
import sys
import dl
sys.setdlopenflags(sys.getdlopenflags() | dl.RTLD_GLOBAL)
%}

%include "std_list.i"
%include "std_map.i"
%include "std_vector.i"
%include "std_string.i"

%import "babelconfig.h"

%include "data.h"

%import "base.h"
%import "chains.h"
%import "math/vector3.h"
%import "bitvec.h"
%import "ring.h"
%import "generic.h"
%import "typer.h"
%import "oberror.h"

%include "obconversion.h"
%include "mol.h"
