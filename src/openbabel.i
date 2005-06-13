%module openbabel
%{
#include "mol.h"
#include "obconversion.h"
%}
%include "std_list.i"
%include "std_map.i"
%include "std_vector.i"
%include "std_string.i"

%import "config.h"

%import "base.h"
%import "data.h"
%import "chains.h"
%import "math/vector3.h"
%import "bitvec.h"
%import "ring.h"
%import "generic.h"
%import "typer.h"
%import "oberror.h"

%include "obconversion.h"
%include "mol.h"
