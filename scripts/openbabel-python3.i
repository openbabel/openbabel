%module obtemplate

%{
// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/residue.h>
#include <openbabel/ring.h>

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
%template (vectorString)		  vector<std::string>;
%template (vVector3)		  vector<OpenBabel::vector3>;

%template (vectorMol)		  vector<OpenBabel::OBMol>;
%template (vectorBond)		vector<OpenBabel::OBBond>;
%template (vectorResidue)	vector<OpenBabel::OBResidue>;
%template (vectorRing)		vector<OpenBabel::OBRing>;
%template (vectorpRing)		vector<OpenBabel::OBRing*>;
%template (vectorData)    vector<OpenBabel::OBGenericData*>;
}

// These methods are renamed to valid Python method names, as otherwise
// they cannot be used from Python
%rename(inc)   *::operator++;
%rename(good)  *::operator bool;

%import <openbabel/babelconfig.h>

%import <openbabel/base.h>
%import <openbabel/residue.h>
%import <openbabel/atom.h>
%import <openbabel/bond.h>
%import <openbabel/mol.h>
%import <openbabel/ring.h>

%ignore *::operator=;

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

