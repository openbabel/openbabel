%module obconversion

%{
// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/base.h>

#include <openbabel/generic.h>
#include <openbabel/griddata.h>
#include <openbabel/math/vector3.h>
#include <openbabel/bitvec.h>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/oberror.h>
#include <openbabel/plugin.h>
#include <openbabel/fingerprint.h>
#include <openbabel/descriptor.h>
#include <openbabel/format.h>
#include <openbabel/forcefield.h>
#include <openbabel/op.h>

%}

// These methods are renamed to valid Python method names, as otherwise
// they cannot be used from Python
%rename(inc)   *::operator++;
%rename(good)  *::operator bool;

%import <openbabel/babelconfig.h>

%import <openbabel/base.h>
%import <openbabel/generic.h>
%import <openbabel/griddata.h>

%import <openbabel/math/vector3.h>
%import <openbabel/bitvec.h>

%import <openbabel/chains.h>
%import <openbabel/typer.h>

%include <openbabel/plugin.h>
%include <openbabel/oberror.h>
%include <openbabel/format.h>
%include <openbabel/obconversion.h>
%import <openbabel/residue.h>
%import <openbabel/internalcoord.h>
%import <openbabel/atom.h>
%import <openbabel/bond.h>
%import <openbabel/mol.h>
%import <openbabel/ring.h>
%import <openbabel/parsmart.h>

%include <openbabel/fingerprint.h>
%include <openbabel/descriptor.h>
%include <openbabel/forcefield.h>

%include <openbabel/op.h>
%include <openbabel/bitvec.h>

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

# Functions to set the log file to std::cout and std::cerr
       
%ignore OBForceField::SetLogFile(std::ostream *pos);
%extend OpenBabel::OBForceField {
  void SetLogToStdOut() 
  {
    self->SetLogFile(&std::cout);
  }

  void SetLogToStdErr() 
  {
    self->SetLogFile(&std::cerr);
  }
};    


