%module openbabel

// For Windows, if not already set, set BABEL_DATADIR to
//     $PYTHONDIR/Lib/site-packages/openbabel_data
// (This is necessary for using FP3, for example)
%pythoncode %{
import os, sys
if sys.platform=="win32":
    if not os.environ.has_key('BABEL_DATADIR'):
        os.environ['BABEL_DATADIR'] = os.path.join(sys.prefix, "Lib",
                                                   "site-packages",
                                                   "openbabel_data")
%}

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

# These classes are renamed so that they can be replaced by Python
# classes of the same name which provide Pythonic iterators
# (see %pythoncode section below)

%rename(_OBAtomAtomIter) OpenBabel::OBAtomAtomIter;
%rename(_OBAtomBondIter) OpenBabel::OBAtomBondIter;
%rename(_OBMolAngleIter) OpenBabel::OBMolAngleIter;
%rename(_OBMolAtomIter) OpenBabel::OBMolAtomIter;
%rename(_OBMolAtomBFSIter) OpenBabel::OBMolAtomBFSIter;
%rename(_OBMolAtomDFSIter) OpenBabel::OBMolAtomDFSIter;
%rename(_OBMolBondIter) OpenBabel::OBMolBondIter;
%rename(_OBMolPairIter) OpenBabel::OBMolPairIter;
%rename(_OBMolRingIter) OpenBabel::OBMolRingIter;
%rename(_OBMolTorsionIter) OpenBabel::OBMolTorsionIter;
%rename(_OBResidueIter) OpenBabel::OBResidueIter;
%rename(_OBResidueAtomIter) OpenBabel::OBResidueAtomIter;

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

# The following class, OBiter, is subclassed to provide Python iterators
# equivalent to the C++ iterators in obiter.h

%pythoncode %{
class OBIter(object):
    OBiterator = None # This is defined by the subclasses

    def __init__(self, mol):
        self.iter = self.OBiterator(mol)
        self.finished = False

    def __iter__(self):
        return self

    def next(self):
        if not self.finished:
            b = self.iter.deref()
            self.iter.inc()
            if not self.iter.good():
                # There is nothing left to iterate over
                self.finished = True
            return b
        else:
            raise StopIteration

class OBAtomAtomIter(OBIter):
    """Iterator over the atoms attached to an atom."""
    OBiterator = _OBAtomAtomIter
class OBAtomBondIter(OBIter):
    """Iterator over the bonds attached to an atom."""
    OBiterator = _OBAtomBondIter
class OBMolAngleIter(OBIter):
    """Iterator over the angles in a molecule."""
    OBiterator = _OBMolAngleIter
class OBMolAtomIter(OBIter):
    """Iterator over the atoms in a molecule."""
    OBiterator = _OBMolAtomIter
class OBMolAtomBFSIter(OBIter):
    """Iterator over the atoms in a molecule in a breadth-first manner."""
    OBiterator = _OBMolAtomBFSIter
class OBMolAtomDFSIter(OBIter):
    """Iterator over the atoms in a molecule in a depth-first manner."""
    OBiterator = _OBMolAtomDFSIter
class OBMolBondIter(OBIter):
    """Iterator over the bonds in a molecule."""
    OBiterator = _OBMolBondIter
class OBMolPairIter(OBIter):
    """Iterator over pairs of atoms in a molecule."""
    OBiterator = _OBMolPairIter
class OBMolRingIter(OBIter):
    """Iterator over the rings in a molecule."""
    OBiterator = _OBMolRingIter
class OBMolTorsionIter(OBIter):
    """Iterator over the torsion angles in a molecule."""
    OBiterator = _OBMolTorsionIter
class OBResidueIter(OBIter):
    """Iterator over the residues in a molecule."""
    OBiterator = _OBResidueIter
class OBResidueAtomIter(OBIter):
    """Iterator over the atoms in a residue."""
    OBiterator = _OBResidueAtomIter
%}

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
        
