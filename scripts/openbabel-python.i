%module openbabel

%begin %{
#define SWIG_PYTHON_2_UNICODE
%}

%{
// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/obutil.h>
#include <openbabel/math/vector3.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/math/transform3d.h>
#include <openbabel/math/spacegroup.h>

#include <openbabel/generic.h>
#include <openbabel/griddata.h>
#include <openbabel/elements.h>

#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/reaction.h>
#include <openbabel/reactionfacade.h>
#include <openbabel/residue.h>
#include <openbabel/internalcoord.h>
#include <openbabel/bondtyper.h>

#include <openbabel/ring.h>
#include <openbabel/obconversion.h>
#include <openbabel/obfunctions.h>
#include <openbabel/oberror.h>
#include <openbabel/plugin.h>
#include <openbabel/fingerprint.h>
#include <openbabel/descriptor.h>
#include <openbabel/format.h>

#include <openbabel/forcefield.h>
#include <openbabel/builder.h>
#include <openbabel/op.h>

#include <openbabel/bitvec.h>
#include <openbabel/data.h>
#include <openbabel/parsmart.h>
#include <openbabel/alias.h>

#include <openbabel/kinetics.h>
#include <openbabel/rotor.h>
#include <openbabel/rotamer.h>
#include <openbabel/spectrophore.h>

#include <openbabel/chargemodel.h>
#include <openbabel/phmodel.h>
#include <openbabel/graphsym.h>
#include <openbabel/isomorphism.h>
#include <openbabel/query.h>
#include <openbabel/canon.h>

#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/squareplanar.h>
#include <openbabel/stereo/bindings.h>

#include <openbabel/chains.h>
#include <openbabel/obiter.h>
%}

// Set and reset dlopenflags so that plugin loading works fine for "import _openbabel"
%pythonbegin %{
import sys
if sys.platform.find("linux") != -1:
    dlflags = sys.getdlopenflags()
    import ctypes
    sys.setdlopenflags(dlflags | ctypes.RTLD_GLOBAL)
%}
%pythoncode %{
if sys.platform.find("linux") != -1:
    sys.setdlopenflags(dlflags)
%}


// Ignore methods that require std::vector of OBAtom.
%ignore OpenBabel::OBMol::FindChildren(std::vector< OBAtom * > &, OBAtom *, OBAtom *);
%ignore OpenBabel::OBResidue::GetAtoms;

#ifdef HAVE_EIGEN
%{
#include <openbabel/conformersearch.h>
#include <openbabel/math/align.h>
%}
#else
%ignore OpenBabel::OBForceField::FastRotorSearch;
%ignore OpenBabel::OBForceField::DiverseConfGen;
#endif

%include "std_list.i"
%include "std_map.i"
%include "std_vector.i"
%include "std_string.i"
%include "std_pair.i"

namespace std {

%define VVTEMPLATE_WRAP(name, T)
%feature("ignore") vector< vector<T> >::append;
%feature("ignore") vector< vector<T> >::assign;
%feature("ignore") vector< vector<T> >::back;
%feature("ignore") vector< vector<T> >::begin;
%feature("ignore") vector< vector<T> >::capacity;
%feature("ignore") vector< vector<T> >::empty;
%feature("ignore") vector< vector<T> >::end;
%feature("ignore") vector< vector<T> >::erase;
%feature("ignore") vector< vector<T> >::front;
%feature("ignore") vector< vector<T> >::get_allocator;
%feature("ignore") vector< vector<T> >::insert;
%feature("ignore") vector< vector<T> >::pop;
%feature("ignore") vector< vector<T> >::pop_back;
%feature("ignore") vector< vector<T> >::push_back;
%feature("ignore") vector< vector<T> >::rbegin;
%feature("ignore") vector< vector<T> >::rend;
%feature("ignore") vector< vector<T> >::reserve;
%feature("ignore") vector< vector<T> >::resize;
%feature("ignore") vector< vector<T> >::size;
%feature("ignore") vector< vector<T> >::swap;
%template(vectorv ## name) vector< vector<T> >;
%enddef

%define VECTORTEMPLATE_WRAP(vectorname, T)
%feature("ignore") vector<T>::append;
%feature("ignore") vector<T>::assign;
%feature("ignore") vector<T>::back;
%feature("ignore") vector<T>::begin;
%feature("ignore") vector<T>::capacity;
%feature("ignore") vector<T>::empty;
%feature("ignore") vector<T>::end;
%feature("ignore") vector<T>::erase;
%feature("ignore") vector<T>::front;
%feature("ignore") vector<T>::get_allocator;
%feature("ignore") vector<T>::insert;
%feature("ignore") vector<T>::pop;
%feature("ignore") vector<T>::pop_back;
%feature("ignore") vector<T>::push_back;
%feature("ignore") vector<T>::rbegin;
%feature("ignore") vector<T>::rend;
%feature("ignore") vector<T>::reserve;
%feature("ignore") vector<T>::resize;
%feature("ignore") vector<T>::size;
%feature("ignore") vector<T>::swap;
%template(vector ## vectorname) vector<T>;
%enddef

%define VECTORPAIRTEMPLATE_WRAP(vectorname, T1, T2)
%feature("ignore") vector< pair<T1, T2> >::append;
%feature("ignore") vector< pair<T1, T2> >::assign;
%feature("ignore") vector< pair<T1, T2> >::back;
%feature("ignore") vector< pair<T1, T2> >::begin;
%feature("ignore") vector< pair<T1, T2> >::capacity;
%feature("ignore") vector< pair<T1, T2> >::empty;
%feature("ignore") vector< pair<T1, T2> >::end;
%feature("ignore") vector< pair<T1, T2> >::erase;
%feature("ignore") vector< pair<T1, T2> >::front;
%feature("ignore") vector< pair<T1, T2> >::get_allocator;
%feature("ignore") vector< pair<T1, T2> >::insert;
%feature("ignore") vector< pair<T1, T2> >::pop;
%feature("ignore") vector< pair<T1, T2> >::pop_back;
%feature("ignore") vector< pair<T1, T2> >::push_back;
%feature("ignore") vector< pair<T1, T2> >::rbegin;
%feature("ignore") vector< pair<T1, T2> >::rend;
%feature("ignore") vector< pair<T1, T2> >::reserve;
%feature("ignore") vector< pair<T1, T2> >::resize;
%feature("ignore") vector< pair<T1, T2> >::size;
%feature("ignore") vector< pair<T1, T2> >::swap;
%template(vpair ## vectorname) vector< pair<T1, T2> >;
%enddef

VECTORTEMPLATE_WRAP(Int, int)
VECTORTEMPLATE_WRAP(UnsignedInt, unsigned int)
VVTEMPLATE_WRAP(Int, int)
VECTORTEMPLATE_WRAP(Double, double)
VECTORTEMPLATE_WRAP(ULong, unsigned long)
VECTORTEMPLATE_WRAP(String, std::string)
VECTORTEMPLATE_WRAP(Vector3, OpenBabel::vector3)
VVTEMPLATE_WRAP(Vector3, OpenBabel::vector3)
VECTORTEMPLATE_WRAP(OBMol, OpenBabel::OBMol)
VECTORTEMPLATE_WRAP(OBBond, OpenBabel::OBBond)
VECTORTEMPLATE_WRAP(OBResidue, OpenBabel::OBResidue)
VECTORTEMPLATE_WRAP(OBRing, OpenBabel::OBRing)
VECTORTEMPLATE_WRAP(pOBRing, OpenBabel::OBRing*)
VECTORTEMPLATE_WRAP(pOBGenericData, OpenBabel::OBGenericData*)
VECTORTEMPLATE_WRAP(pOBInternalCoord, OpenBabel::OBInternalCoord*)

%template(pairUIntUInt) pair<unsigned int, unsigned int>;
VECTORPAIRTEMPLATE_WRAP(UIntUInt, unsigned int, unsigned int);
%template(vvpairUIntUInt) vector< vector< pair<unsigned int, unsigned int> > >;
}

%define CAST_GENERICDATA_TO(subclass)
%inline %{
OpenBabel::OB ## subclass *to ## subclass(OpenBabel::OBGenericData *data) {
    return (OpenBabel::OB ## subclass *) data;
}
%}
%enddef
%inline %{ /* can't use macro -- AliasData not OBAliasData */
OpenBabel::AliasData *toAliasData(OpenBabel::OBGenericData *data) {
    return (OpenBabel::AliasData*) data;
}
%}
CAST_GENERICDATA_TO(AngleData)
CAST_GENERICDATA_TO(CommentData)
CAST_GENERICDATA_TO(ConformerData)
CAST_GENERICDATA_TO(ExternalBondData)
CAST_GENERICDATA_TO(GridData)
CAST_GENERICDATA_TO(MatrixData)
CAST_GENERICDATA_TO(NasaThermoData)
CAST_GENERICDATA_TO(PairData)
CAST_GENERICDATA_TO(PairInteger)
CAST_GENERICDATA_TO(PairFloatingPoint)
// CAST_GENERICDATA_TO(PairTemplate)
CAST_GENERICDATA_TO(RateData)
CAST_GENERICDATA_TO(RotamerList)
CAST_GENERICDATA_TO(RotationData)
CAST_GENERICDATA_TO(SerialNums)
CAST_GENERICDATA_TO(SetData)
CAST_GENERICDATA_TO(SymmetryData)
CAST_GENERICDATA_TO(TorsionData)
CAST_GENERICDATA_TO(UnitCell)
CAST_GENERICDATA_TO(VectorData)
CAST_GENERICDATA_TO(VibrationData)
CAST_GENERICDATA_TO(VirtualBond)
CAST_GENERICDATA_TO(TetrahedralStereo)
CAST_GENERICDATA_TO(CisTransStereo)
CAST_GENERICDATA_TO(SquarePlanarStereo)

// This method is renamed to a valid Python method name, as otherwise
// it cannot be used from Python
%rename(inc)   *::operator++;
%ignore *::operator=;
%ignore *::operator[];

%import <openbabel/babelconfig.h>

%include <openbabel/data.h>
%include <openbabel/obutil.h>
%include <openbabel/math/vector3.h>
%warnfilter(503) OpenBabel::matrix3x3; // Not wrapping any of the overloaded operators
%include <openbabel/math/matrix3x3.h>
%include <openbabel/math/transform3d.h>
%include <openbabel/math/spacegroup.h>
%warnfilter(503) OpenBabel::OBBitVec; // Not wrapping any of the overloaded operators
%include <openbabel/bitvec.h>

// CloneData should be used instead of the following method
%ignore OpenBabel::OBBase::SetData;
%include <openbabel/base.h>

%include <openbabel/generic.h>
%template(OBPairInteger) OpenBabel::OBPairTemplate<int>;
%template(OBPairFloatingPoint) OpenBabel::OBPairTemplate<float>;
%include <openbabel/griddata.h>

%include <openbabel/chains.h>
%include <openbabel/typer.h>

// To avoid warning in plugin.h about "Nothing known about std::binary_function"
namespace std {
        template <T1, T2, T3>
        class binary_function {};
}
%template(dummy) std::binary_function <const char *, const char *, bool>;
%include <openbabel/plugin.h>

// To avoid warning in oberror.h about "Nothing known about std::stringbuf"
namespace std { class stringbuf {}; }
%warnfilter(503) OpenBabel::OBError; // Not wrapping any of the overloaded operators
%include <openbabel/oberror.h>
%include <openbabel/format.h>
%include <openbabel/obconversion.h>
%include <openbabel/obfunctions.h>
%include <openbabel/residue.h>
%include <openbabel/internalcoord.h>
%include <openbabel/atom.h>
%include <openbabel/bond.h>
%include <openbabel/reaction.h>
%include <openbabel/reactionfacade.h>

// Remove C++ iterators
%pythoncode %{
def exceptionIter(*args):
    raise Exception("""\nThis method can only be used from C++. To iterate from Python
use the Iter classes (OBMolAtomIter, etc.) as described at
http://openbabel.org/wiki/Python""")
%}
%define IGNORE_ITER(parent, iteree)
%ignore OpenBabel::parent::Begin ## iteree ## s;
%ignore OpenBabel::parent::End ## iteree ## s;
%ignore OpenBabel::parent::Begin ## iteree;
%ignore OpenBabel::parent::Next ## iteree;
%enddef
IGNORE_ITER(OBMol, Bond)
IGNORE_ITER(OBMol, Atom)
IGNORE_ITER(OBMol, Residue)
%include <openbabel/mol.h>
%pythoncode %{
OBMol.BeginAtoms = OBMol.EndAtoms = OBMol.BeginAtom = OBMol.EndAtom = exceptionIter
OBMol.BeginBonds = OBMol.EndBonds = OBMol.BeginBond = OBMol.EndBond = exceptionIter
OBMol.BeginResidues = OBMol.EndResidues = OBMol.BeginResidue = OBMol.EndResidue = exceptionIter
%}

%include <openbabel/ring.h>
%include <openbabel/parsmart.h>
%include <openbabel/alias.h>
%ignore OpenBabel::FptIndex;
%include <openbabel/fingerprint.h>
%ignore OpenBabel::OBDescriptor::LessThan;
%include <openbabel/descriptor.h>
// wrap GetRGB parameters
%include "typemaps.i"
%apply double *OUTPUT { double *r, double *g, double *b };
%include <openbabel/elements.h>
// void GetRGB(unsigned int atomic_number, double *r, double *g, double *b);
%clear double *r, double *g, double *b;

// Ignore shadowed methods
%ignore OpenBabel::OBForceField::VectorSubtract(const double *const, const double *const, double *);
%ignore OpenBabel::OBForceField::VectorMultiply(const double *const, const double, double *);
%include <openbabel/forcefield.h>

%include <openbabel/builder.h>
%include <openbabel/op.h>

%include <openbabel/chargemodel.h>
%apply std::string& INPUT { std::string &start } // Required for OBChemTsfm.Init
%apply std::string& INPUT { std::string &end }   // Required for OBChemTsfm.Init
%include <openbabel/phmodel.h>
%include <openbabel/graphsym.h>
%include <openbabel/isomorphism.h>
%include <openbabel/query.h>
%include <openbabel/canon.h>

%include <openbabel/stereo/stereo.h>

// Ignore shadowed method
%ignore OpenBabel::OBRotor::GetRotAtoms() const;
%include <openbabel/rotor.h>
%ignore OpenBabel::Swab;
%include <openbabel/rotamer.h>
%include <openbabel/spectrophore.h>
#ifdef HAVE_EIGEN
%include <openbabel/conformersearch.h>
%include <openbabel/math/align.h>
#endif

// The following %ignores avoid warning messages due to shadowed classes.
// This does not imply a loss of functionality as (in this case)
// the shadowed class is identical (from the point of view of SWIG) to
// the shadowing class.
// This is because C++ references (&) are transformed by SWIG back into
// pointers, so that OBAtomIter(OBMol &) would be treated the same as
// OBAtomIter(OBMol *).

%ignore OBAtomAtomIter(OBAtom &);
%ignore OBAtomBondIter(OBAtom &);
%ignore OBMolAngleIter(OBMol &);
%ignore OBMolAtomIter(OBMol &);
%ignore OBMolAtomBFSIter(OBMol &);
%ignore OBMolAtomDFSIter(OBMol &);
%ignore OBMolAtomBFSIter(OBMol &, int);
%ignore OBMolAtomDFSIter(OBMol &, int);
%ignore OBMolBondIter(OBMol &);
%ignore OBMolBondBFSIter(OBMol &);
%ignore OBMolBondBFSIter(OBMol &, int);
%ignore OBMolPairIter(OBMol &);
%ignore OBMolRingIter(OBMol &);
%ignore OBMolTorsionIter(OBMol &);
%ignore OBResidueIter(OBMol &);
%ignore OBResidueAtomIter(OBResidue &);

// SWIG treats operator-> specially (see 6.24 "Smart pointers and operator->()").
// If we leave this in, it adds
// all of the methods of the underlying object (e.g. OBAtom) to the
// iterator object, causing bloat.
%ignore OpenBabel::OBAtomAtomIter::operator->;
%ignore OpenBabel::OBAtomBondIter::operator->;
%ignore OpenBabel::OBMolAngleIter::operator->;
%ignore OpenBabel::OBMolAtomIter::operator->;
%ignore OpenBabel::OBMolAtomBFSIter::operator->;
%ignore OpenBabel::OBMolAtomDFSIter::operator->;
%ignore OpenBabel::OBMolAtomBFSIter::operator->;
%ignore OpenBabel::OBMolAtomDFSIter::operator->;
%ignore OpenBabel::OBMolBondIter::operator->;
%ignore OpenBabel::OBMolBondBFSIter::operator->;
%ignore OpenBabel::OBMolBondBFSIter::operator->;
%ignore OpenBabel::OBMolPairIter::operator->;
%ignore OpenBabel::OBMolRingIter::operator->;
%ignore OpenBabel::OBMolTorsionIter::operator->;
%ignore OpenBabel::OBResidueIter::operator->;
%ignore OpenBabel::OBResidueAtomIter::operator->;

// These classes are renamed so that they can be replaced by Python
// classes of the same name which provide Pythonic iterators
// (see %pythoncode section below)

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
%rename(_OBFingerprintIter) OpenBabel::PluginIter<OBFingerprint>;


%include <openbabel/obiter.h>

// The following class, OBiter, is subclassed to provide Python iterators
// equivalent to the C++ iterators in obiter.h and the plugin iterators

%pythoncode %{
class OBIter(object):
    OBiterator = None # This is defined by the subclasses

    def __init__(self, *params):
        self.iter = self.OBiterator(*params)
        self.finished = False
        if not self.iter.__bool__():
            self.finished = True

    def __iter__(self):
        return self

    def next(self):
        if not self.finished:
            b = self.iter.__ref__()
            self.iter.inc()
            if not self.iter.__bool__():
                # There is nothing left to iterate over
                self.finished = True
            return b
        else:
            raise StopIteration

    __next__ = next

class OBIterWithDepth(OBIter):
    def next(self):
        if not self.finished:
            b = self.iter.__ref__()
            depth = self.iter.CurrentDepth()
            self.iter.inc()
            if not self.iter.__bool__():
                # There is nothing left to iterate over
                self.finished = True
            return b, depth
        else:
            raise StopIteration

    __next__ = next

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
class OBMolAtomBFSIter(OBIterWithDepth):
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

// Copy some of the global variables in cvar into the openbabel namespace

%pythoncode %{
obErrorLog = cvar.obErrorLog
ttab = cvar.ttab
atomtyper = cvar.atomtyper
aromtyper = cvar.aromtyper
%}

// Functions to set the log file to std::cout and std::cerr

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

%extend OpenBabel::OBMol {
  void SetTorsion(int i, int j, int k, int l, double ang)
  {
    self->SetTorsion(self->GetAtom(i), self->GetAtom(j),
                     self->GetAtom(k), self->GetAtom(l), ang);
  }
};


%pythoncode %{
def exception(*args):
    raise Exception("Use OBMol.CloneData instead. OBMol.SetData is only for use from C++.")
OBMol.SetData = exception
%}

%include "stereo.i"
