%module OpenBabel
%include "std_string.i"

// These fields are renamed to valid constant names
%rename(U1MA) OpenBabel::OBResidueIndex::_1MA;
%rename(U1MG) OpenBabel::OBResidueIndex::_1MG;
%rename(U2MG) OpenBabel::OBResidueIndex::_2MG;
%rename(U7MG) OpenBabel::OBResidueIndex::_7MG;
%rename(U5MU) OpenBabel::OBResidueIndex::_5MU;
%rename(U5MC) OpenBabel::OBResidueIndex::_5MC;

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
#include <openbabel/kekulize.h>

#include <openbabel/base.h>
#include <openbabel/elements.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/reaction.h>
#include <openbabel/reactionfacade.h>
#include <openbabel/residue.h>
#include <openbabel/internalcoord.h>

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
#include <openbabel/kekulize.h>
%}

%include "std_map.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_pair.i"
%include "std_shared_ptr.i"

namespace std {

%define VVTEMPLATE_WRAP(name, T) 
%feature("ignore") vector< vector<T> >::append;
%feature("ignore") vector< vector<T> >::assign;
%feature("ignore") vector< vector<T> >::back;
%feature("ignore") vector< vector<T> >::begin;
%feature("ignore") vector< vector<T> >::capacity;
%feature("ignore") vector< vector<T> >::clear;
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
%template(VectorV ## name) vector< vector<T> >;
%enddef

%define VECTORTEMPLATE_WRAP(vectorname, T) 
%feature("ignore") vector<T>::append;
%feature("ignore") vector<T>::assign;
%feature("ignore") vector<T>::back;
%feature("ignore") vector<T>::begin;
%feature("ignore") vector<T>::capacity;
%feature("ignore") vector<T>::clear;
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
%template(Vector ## vectorname) vector<T>;
%enddef

VECTORTEMPLATE_WRAP(Int, int)
VECTORTEMPLATE_WRAP(UnsignedInt, unsigned int)
VVTEMPLATE_WRAP(Int, int)
VECTORTEMPLATE_WRAP(Double, double)
VECTORTEMPLATE_WRAP(String, std::string)
VECTORTEMPLATE_WRAP(Vector3, OpenBabel::vector3)
VECTORTEMPLATE_WRAP(OBMol, OpenBabel::OBMol)
VECTORTEMPLATE_WRAP(OBBond, OpenBabel::OBBond)
VECTORTEMPLATE_WRAP(OBResidue, OpenBabel::OBResidue)
VECTORTEMPLATE_WRAP(OBRing, OpenBabel::OBRing)
VECTORTEMPLATE_WRAP(pOBRing, OpenBabel::OBRing*)
VECTORTEMPLATE_WRAP(pOBGenericData, OpenBabel::OBGenericData*)
VECTORTEMPLATE_WRAP(pOBTetrahedralStereo, OpenBabel::OBTetrahedralStereo*)
}

%define CAST_GENERICDATA_TO(subclass)
%inline %{
OpenBabel::OB ## subclass *to ## subclass(OpenBabel::OBGenericData *data) {
    return (OpenBabel::OB ## subclass *) data;
}
%}
%enddef
CAST_GENERICDATA_TO(AngleData)
CAST_GENERICDATA_TO(CommentData)
CAST_GENERICDATA_TO(ConformerData)
CAST_GENERICDATA_TO(ExternalBondData)
CAST_GENERICDATA_TO(GridData)
CAST_GENERICDATA_TO(MatrixData)
CAST_GENERICDATA_TO(NasaThermoData)
CAST_GENERICDATA_TO(PairData)
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

// These methods are renamed to valid method names
%rename(inc)   *::operator++;
%rename(good)  *::operator bool;
%rename(deref) *::operator->;
%rename(add)  *::operator+=;
%rename(idx)  *::operator[];
%ignore *::DescribeBits;
%ignore *::operator=;
%ignore *::operator*=;
%ignore *::operator/=;
%ignore *::operator-=;
%ignore *::operator!=;
%ignore *::operator&=;
%ignore *::operator^=;
%ignore *::operator|=;

%import <openbabel/babelconfig.h>

%include <openbabel/data.h>
%include <openbabel/obutil.h>
%include <openbabel/math/vector3.h>
%warnfilter(503) OpenBabel::matrix3x3; // Not wrapping any of the overloaded operators
%include <openbabel/math/matrix3x3.h>

%import <openbabel/math/spacegroup.h>
%warnfilter(503) OpenBabel::OBBitVec; // Not wrapping any of the overloaded operators
%include <openbabel/bitvec.h>

// CloneData should be used instead of the following method
%ignore OpenBabel::OBBase::SetData;
%include <openbabel/base.h>
%include <openbabel/generic.h>
%include <openbabel/griddata.h>
%include <openbabel/kekulize.h>


%import <openbabel/chains.h>
%import <openbabel/typer.h>
%include <openbabel/plugin.h>

// To avoid warning in oberror.h about "Nothing known about std::stringbuf"
namespace std { class stringbuf {}; }
%warnfilter(503) OpenBabel::OBError; // Not wrapping any of the overloaded operators
%include <openbabel/oberror.h>
%include <openbabel/format.h>
%include <openbabel/obconversion.h>
%include <openbabel/residue.h>
%include <openbabel/internalcoord.h>
%include <openbabel/elements.h>
%include <openbabel/atom.h>
%include <openbabel/bond.h>
%include <openbabel/mol.h>
%include <openbabel/reaction.h>
%include <openbabel/reactionfacade.h>
%include <openbabel/ring.h>
%include <openbabel/parsmart.h>
%include <openbabel/alias.h>
%ignore OpenBabel::FptIndex;
%include <openbabel/fingerprint.h>
%include <openbabel/descriptor.h>

// Ignore shadowed methods
%ignore OpenBabel::OBForceField::VectorSubtract(const double *const, const double *const, double *);
%ignore OpenBabel::OBForceField::VectorMultiply(const double *const, const double, double *);
#ifdef HAVE_EIGEN3
%{
#include <openbabel/conformersearch.h>
#include <openbabel/math/align.h>
%}
#else
%ignore OpenBabel::OBForceField::FastRotorSearch;
%ignore OpenBabel::OBForceField::DiverseConfGen;
#endif
%include <openbabel/forcefield.h>

%include <openbabel/builder.h>
%include <openbabel/op.h>

%ignore OpenBabel::OBRotor::GetRotAtoms() const;
%include <openbabel/rotor.h>
%ignore OpenBabel::Swab;
%include <openbabel/rotamer.h>

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


%ignore *::operator=;

%include <openbabel/obiter.h>

%include <openbabel/chargemodel.h>
%include <openbabel/phmodel.h>
%include <openbabel/graphsym.h>
%include <openbabel/isomorphism.h>
%include <openbabel/query.h>
%include <openbabel/canon.h>

%include <openbabel/stereo/stereo.h>

%include <openbabel/kekulize.h>
%include <openbabel/obfunctions.h>

%include "stereo.i"

%{
#include <openbabel/tautomer.h>
using namespace OpenBabel;

struct RubyCallbackArgs {
    VALUE callback;
    VALUE argument;
};

static VALUE call_ruby_callback(VALUE data) {
    RubyCallbackArgs *args = reinterpret_cast<RubyCallbackArgs *>(data);
    return rb_funcall(args->callback, rb_intern("call"), 1, args->argument);
}

static OBMol *unwrap_obmol(VALUE self) {
    void *argp = 0;
    int res = SWIG_ConvertPtr(self, &argp, SWIGTYPE_p_OpenBabel__OBMol, 0);
    if (!SWIG_IsOK(res) || !argp) {
        rb_raise(rb_eTypeError, "expected OpenBabel::OBMol");
    }
    return reinterpret_cast<OBMol *>(argp);
}

// A functor class that wraps the Ruby block as the callback
class RubyTautomerFunctor : public UniqueTautomerFunctor {
private:
    VALUE ruby_callback;
    int exception_state;

public:
    // Constructor that takes the Ruby callback block (Proc or lambda)
    RubyTautomerFunctor(VALUE callback) : ruby_callback(callback), exception_state(0) {}

    // Override the operator() to call the Ruby block
    void operator()(OBMol *mol, const std::string &) override {
        if (exception_state) {
            return;
        }

        RubyCallbackArgs args = {
            ruby_callback,
            SWIG_NewPointerObj(new OBMol(*mol), SWIGTYPE_p_OpenBabel__OBMol, SWIG_POINTER_OWN)
        };
        rb_protect(call_ruby_callback, reinterpret_cast<VALUE>(&args), &exception_state);
    }

    int state() const {
        return exception_state;
    }
};

// The method to be defined on OpenBabel::OBMol in Ruby
static VALUE enumerate_tautomers(VALUE self) {
    // Return enumerator if no block provided
    if (!rb_block_given_p()) {
        return rb_enumeratorize(self, ID2SYM(rb_intern("tautomers")), 0, NULL);
    }

    RubyTautomerFunctor functor(rb_block_proc());
    EnumerateTautomers(unwrap_obmol(self), functor);

    if (functor.state()) {
        rb_jump_tag(functor.state());
    }

    return Qnil;
}

static VALUE canonical_tautomer(VALUE self) {
    OBMol *mol = new OBMol(*unwrap_obmol(self));
    CanonicalTautomer(mol);
    return SWIG_NewPointerObj(mol, SWIGTYPE_p_OpenBabel__OBMol, SWIG_POINTER_OWN);
}
%}

%init %{
    rb_define_method(SwigClassOBMol.klass, "tautomers", VALUEFUNC(enumerate_tautomers), 0);
    rb_define_method(SwigClassOBMol.klass, "canonical_tautomer", VALUEFUNC(canonical_tautomer), 0);
%}
