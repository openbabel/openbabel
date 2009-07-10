// needed to work around bug in SWIG -- can't completely override module from command-line
%module "Chemistry::OpenBabel"

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
#include <openbabel/griddata.h>

#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/residue.h>
#include <openbabel/internalcoord.h>

#include <openbabel/ring.h>
#include <openbabel/obconversion.h>
#include <openbabel/oberror.h>
#include <openbabel/plugin.h>
#include <openbabel/fingerprint.h>
#include <openbabel/descriptor.h>
#include <openbabel/format.h>

#include <openbabel/forcefield.h>
#include <openbabel/op.h>

#include <openbabel/bitvec.h>
#include <openbabel/data.h>
#include <openbabel/parsmart.h>
#include <openbabel/alias.h>
#include <openbabel/atomclass.h>
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


%inline %{
OpenBabel::OBPairData *toPairData(OpenBabel::OBGenericData *data) {
	return (OpenBabel::OBPairData *) data;
}

OpenBabel::OBUnitCell *toUnitCell(OpenBabel::OBGenericData *data) {
	return (OpenBabel::OBUnitCell *) data;
}
%}

%rename(inc)   *::operator++;
%rename(good)  *::operator bool;
%rename(deref) *::operator->;

%import <openbabel/babelconfig.h>

%include <openbabel/data.h>
%include <openbabel/rand.h>
%include <openbabel/obutil.h>
%include <openbabel/math/vector3.h>
%import <openbabel/math/matrix3x3.h>

%import <openbabel/math/spacegroup.h>
%include <openbabel/base.h>
%include <openbabel/generic.h>
%include <openbabel/griddata.h>

%import <openbabel/chains.h>
%import <openbabel/bitvec.h>
%import <openbabel/typer.h>

%import <openbabel/plugin.h>

%include <openbabel/oberror.h>
%include <openbabel/format.h>
%include <openbabel/obconversion.h>
%include <openbabel/residue.h>
%include <openbabel/internalcoord.h>
%include <openbabel/atom.h>
%include <openbabel/bond.h>
%include <openbabel/mol.h>
%include <openbabel/ring.h>
%include <openbabel/parsmart.h>
%include <openbabel/alias.h>
%include <openbabel/atomclass.h>

%include <openbabel/fingerprint.h>
%include <openbabel/descriptor.h>
%include <openbabel/forcefield.h>

%include <openbabel/op.h>

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
%ignore OBMolAtomBFSIter(OBMol &, int);
%ignore OBMolAtomDFSIter(OBMol &);
%ignore OBMolAtomDFSIter(OBMol &, int);
%ignore OBMolBondIter(OBMol &);
%ignore OBMolPairIter(OBMol &);
%ignore OBMolRingIter(OBMol &);
%ignore OBMolTorsionIter(OBMol &);
%ignore OBResidueIter(OBMol &);
%ignore OBResidueAtomIter(OBResidue &);

%ignore *::operator=;

%perlcode %{
# http://www.perl.com/pub/a/2005/06/16/iterators.html
sub gen_iterator {
    my $iterator_function = shift;
    return sub {
        my $iter = new $iterator_function(shift);

        my ($current_state, $done);

        return sub {
            # code to calculate $next_state or $done;
            return undef if $done;
            my $next_state = $iter->deref();
            $iter->inc();
            if (!$iter->good()) {$done = 1;};
            return $current_state = $next_state;
        };
    };
}
%}

%define MAKE_ITER_CLASS(itername, renamediter)
%rename itername renamediter;
%perlcode %{
package Chemistry::OpenBabel::itername;

sub new {
   my $pkg = shift;
   my $self = Chemistry::OpenBabel::gen_iterator(Chemistry::OpenBabel::renamediter)->(shift);
   bless($self, $pkg);
}
%}
%enddef

MAKE_ITER_CLASS(OBAtomAtomIter, _OBAtomAtomIter)
MAKE_ITER_CLASS(OBAtomBondIter, _OBAtomBondIter)
MAKE_ITER_CLASS(OBMolAngleIter, _OBMolAngleIter)
MAKE_ITER_CLASS(OBMolAtomIter, _OBMolAtomIter)
MAKE_ITER_CLASS(OBMolAtomBFSIter, _OBMolAtomBFSIter)
MAKE_ITER_CLASS(OBMolAtomDFSIter, _OBMolAtomDFSIter)
MAKE_ITER_CLASS(OBMolBondIter, _OBMolBondIter)
MAKE_ITER_CLASS(OBMolPairIter, _OBMolPairIter)
MAKE_ITER_CLASS(OBMolRingIter, _OBMolRingIter)
MAKE_ITER_CLASS(OBMolTorsionIter, _OBMolTorsionIter)
MAKE_ITER_CLASS(OBResidueAtomIter, _OBResidueAtomIter)
MAKE_ITER_CLASS(OBFingerprintIter, _OBFingerprintIter)

%include <openbabel/obiter.h>
