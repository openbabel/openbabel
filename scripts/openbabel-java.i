%module openbabel_java

%{
// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif


#include <openbabel/obutil.h>
#include <openbabel/rand.h>
#include <openbabel/math/vector3.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/math/transform3d.h>
#include <openbabel/math/spacegroup.h>

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

#include <openbabel/data.h>
#include <openbabel/parsmart.h>




%}

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

%ignore *::operator=;
%ignore *::operator*;
%ignore *::operator*=;
%ignore *::operator+;
%ignore *::operator+=;
%ignore *::operator-;
%ignore *::operator-=;
%ignore *::operator++;
%ignore *::operator--;
%ignore *::operator/;
%ignore *::operator/=;
%ignore *::operator bool;

%import <openbabel/babelconfig.h>

%include <openbabel/data.h>
%include <openbabel/rand.h>
%include <openbabel/obutil.h>
%include <openbabel/math/vector3.h>
%include <openbabel/math/matrix3x3.h>
%include <openbabel/math/transform3d.h>
%include <openbabel/math/spacegroup.h>

%include <openbabel/base.h>
%include <openbabel/generic.h>
%include <openbabel/griddata.h>

%include <openbabel/chains.h>
%include <openbabel/bitvec.h>
%include <openbabel/typer.h>

%include <openbabel/plugin.h>

%include <openbabel/oberror.h>
%include <openbabel/format.h>
%include <openbabel/obconversion.h>
%include <openbabel/residue.h>
%include <openbabel/internalcoord.h>

%typemap(javacode) OpenBabel::OBAtom
%{
    private int currentDepth = 0;
    public void SetCurrentDepth(int d) {currentDepth = d;}
    public int GetCurrentDepth() {return currentDepth;}
%}
%include <openbabel/atom.h>
%include <openbabel/bond.h>
%include <openbabel/mol.h>
%include <openbabel/ring.h>
%include <openbabel/parsmart.h>

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
%ignore OBMolAtomBFSIter(OBMol &, int);
%ignore OBMolAtomBFSIter(OBMol &);
%ignore OBMolAtomDFSIter(OBMol &, int);
%ignore OBMolAtomDFSIter(OBMol &);
%ignore OBMolBondIter(OBMol &);
%ignore OBMolPairIter(OBMol &);
%ignore OBMolRingIter(OBMol &);
%ignore OBMolTorsionIter(OBMol &);
%ignore OBResidueIter(OBMol &);
%ignore OBResidueAtomIter(OBResidue &);

%ignore SpaceGroup::RegisterSpaceGroup;

%define WRAPITERATOR(NAME, RETURNS, OPTIONAL)
 %rename(_next) OpenBabel::NAME::next;
 %rename(inc) OpenBabel::NAME::operator++;
 %rename(HasMore) OpenBabel::NAME::operator bool;
 %typemap(javainterfaces) OpenBabel::NAME "Iterable<RETURNS>, Iterator<RETURNS>";
 %typemap(javaimports) OpenBabel::NAME "  import org.openbabel.RETURNS;
  import java.util.Iterator;"
 %typemap(javacode) OpenBabel::NAME
 %{
	public Iterator<RETURNS> iterator() {
		return this;
	}

	public boolean hasNext() {
		return HasMore();
	}

	public RETURNS next() {
		RETURNS atom = __ref__();
		OPTIONAL
		inc();
		return atom;
	}
	
	public void remove() {
	}	
 %}
%enddef

WRAPITERATOR(OBMolAtomIter, OBAtom, );
WRAPITERATOR(OBMolAtomDFSIter, OBAtom, );
WRAPITERATOR(OBMolAtomBFSIter, OBAtom, atom.SetCurrentDepth(CurrentDepth()););
WRAPITERATOR(OBMolBondIter, OBBond, );
WRAPITERATOR(OBMolAngleIter, vectorUnsignedInt, );
WRAPITERATOR(OBAtomAtomIter, OBAtom, )
WRAPITERATOR(OBAtomBondIter, OBBond, );
WRAPITERATOR(OBMolRingIter, OBRing, );
WRAPITERATOR(OBMolTorsionIter, vectorUnsignedInt, );
WRAPITERATOR(OBResidueIter, OBResidue, );
WRAPITERATOR(OBResidueAtomIter, OBAtom, );
WRAPITERATOR(OBMolPairIter, vectorUnsignedInt, );

%include <openbabel/obiter.h>

