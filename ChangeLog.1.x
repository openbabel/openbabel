2004-02-10  Geoff Hutchison  <ghutchis@wso.williams.edu>

	* configure, configure.in: Update with notes on versioning for
	shared library support.

	* src/Makefile.ac, src/Makefile.in: Apply specific 0:0:0 current
	version for shared libopenbabel -- change this as appropriate for
	new versions. See libtool manual for more information on what
	versions are needed for revisions, backwards incompatibility, etc.
	
2004-02-18  Michael Banck  <mbanck@gmx.net>

	* src/pqs.cpp: Include <strings.h>, if available.

2004-02-17  Michael Banck  <mbanck@gmx.net>

	* configure.in (AC_INIT): Changed PACKAGE_NAME to 'OpenBabel'.
	Additionally specifiy PACKAGE_TARNAME as 'openbabel'.
	(AM_INIT_AUTOMAKE): Removed deprecated arguments PACKAGE and
	VERSION.
	(AC_DEFINE_DIR): Use AC_PACKAGE_TARNAME instead of
	AC_PACKAGE_NAME.
	* configure: Regenerated.

2004-02-10  Geoff Hutchison  <ghutchis@wso.williams.edu>

	* src/mol.cpp: Documentation updates.

	* NEWS, THANKS: Updates for 1.100.2 release.

	* aromatic.txt, aromatic.txt: Minor updates for aromaticity of
	various heteroatoms.

	* element.txt: * Minor documentation update.

2004-02-14  Michael Banck  <mbanck@gmx.net>

	* src/base.cpp [HAVE_CONFIG_H]: Removed the conditional around
	#include babelconfig.h (reverted patch).
	* src/bgf.cpp, src/data.cpp, src/gromos96.cpp, src/main.cpp,
	src/newmain.cpp, src/oberror.cpp, src/obutil.cpp, src/patty.cpp,
	src/pdb.cpp, src/phmodel.cpp, src/povray.cpp, src/residue.cpp,
	src/rotor.cpp, src/typer.cpp, src/base.h, src/bitvec.h, src/crk.h,
	src/data.h, src/fileformat.h, src/grid.h, src/mol.h,
	src/obifstream.h, src/obutil.h: Likewise.
	
2004-02-13  Michael Banck  <mbanck@gmx.net>

	* configure.in (AM_INIT_AUTOMAKE): Added 'no-define', to prevent
	defines for PACKAGE and VERSION in src/babelconfig.h.in, which
	would confuse other projects.
	* configure: Regenerated.
	* src/babelconfig.h.in: Regenerated.

2004-02-12 Fabien Fontaine <ffontaine@imim.es>

	* tools/obgrep.cpp: -t NUM option added to control 
	the number of match in each molecule
	*tools/obprop.cpp: new tools added to compute and print 
	simple molecular properties for non C++ programmers
	
2004-02-10  Geoff Hutchison  <ghutchis@wso.williams.edu>

	* configure, configure.in: Change --enable-shared to the
	default. Build shared libraries except where specified by the
	user.
	
2004-02-10  Michael Banck  <mbanck@gmx.net>

	* test/cmltest/Makefile.am (EXTRA_DIST): Removed auto-generated 
	.cml files.
	* test/cmltest/Makefile: Regenerated.

2004-02-09  Michael Banck  <mbanck@gmx.net>

	* src/base.cpp [HAVE_CONFIG_H]: New conditional, wrap babelconfig.h
	#include into it.
	* src/bgf.cpp, src/data.cpp, src/gromos96.cpp, src/main.cpp,
	src/newmain.cpp, src/oberror.cpp, src/obutil.cpp, src/patty.cpp,
	src/pdb.cpp, src/phmodel.cpp, src/povray.cpp, src/residue.cpp,
	src/rotor.cpp, src/typer.cpp, src/base.h, src/bitvec.h, src/crk.h,
	src/data.h, src/fileformat.h, src/grid.h, src/mol.h,
	src/obifstream.h, src/obutil.h: Likewise.

2004-02-09  Michael Banck  <mbanck@gmx.net>

	* src/base.cpp: Include babelconfig.h.
	* src/base.h: Likewise.
	* src/crk.h: Likewise.
	* src/data.h: Likewise.
	* src/fileformat.h: Likewise.
	* src/main.cpp: Likewise.
	* src/oberror.cpp: Likewise.
	* src/residue.cpp: Likewise.

2004-02-08  Geoff Hutchison  <ghutchis@wso.williams.edu>

	* src/mol.cpp, src/qchem.cpp, src/gaussian.cpp, src/report.cpp:
	Improved support for spin multiplicity and total charge --
	calculate from atomic formal charges and spin multiplicities where
	not specified.

	* Doxyfile, src/base.cpp, src/base.h, src/binary.cpp,
	src/binary.h, src/chains.cpp, src/chains.h, src/data.h,
	src/generic.cpp, src/generic.h, src/grid.h, src/mol.h,
	src/oberror.h, src/obifstream.h, src/obutil.h, src/parsmart.h,
	src/parsmi.cpp, src/pdb.cpp, src/phmodel.cpp, src/phmodel.h,
	src/residue.cpp, src/rotor.cpp, src/typer.cpp, src/typer.h:
	Updated documentation for generation by doxygen.

2004-02-07  Michael Banck  <mbanck@gmx.net>

	* src/math/matrix3x3.h [ __sgi]: Removed and...
	[HAVE_IOSTRAM], [HAVE_FSTREAM]: Replaced with this.
	* src/math/vector3.h: Likewise.
	* src/base.cpp [#include <iostream>]: Removed and...
	[HAVE_IOSTREAM, HAVE_IOSTREAM_H]: Replaced with this.
	* src/base.h: Likewise.
	* src/crk.h: Likewise.
	* src/data.h: Likewise.
	* src/fileformat.h: Likewise.
	* src/main.cpp: Likewise.
	* src/mol.h: Likewise.
	* src/oberror.cpp: Likewise.
	* src/oberror.h: Likewise.
	* src/obutil.h: Likewise.
	* src/residue.cpp: Likewise.
	* src/base.cpp [#include <fstream>]: Removed and...
	[HAVE_FSTREAM, HAVE_FSTREAM_H]: Replaced with this.
	* src/binary.h: Likewise.
	* src/crk.h: Likewise.
	* src/data.h: Likewise.
	* src/fileformat.h: Likewise.
	* src/main.cpp: Likewise.
	* src/mol.h: Likewise.
	* src/residue.cpp: Likewise.

2004-02-06  Geoff Hutchison  <ghutchis@wso.williams.edu>

	* Doxyfile, doc/Doxyfile-man: Updated for 1.100.2 release and for
	newer doxygen versions.

	* src/base.cpp: Doxygen doc updates. Should produce links to key
	classes in intro page.

	* src/base.h, src/generic.cpp, src/mol.cpp, src/mol.h,
	src/math/matrix3x3.cpp, src/math/vector3.cpp: Documentation updates.
	
	* src/data.cpp, src/data.h, src/extable.txt, src/extable.h,
	src/fileformat.cpp, src/fileformat.h: Updates for PQS format,
	contributed by Pawel Wolinski <pwolinsk@uark.edu>.

	* src/pqs.cpp: New code, as above.

	* src/typer.h: Turn on root-atom detection via SelectRootAtoms.
	
	* test/smartstest.txt, test/attype.00.smi: Add new tests for ring
	and SMARTS testing.

	* test/ringresults.txt, test/smartsresults.txt: Regenerated.

2004-02-06  Michael Banck  <mbanck@gmx.net>

	* src/cml.cpp (debug): Fix loop control variable j.
	Contributed by Francesco Bresciani <fbresciani@users.sourceforge.net>.

2004-02-04 Fabien Fontaine <ffontaine@imim.es>

	* tools/obgrep.cpp: -n option added to print 
	the name of the matched molecules only

2004-02-02  Geoff Hutchison  <ghutchis@wso.williams.edu>

	* src/mdl.cpp (WriteSDFile): Add support for wedge/hatch, as
	contributed by Eugen Leitl. 

	* src/typer.cpp: Remove debugging information used for
	SelectRootAtoms development. (oops!)

	* configure, aclocal.m4, */Makefile.am, */Makefile.in: Regenerated
	using autoconf 2.59, automake 1.8.2, libtool 1.5.2 to prevent recent
	problems linking on Mac OS X.

2004-01-31  Michael Banck  <mbanck@gmx.net>

	* src/crk.cpp (WriteCRK): Rename loop variable from n to m.
	Contributed by Chris Morley <c.morley@gaseq.co.uk>.

2004-01-27  Michael Banck  <mbanck@gmx.net>

	* configure.in: Rename defined directory DATADIR to BABEL_DATADIR,
	to avoid problems with other programs. Closes: #885828.
	* configure: Regenerated.
	* src/babelconfig.h.in: Regenerated.
	* src/bondtyper.cpp: Rename _dir from DATADIR to BABEL_DATADIR.
	* src/data.cpp: Likwise.
	* src/pdb.cpp: Likewise.
	* src/phmodel.cpp: Likewise.
	* src/rotor.cpp: Likewise.
	* src/typer.cpp: Likewise.

2004-01-16  Michael Banck  <mbanck@gmx.net>

	* test/Makefile.am (TESTS): Run test-scripts with absolute path
	$(top_srcdir)/test/, as they are not present in the current
	directory when building in a builddir.
	(SUBDIRS): New directive, including cmltest.
	* test/Makefile.in: Regenerated.
	* test/cmltest/Makefile.am: New file.
	* test/cmltest/Makefile.in: New file, autogenerated.
	* configure.in (AC_OUTPUT): Added test/cmltest/Makefile.
	* configure: Regenerated.

2004-01-16  Michael Banck  <mbanck@gmx.net>

	Second and last pass of test-suite-in-builddir handling, fixing
	the shell scripts.

	* test/cml.sh: Export $srcdir and $builddir for use in the called
	scripts. Prepend $cmltestdir to avoid absolute pathnames.
	* test/cmltest/test.sh: Define $builddir and specify path to the
	babel executable from there. Replace every absolute pathname
	with relative ones.
	* test/roundtrip.sh: Likewise.

2004-01-13  Michael Banck  <mbanck@gmx.net>

	First shot at making a test-suite run from a build-directory
	possible.

	* test/Makefile.am (AM_CPPFLAGS): Added -DTESTDATADIR, pointing
	to $(srcdir)/test, where the test-suite data files are located.
	* test/Makefile.in: Regenerated.
	* test/ringtest.cpp: Use strings to handle the filenames. If
	TESTDATADIR is defined, prepend it to the filename. Transform
	the string to a char* when calling SafeOpen().
	* test/smartstest.cpp: Likewise.
	* test/unitcell.cpp: Likewise.

2003-12-17  Fabien Fontaine <ffontaine@imim.es>
	* src/typer.cpp (AssignImplicitValence)
	Fix protonation of amines bug due to spin multiplicity assignment 
	move spin multiplicity lines of CM to OBMol::AssignSpinMultiplicity();
	Add spin multiplicity flag  and related functions in mol.h
	
2003-12-08  Fabien Fontaine <ffontaine@imim.es>
	* src/mol.cpp (DeleteNonPolarHydrogens): uncomment "IncrementMod();"
	to avoid _c initialization when Deleting atoms
	
2003-12-01  Michael Banck  <mbanck@gmx.net>

	* src/data.cpp (GetExtension): Make sure that the returned long
	description stays defined by copying it into a static string.
	Contributed by David Mathog <mathog@mendel.bio.caltech.edu>.
	* (GetDescription): Likewise.

2003-12-01  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/pdb.cpp, src/mol2.cpp, src/gromos96.cpp:
	Residue perception via OBChainsParser turned on, contributed by
	Malcolm Gillies.

	* src/atom.cpp(GetResidue): Call SetChainsPerceived--ensures that
	chain perception code will not be called repeatedly on failure,
	contributed by Malcolm Gillies.

	* src/crk.h, src/crk.cpp: New files for Chemical Resource Kit
	support, contriubted by Alex Clark <a.clark@bendcable.com>. 

	* src/Makefile.am, src/Makefile.in: As above.

	* src/data.cpp, src/data.h, src/extable.txt, src/fileformat.cpp,
	src/fileformat.h: Changes to support read/write CRK files.

	* src/mol.h: Fix mistake in HasAromaticBond() -- aromatic bonds
	are "5", not "4"!

	* src/bond.cpp(IsTriple): New method for triple bonds.
	
	* src/typer.h: Add new SelectRootAtoms method (defaults to off)
	contributed from JOELib.

	* src/typer.cpp(SelectRootAtoms): Implement it, as ported from
	Joerg Wegner's JOELib Java code.
	(AssignAromaticFlag): Use SelectRootAtoms -- if
	avoidInnerRingAtoms is set to false, same as previous code.
	
	* test/attype.00.smi: Add new test patterns from JOELib tests.

	* test/smartstest.txt: Add new pattern from JOELib tests.

	* test/ringresults.txt, test/smartsresults.txt: Regenerate using
	current code (i.e. do not use SelectRootAtoms code yet).

2003-11-22  Michael Banck  <mbanck@gmx.net>

	* src/povray.cpp (CalcBoundingBox): Decouple declaration of 
	`unsigned int i' from the for-loop. Contributed by Chris 
	Morley.
	
	* src/Makefile.am (libopenbabelinclude_HEADERS): Added 
	babelconfig.h

	* src/Makefile.am (EXTRA_DIST): Removed newmain.cpp, bondtyper.cpp
	and bondtyper.h for now.

	* configure.in (AC_CHECK_HEADERS): Added iostream and fstream. 
	They seem to be not always present on Irix.
	(AC_HEADER_TIME): New macro.
	(AC_CHECK_TYPES): New macro, checking for clock_t.
	(AC_CHECK_FUNCS): New macro, checking for rint.
	* configure: Regenerated.
	* src/babelconfig.h.in: Regenerated.
	* src/bitvec.h, src/grid.h, src/obifstream.h: Include 
	babelconfig.h.
	[__sgi]: Removed, and ...
	[HAVE_IOSTREAM]: Replaced with this.
	* src/main.cpp [__BORLANDC__]: Removed strncasecmp function.
	* src/mol.h [WIN32]: Removed, and ...
	[!HAVE_RINT]: Replaced with this.
	* src/obutil.h: Include babelconfig.h.
	[WIN32] (include time.h): Removed, and ...
	[TIME_WITH_SYS_TIME, HAVE_SYS_TIME_H]: Replaced with this.
	[WIN32] (clock_t): Removed, and ...
	[HAVE_CLOCK_T]: Replaced with this.
	* src/obutil.cpp: Include babelconfig.h.
	[WIN32]: Removed, and ...
	[HAVE_CONIO_H]: Replaced with this.
	
	* src/mol.h: Move #include for babelconfig.h to the top.

2003-11-21  Fabien Fontaine  <ffontaine@imim.es>

	* src/typer.cpp (AssignImplicitValence): SetAromaticPerceived flag on 
	to ensure that the aromatic typer is not called when assigning 
	implicit valence, and check if implicit valence has already been 
	assigned
	*src/mol.h SetFlags function added

2003-11-17  Michael Banck  <mbanck@gmx.net>

	* src/chemtool.cpp (WriteGX): Replaced round(x) with 
	floor(x + 0.5) as round() does not seem to be available on all
	systems.

2003-11-13  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* test/attype.00.smi, test/smartstest.txt: Updated with tests from
	SMARTS bug reports PR#776094 and 796649. (Not problems in the
	current CVS code, but added as new test cases anyway.)

	* test/ringtest.cpp, test/smartstest.cpp: Allow generation of the
	results files.

	* test/ringresults.txt, test/smartsresults.txt: Generated from
	above. (Current implementation validates previous results as well
	as new patterns above.)

	* src/atomtyp.txt: Add Sybyl typing for phenol oxygens and ester
	sp3 oxygens.

	* src/aromatic.txt: Remove phosphole from aromatic patterns, add
	selenophene. (Phosphole is a classic non-aromatic system due to
	the non-planar P atom.)

	* src/aromatic.h, src/atomtyp.h, src/extable.h: Autogenerated.

2003-11-12  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* configure, */Makefile.in: Regenerated using automake-1.7.8,
	autoconf-2.58, libtool-1.5. (Should solve problems with CVS
	versions needing to call autoconf/automake due to timestamps.)
	Added AM_MAINTAINER_MODE.

	* aclocal.m4: Put this back in CVS so aclocal isn't needed to
	generate from CVS.

	* src/Makefile.am: Add back targets for generating binary .h files
	from the appropriate .txt data files.
	
	* src/parsmi.cpp: Back out Chris Morley's SMILES extensions for
	radicial centers--current causes many SMILES/SMARTS errors.

	* test/main.cpp, test/roundtrip.cpp: Add comments on additional
	tests needed soon. (These are all in the project tracker for
	Release-2.0 roadmap too.)

	* test/smartstest.cpp: Improve debugging information on
	errors--give the first atom matching if we have a match, even if
	it's not the number of matches expected by the results.

	* test/test-set.sh: Only mention running the test if the test
	directory is there.

2003-11-11  Michael Banck  <mbanck@gmx.net>

	* configure.in (AC_CHECK_HEADERS): Removed checks for some C++ 
	Standard Library headers (algorithm vector map list iostream
	fstream deque). They are essential for openbabel and we don't
	even exit if ./configure does not find them, let alone work
	around them missing. So in effect, the checks were just wasting
	CPU cycles until now. stream and strstream are still checked.
	* configure: Regenerated.
	* src/babelconfig.h.in: Regenerated.

2003-11-11  Michael Banck  <mbanck@gmx.net>

	* configure.in (AC_REPLACE_FUNCS): New directive, including
	snprintf and strncasecmp for now.
	* configure: Regenerated.
	* src/babelconfig.h.in: Regenerated.
	* src/Makefile.am (libopenbabel_la_LIBADD): Added @LTLIBOBJS@.
	(babel_LDADD): Added @LIBOBJS@.
	* src/Makefile.in: Regenerated.
	* src/mol.h [!HAVE_SNPRINTF] (snprintf): New declaration.
	[!HAVE_STRNCASECMP] (strncasecmp): Likewise.
	* src/snprintf.c: New file, taken from Mark Martinec.
	* src/snprintf.h: New file, taken from Mark Martinec.
	* src/strncasecmp.c: New file, taken from gnulib.

2003-11-10  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/main.cpp: Silently ignore cases with 0 atoms in a
	molecule. This is frequently due to trailing empty lines in a file.

	* test/Makefile.am: Separate the tests in obtest into distinct
	programs and tests for a more accurate "make check" result.

	* test/main.cpp, test/matrixtest.cpp, test/ringtest.cpp,
	test/smartstest.cpp: Ditto.

	* test/cml.sh: Skip the test in "make check" if the directory
	doesn't exist.

	* test/test-set.sh: If the test-set directory exists, then run it,
	otherwise skip the test.

2003-11-10  Michael Banck  <mbanck@gmx.net>

	* src/math/Makefile.am (AM_CPPFLAGS): New directive, adding
	$(top_srcdir)/src to includes to make building in a seperate
	build-directory possible.
	* test/Makefile.am (AM_CPPFLAGS): Likewise.
	* tools/Makefile.am (AM_CPPFLAGS): Likewise.
	
2003-11-08  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* configure, */Makefile.in: Updated with current autoconf, automake.

	* src/cml.cpp: Add OB_VERSION const to use BABEL_VERSION
	definition (automagically updated from babelconfig.h).

	* src/extable.txt: Add .ins ShelX extension.

	* test/cml.sh, test/cmltest/roundtrip.sh, test/cmltest/test.sh:
	Shell scripts to run the CML test suite.
	
	* test/smartstest.cpp: Add more debugging information if SMARTS
	tests fail.

	* src/shelx.cpp: Changes suggested by Louis Ricard
	<louis.ricard@polytechnique.fr> to fix several ShelX bugs.
	
	* src/atom.cpp, src/mdl.cpp, src/mol.h, src/parsmi.cpp, src/smi.cpp,
	src/typer.cpp: Support for radical centers, contributed by Chris
	Morley <c.morley@gaseq.co.uk>.

	Changes contributed by Malcolm Gillies <malcolm.b.gillies@anu.edu.au>

	* src/pdb.cpp: (ReadPDB) BUG: OBMol::ConnectTheDots() fails inside a
	OBMol::BeginModify()/OBMol::EndModify() pair, as OBMol::Has3D() will
	always return false. Rearrange to avoid this.
	Use ATOM/HETATM record element name information in columns
	76-77 to determine element type when possible and necessary.
	(WritePDB) BUG: writing of ATOM records could fail due to use of
	strncpy to copy atom and residue names (i.e. lack of null termination)
	Write HETATM records instead of ATOM records if
	res->IsHetAtom(atom)==true

2003-10-22  Fabien Fontaine <ffontaine@imim.es>

	Added tools directory: contains utilities using libopenbabel

	* tools/obgrep.cpp: SMART molecule grep
	* tools/obfit.cpp: align molecules on SMART substructures
	* tools/obrotate.cpp: rotate a tortional bond matching a SMART

	* src/mol2.cpp (ReadMol2): Continue reading untill EOF or next
	MOLECULE record to avoid 'error: has zero atoms!' message

2003-10-20  Fabien Fontaine <ffontaine@imim.es>

	* src/patty.h Istype function added  
	* src/patty.cpp (read_rules) (assign_rules) : memory 
	allocation changed
	
2003-10-03  Michael Banck  <mbanck@gmx.net>

	* src/atom.cpp (GetDistance): New method.
	* (GetAngle): Likewise.
	* src/mol.h: Declare them.

2003-09-25  Michael Banck  <mbanck@gmx.net>
	
	* src/mol.h (_internals): New protected class attribute.
	(GetInternalCoord): New declaration.
	(SetInternalCoord): New method.
	(BeginInternalCoord): New method for iterating over _internals.
	(NextInternalCoord): Likewise.

	* src/mol.cpp (GetInternalCoord): New method.

2003-09-21  Peter Murray-Rust <pm286@cam.ac.uk>

	* src/cml.cpp: removed inconsistences between atomRef, atomRef1,
	atomRef2, and atomRefs2 attributes.
	* test/cmltest/cs2a.cml: changed atomRefs1 and atomRefs2 attributes 
	to atomRef
	* test/cmltest/test.bat: minor edits to allow echo
	
2003-09-14  Michael Banck <mbanck@gmx.net>

	* src/povray.cpp: Removed getlogin() declaration.
	(OutputHeader): Removed log_name, don't output login name in
	Povray header.

2003-09-13  Michael Banck <mbanck@gmx.net>

	* src/chemtool.cpp: New file.
	* src/Makefile.am (libopenbabel_la_SOURCES): Added chemtool.cpp.
	* src/Makefile.in: Regenerated.
	* src/data.cpp (TextToType): Added CHEMTOOL typestring.
	* src/data.h (io_type): Added CHEMTOOL.
	* src/extable.txt: Added Chemtool writing support.
	* src/fileformat.cpp (OBFileFormat::WriteMolecule): Added
	WriteCHT function call for CHEMTOOL io_type.
	* src/fileformat.h (OBFileFormat): Added WriteCHT() method.

2003-07-02  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/data.cpp (GetAtomicNum): Fix mistake in calls to strcasecmp,
	as pointed out by Chris Morley.

	* src/math/Makefile.am: Make sure math/matrix3x3.h and
	math/vector3.h are installed with a "make install" as reported by
	Jean. Fixes PR#761644.

	* Makefile.in, src/math/Makefile.in: Regenerate using automake.

2003-06-26  Michael Banck <mbanck@gmx.net>

	* NEWS: New file.

2003-06-23  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/main.cpp: Move man documentation to end of file (right
	after usage() code) for better "readability". Also make sure that
	we quit attempting to read molecules from a file when we find no
	atoms. Fixes some problems seen when roundtripping.
	
2003-06-23  Vincent Favre-Nicolin  <vincefn@users.sourceforge.net>

	* doc/Doxyfile-man, src/main.cpp, Doxyfile: Added the man page in
	doxygen format in main.cpp, for automatic generation of the man
	page from the source code, using 'Doxygen Doxyfile-man' in the
	doc/ directory. Normal Doxygen run now skips src/main.cpp

2003-06-23  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* doc/babel.1: Update with current command-line options and file
	formats.

	* doc/FAQ.html, doc/Migration.html: Use HTML versions from the
	website rather than text versions.

	* doc/Makefile.in, doc/Makefile.am: Update accordingly.

	* src/cml.cpp(CleanUp): New method to clean up global vector
	variables after ReadCML() and WriteCML() in attempt to fix
	PR#736001. Does not change normal command-line path, but should
	fix library useage when multiple files may be loaded.
	
2003-06-22  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/main.cpp: Updated to use OBMolVector for importing multiple
	molecules in a file. Emulates behavior of babel, but improvements
	are needed. (e.g. some formats like separators if it's a
	multi-molecule file.) Also added -f (first) and -l (last)
	flags to pick particular molecules. Fixes feature request
	PR#575243.

	* src/molvector.h, src/molvector.cpp: Add new PushMol method for
	adding new molecules to the end of the vector. Improve Write
	support and clean up overloaded methods--instead use default
	values.

	* src/*.txt: Updated copyright information.
	* src/*.cpp, src/*.h: Ditto. (As needed)

2003-06-15  Michael Banck  <mbanck@gmx.net>

	* src/pdb.cpp (WritePDB): Cast atom->GetValence() to int to avoid
	signed/unsigned comparison.
	* src/ring.cpp (OBRingSearch::RemoveRedundant): Remove register
	attribute for i and j to suppress a warning.
	(OBRingSearch::AddRingFromClosure): Remove unused parameter 
	`int level'.
	* src/ring.h (AddRingFromClosure): Remove unused parameter `int'.  
	* src/povray.cpp (OutputMoleculeBonds): Remove unused parameter
	`OpenBabel::OBMol mol'.

2003-06-15  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	Updated to latest autoconf, libtool, automake versions to fix a
	few minor bugs and add libtool-1.5 Mac OS X support.

	* depcomp, ltmain.sh, missing, config.guess, config.sub: Update
	with latest GNU versions.
	
	* src/Makefile.am: Remove main.cpp from targets for library and
	make sure only one set of headers exists.

	* configure: Regenerated with autoreconf.
	* */Makefile.in: Regenerated with autoreconf.
	
	* openbabel.pc.in: Fix minor typo.

	* *.cvsignore: Add new libtool-generated files.

2003-06-15  Michael Banck  <mbanck@gmx.net>

	Added pkg-config support, thanks to Jean Bréfort.

	* openbabel.pc.in: New file.
	* configure.in: Add openbabel.pc to AC_OUTPUT.
	* configure: Regenerated.

	* Makefile.am: Install openbabel.pc in $(libdir)/pkgconfig.
	* Makefile.in, doc/Makefile.in, src/Makefile.in, src/math/Makefile.in,
	src/windows/Makefile.in, test/Makefile.in: Regenerated.

	* ChangeLog: Fixed a wrong filename in the last entry, minor
	corrections.

2003-06-15  Michael Banck  <mbanck@gmx.net>

	Switch to automake.
	
	* configure.in: Add automake macros, remove some obsolete legacy 
	stuff (AC_PATH_PROG(AR), $top_builddir, AC_PROG_RANLIB), enable
	libtool.
	* configure: Regenerated.

	* Makefile.am, test/Makefile.am, doc/Makefile.am, 
	src/Makefile.am, src/windows/Makefile.am, src/math/Makefile.am:
	New files.
	* doc/Makefile.in, src/windows/Makefile.in: New files, 
	autogenerated.
	* Makefile.in, test/Makefile.in, src/Makefile.in, 
	src/math/Makefile.in: Regenerated.

	* src/babelconfig.h.in: Regenerated.

2003-06-10  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/Makefile.in: Remove cwrap.h and cwrap.o as these have been
	removed.

	* src/mol.cpp (ConnectTheDots, PerceiveBondOrders): Make sure
	these don't try to run on 2D structures, fixing bug reported by Jean.

	Remove a variety of compilation warnings.
	
	* src/cacao.cpp, src/data.cpp, src/mpqc.cpp, src/povray.cpp:
	Remove unused variables.

	* src/mopac.cpp: Switch to unsigned int to remove comparison warning.

	* src/obutil.h, src/obutil.cpp: New IsNear and IsNearZero
	functions to compare floating-point numbers instead of unreliable
	== operator. Ensures numbers are closer than a user-specified epsilon.
	
	* src/mol.cpp, src/mol2.cpp, src/molchrg.cpp,
	src/math/vector3.cpp: Improve reliability of floating point
	comparisons. Replace == operator compares with new IsNear and
	IsNearZero functions.
	(vector3.cpp probably should replace == operators as well.)

2003-06-09  Michael Banck  <mbanck@gmx.net>

	* src/main.cpp: 
	- Renamed usage() to help() and added a new method
	  usage() with a terse usage output and a pointer to help().
	- Added -H switch for help().
	- Reformatted usage() and help() output a bit.
	- Print out an 'unrecognized option' error message when 
	  option parsing fails, in addition to usage().
	
2003-06-09  Michael Banck  <mbanck@gmx.net>

	* src/main.cpp: Added a pointer to the Mailing-List for sending
	bug-reports to the end of usage().

2003-06-03  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/mol.h: Add OBBond::UnsetKekule for PerceiveBondOrders support.

	* src/mol.cpp (PerceiveBondOrders): Add aromatization pass to
	complete the implementation of Roger Sayle's algorithm. This
	implementation is rather naive--it types rings, then passes back
	to Kekulize. Still works pretty well.

2003-06-02  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	Lots of minor bugfixes from the past few weeks.
	
	* src/element.txt, src/element.h: Update phosphorus covalent
	radius as suggested by Jean Brefort (PR# 743730).

	* configure.in: Set "bug report" e-mail address, remove custom
	openbabel_version variable and change AC_CHECK_HEADER calls to
	AC_CHECK_HEADERS to gracefully fail and set appropriate HAVE_
	macros in babelconfig.h.

	* configure, src/babelconfig.h.in: Updated using autoreconf.

	* src/data.cpp, src/povray.cpp: Update to use sstream if present
	and strstream if not.

	* src/main.cpp: Add support for centering coordinates as well as
	checking for 2D vs. 3D coordinates when calling WriteMolecule.

	* src/mol.cpp (Kekulize): Remove 255 atom limit discovered by
	Frank Schmitt. (Some limit is reasonable, but SMARTS matching
	should be better about matching non-Kekulize'd structures.)
	
	* Makefile.in: Use standard PACKAGE_VERSION instead of
	openbabel_version.

	* src/Makefile.in: Get rid of bitgrid.[cpp,h] and version.h.

	* src/base.cpp: Minor documentation updates.
	
	* src/extable.txt, src/extable.h, src/fileformat.cpp,
	src/fileformat.h: Add in preliminary hooks for ShelX support.

	* src/shelx.cpp, src/bondtyper.h, src/bondtyper.cpp: New files.

	* test/unitcell.cpp: Make sure to initialize variables.

	* test/roundtrip.cpp: Use babelconfig.h instead of version.h.
	
2003-04-22  Michael Banck  <mbanck@gmx.net>

	Switch to autoheader.

	* configure.in: 
	- Added AC_CONFIG_HEADER(src/babelconfig.h), along with
		  two cosmetic changes.
	- Define BABEL_VERSION and DATADIR directly. DATADIR is
	  $datadir/openbabel (i.e. $pkgdatadir in automake) as before,
	  $datadir can be specified normally as option to ./configure.
	* configure: Regenerated

	* acinclude.m4: New file for local m4-files, it currently has
	  AC_DEFINE_DIR (used in order to define DATADIR) as content.

	* src/Makefile.in: Removed DEFS=-DDATADIR.
	* src/math/Makefile.in: Likewise.

	* src/babelconfig.h.in: Initial autogenerated version

	* src/bgf.cpp, src/gromos96.cpp, src/newmain.cpp, src/patty.cpp,
	  src/pdb.cpp, src/phmodel.cpp, src/rotor.cpp, src/typer.cpp: 
	  Added #include "babelconfig.h", removed #include "version.h" where
	  applicable

2003-04-13  Michael Banck  <mbanck@gmx.net>

	* config.guess, config.sub, ltmain.sh: Updated to current FSF
	versions

2003-03-26  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/bondtyp.txt, src/bondtyp.h: Files for the PerceiveBondOrder
	functional group typing.

	* Src/Makefile.in: Add rules for generating bondtyp.h from the
	text file.
	
2003-03-18  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/configure.in, src/configure: Update with "1.100.1" version
	number for now.

	* test/roundtrip.cpp (main): Improve usability by making BOX and
	SMI file types special cases. (Hard to compare atoms and coords
	for them.)

	* src/pdb.cpp (WritePDB): Make sure two-character elements end up
	in the correct columns even with residue information.
	(ReadPDB): Refer to the original text (with spaces) before
	deciding if an atom type is a two-character element or a one-char
	element.

	* src/*.h, *.cpp: Switch from float to double throughout. Improves
	accuracy significantly in coordinate transformations.

	* src/binary.cpp (PackCoordinate, UnpackCoordinate): Make sure to
	cast to float before packing for compatibility with file format.

2003-03-17  Geoff Hutchison  <hutchisn@chem.northwestern.edu>
	
	* src/cml.cpp: Fix unused declarations that cause problems with
	the xlC compiler reported by Markku Laukkanen.

	* src/pdb.cpp (WritePDB): Make sure strings are null-terminated
	and convert element name to upper-case for the column 77-78 entry.

2003-03-15  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/cacao.cpp (SetHilderbrandt): Fix one coredump found on
	roundtrip testing and add notes on what needs to be fixed for an
	iron-clad version of this.

	* src/hin.cpp (ReadHIN): Fix bugs uncovered in latest roundtrip
	testing. Some lines are exactly 11 tokens long and this should not
	stop the parsing.

	* src/viewmol.cpp (ReadViewMol): Remove loop attempting to find
	the first part of the file -- causes problems since two lines
	instead of one are read off the top.

2003-03-14  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* test/roundtrip.cpp: Detect cases where files have no atoms.

	* test/Makefile.in: Add targets to create roundtrip program.

	* src/cacao.cpp (ReadCaccrt): Add an obUnitCell object to the
	molecule for keeping track of unit cell data. Call
	ConnectTheDots() and PerceiveBondOrders().
	(WriteCaccrt): If obUnitCell data exists, output unit cell data.

	* src/dmol.cpp (ReadDMol): Read $cell vector section if it exists.
	(WriteDMol): Write unit cell data if it exists.

	* src/parsmi.cpp (ParseComplex): Add isotope information where
	available.
	(ParseSmiles): Fix segfault in test file caused
	by white-space character ending SMILES string. Ignore them and
	continue.
	
	* src/generic.cpp (OBUnitCell::ctor): Make sure to set the generic
	data type and string.

	* src/generic.h (OBUnitCell): Add space group storage as a string.

	* src/mopac.cpp (ReadMOPAC): Remove debugging output left
	inadvertantly.

2003-03-13  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/generic.cpp, src/generic.h: Add support for retrieving
	OBUnitCell data by a matrix3x3 and add a new method
	GetOrthoMatrix() for translating from fractional coordinates.

	* src/pdb.cpp (WritePDB): Make sure that residue and atom names
	are only a maximum length to format columns correctly. Add element
	symbol to columns 77 and 78 to better support the current standard.

	* src/math/matrix3x3.h, src/math/matrix3x3.cpp (SetRow): New
	method to set a row from a vector3.

	* src/mpqc.cpp (ReadMPQC): Fix bugs in MPQC output reading,
	uncovered by new roundtrip testing.

	* src/cacao.cpp (WriteCaccrt): Make sure the buffer is actually
	written to the file. Previous output passed roundtrip testing (no
	coredumps) but nothing was actually written to the file!

	* test/roundtrip.cpp: New file to compare two molecules after
	roundtripping. (Perhaps it should be called molcmp?)

2003-03-12  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/pdb.cpp (WritePDB): Make sure that atom type definitions are
	written correctly as per the PDB standard. Fixes PR#701620, as
	suggested by Christian Hofbauer.

	* src/report.cpp: Add support for listing torsion angles. Angle
	and torsion data is not the same as in Babel 1.6, but is
	coherent. Fixes feature request PR#568088.

2003-03-11  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/generic.h, src/generic.cpp: Add new OBUnitCell generic data
	type for storing and converting unit cell parameters.

	* test/Makefile.in, test/main.cpp, test/unitcell.cpp,
	test/unitcell.txt: New test to confirm that unit cell conversions
	occur accurately.

	* src/atom.cpp (GetAtomicMass): If an isotope has been set, use
	that as the atomic mass. Otherwise stick to the normal average
	mass.

	* src/isotope-small.txt, src/Makefile.in, src/isotope.h: Smaller
	set of isotopes for generating default data to minimize growth in
	executable size.

2003-03-10  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/generic.h: Add new data type "slots" to obDataType enum.

	* src/pdb.cc (class OBSerialNums): Use new obSerialNums enum,
	rather than overriding user-defined obData# slots. These should
	only be used by user-level code, not by Open Babel internal code.

2003-03-10  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/atom.cpp (GetType): Fix compilation error when copying in
	"D" atom type for a deuterium atom.

	* src/mol2.cpp (ReadMol2): Read residue information if available,
	along with atom IDs and residue number. Should address feature
	request PR#697614 submitted by Claudio Cavasotto.

	* src/mol.h, src/mol.cpp: Add API support for a molecular total
	charge.
	
2003-03-06  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/Makefile.in: Add support for the isotopes.txt -> .h
	conversion.

	* src/mol.h: Add support for storing atomic isotopes and
	"exact mass" methods.

	* src/mol.cpp (GetExactMass): New method, using isotope table.
	
	* src/atom.cpp: Add support for atomic isotopes, including
	modifications to SetType methods to ensure deuterium isotopes
	translate into the correct type.

	* src/data.h, src/data.cpp: Add new OBIsotopeTable for obtaining
	exact masses of given element/isotope combinations.
	(OBElementTable::GetAtomicNum): Expand to return isotope for the
	case where a "D" or "T" symbol is passed to the table.

	* src/isotope.txt, src/isotope.h: Update to include atomic numbers
	at the beginning of each line so ParseLine() knows "where it is."

	* src/report.cpp: Add output of molecular weight and exact mass.
	
2003-03-02  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/mopac.cpp (ReadMOPAC): Fix problems with MOPAC partial
	charges not translating to other formats. (PR#660364)

2003-02-28  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/mol.h: Add new method OBAtom::SamllestBondAngle().

	* src/atom.cpp: Implement it.

	* src/mol.cpp (ConnectTheDots): Improve bonding by checking for
	small bond angles while removing long bonds.

	* src/pdb.cpp: Do not attempt to add alpha-peptide
	bonds. ConnectTheDots is much better at this and doesn't attempt
	to add strange very long-range bonds.

2003-01-17  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/base.cpp: Add introduction to documentation.

	* src/mol.cpp: Don't attempt to use Kekulize or Begin/EndModify to
	detect aromaticity.

	* src/phmodel.txt, src/phmodeldata.h: Add histidine transformation
	and unused tryptophan transformation to help with problems
	observed by Richard Gillian (PR#575964).

	* src/resdata.txt src/resdata.h: Fix incorrect bond orders for ARG
	residue. (PR#575964).

	* test/matrixtest.cpp: Make matrix test a bit more verbose to make
	debugging problems on non-Intel platforms easier. (Relates to
	PR#647417). Problem derives from numerical accuracy, but loosening
	tolerance to 2e-6 solves the problem.
	
2003-01-15/16  Peter Murray-Rust <pm286@cam.ac.uk>
	       Vincent Favre-Nicolin <vincefn@users.sourceforge.net>
	
	* Changes to compile using the free Borland C++ compiler 5.5.1:
	a few fixes ISO C++ fixes, some #if's and added makefiles (bc32.mak)

2003-01-06  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	Lots of changes from the last few weeks.

	* Doxyfile: Bump to 1.100.1. Make sure to check .cpp files for
	docs too.

	* Makefile.in (disclean): Remove doc/API for distclean target.

	* src/bin2hex.pl: Add comments to headers about contents
	(i.e. what the headers do, why they're not supposed to be installed).

	* src/chains.cpp: Fix mistake in Pyroglutamate comment from OE->OB
	find/replace.
	
	* src/aromatic.txt, src/atomtype.txt, src/element.txt,
	src/extable.txt, src/phmodel.txt, src/resdata.txt, src/types.txt:
	Updated comments at top of file to give more documentation and
	point to appropriate source files.
	
	* src/atom.cpp, src/base.h, src/bitvec.cpp, src/bitvec.h,
	src/bond.cpp, src/data.cpp, src/data.h, src/fileformat.h,
	src/generic.h, src/grid.h, src/mol.cpp, src/mol.cpp,
	src/molchrg.cpp, src/molchrg.h src/obutil.cpp, src/obutil.h,
	src/parsmart.cpp, src/parsmart.h, src/patty.cpp, src/patty.h,
	src/residue.cpp, src/ring.cpp, src/ring.h, src/rotor.h,
	src/typer.cpp, src/typer.h, src/math/matrix3x3.cpp,
	src/math/matrix3x3.h, src/math/vector3.cpp, src/math/vector3.h: 
	Updated Doxygen docs (brief comments in headers mostly, more
	detailed documentation if available in source files.)

	* src/mol.h, src/atom.cpp, src/bond.cpp, src/residue.cpp: Added
	support for OBGenericData stored in an OBAtom, OBBond or
	OBResidue.

	* src/isotope.txt, src/bond.txt: New data files, currently unused.

	* src/Makefile.in: Add new zindo target.

	* src/zindo.cpp: New ZINDO input support.

	* src/fileformat.cpp: Call WriteZindo for ZINDO file format.

2002-12-19  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/tinker.cpp: Add comment that we're using MM2 parameters.

	* src/mol.cpp (PerceiveBondOrders): Fix ugly bug causing crashes
	with -NO2 groups, reported by Richard Muller.

2002-12-6  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* Doxyfile, README, configure.in, configure: Rename version to
	1.100.

	* doc/FAQ: Update with a few minor corrections -- synched with
	website.

	* THANKS: A few omissions -- synched with website.

	* src/mdl.cpp, src/smi.cpp: Throw errors if the molecules are too
	large, rather than segfault. In MDL, this seems to be a problem
	inherent in the file format. For SMI, this seems to be a recursion
	problem.

	* src/mol.cpp: Take out the "quick-and-dirty" aromatic pass, it's
	too rough to be useful. I need to figure out how to pass off the
	rough structure to the Kekulize routines.

	* src/types.txt, src/types.h: Make sure there's a default type for
	an "N" atom, caught by roundtrip testing from PDB files. Fix typo
	for O.co2 v. Oco2 external type. Also add XYZ forms for all atoms.

2002-12-5  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* Doxyfile: Update for 1.99.1 release.

	* Makefile.in (dist): Make sure to remove autom4te.cache, etc.

	* THANKS: Plenty of contributors to this next release!

	* configure.in: When searching for the "ar" program, remember that
	Solaris has it in /usr/ccs/bin, so we should check there and fall
	back to the user's $PATH if it's not found.

	* configure: Regeneratre using autoconf-2.56.

	* src/mol.h: Add OBAtom::GetAtomicMass(), which will eventually
	allow for isotope-labeled atoms.

	* src/atom.cpp: Implement it.

	* src/mol.cpp (GetMolWt): Write using above.

	* src/cml.cpp: For the moment, turn off debugging output. (This
	should really go through a filter to allow users to see it if
	needed.)

	* src/data.cpp (GetAtomicNum): Have D and T atoms return for
	hydrogen. Eventually needs to support isotopes!

	* src/element.txt, src/element.h: Update with current IUPAC
	masses.

	* src/extable.txt, src/extable.h: Make sure all extensions are
	lowercased--since all comparisons are done this way. Fixes bug
	with GROMOS96 translation.

	* src/fileformat.h, src/fileformat.cpp (ReadMolecule): Fix thinko
	with BGF format--call the proper procedue.
	(WriteMolecule): Write Box format with a default size. (!!)

	* src/bgf (ReadBGF): Cleanups found from roundtrip testing. Make
	sure bond records are parsed correctly and run atom types through
	the element table rather than the type table.
	(WriteBGF): Change remark to "Open Babel" rather than "Babel" so
	people know who created the bugs. ;-)

	* src/gromos96.cpp: Use "Open Babel" as above.

	* src/pdb.cpp (WriteDelphiPDB, WritePDB): Make sure we have a
	proper line length. Bug discovered via roundtrip testing (and our
	complaints about short CONECT lines).

	* src/povray.cpp: Don't attempt to get a hostname--this introduces
	all sorts of portability problems. Also use the atomic symbol
	rather than the atom type. Less flexible, but POV-Ray doesn't like
	some symbols in atom types (e.g. period in O.co2).

	* src/report.cpp: Add report on chiral atoms. (Probably needs some
	improvement still.)

	* src/types.h, src/types.txt: Add missing types, notably for Ar
	and B2, uncovered via roundtrip testing to the c3d file formats.

	* src/xed.cpp: Fix segfault caused by improper iteration over
	bonds. Uncovered via roundtrip testing.

2002-12-3  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/data.cpp (ParseLine): Fix bug that ignored the first
	non-comment line.

2002-12-2  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* Makefile.in, src/Makefile.in, src/math/Makefile.in,
	test/Makefile.in: Makefile changes to ensure the src/math
	directory is called for all targets (including "clean") and
	testing data files are only copied when needed for $(builddir) !=
	$(srcdir).

	* src/mol.h: Fix typo in OBMol documentation.

	* src/mol.cpp: Add quick-and-dirty aromatic marking.

2002-12-1  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/bond.txt: Data file for functional group recognition in bond
	order typing (currently unused).

	* src/*.txt: Update with headers, including copyright notice.

	* src/data.cpp: Ignore comment characters (including copyright
	notice).

	* src/pdb.cpp: Fix bug with reading in residue data--the source
	file is called resdata.txt, not residue.txt.

2002-11-28  Vincent Favre-Nicolin <vincefn@users.sf.net>

	* src/obutil.cpp(CartesianToInternal): Fix the conversion
	from cartesian to Internal coordinates. The previous code
	prevented atom 3 to be bonded with atom 1, and any
	atom with idx>3 to be bonded with atom 1 or 2, which does
	happen sometimes...

2002-11-27  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/atom.cpp(MemberOfRingSize): New method -- if atom is in a
	ring, return size. (Need to fix to designate *which* ring if in
	multiple.)
	(SetHybAndGeom): Fix to use scaled bondlengths for aromatic,
	double and triple bonds and adjust expected bond angles if in a
	ring.

	* src/mol.cpp: Add << and >> convenience functions based on
	OBFileFormat objects.

	* src/ghemical.cpp: Fix header line for new "gpr" extension.

	* src/povray.cpp: New file, contributed by Steffen Reith to output
	POV-Ray scene files.
	
	* src/data.cpp, src/data.h: Add POV file type.

	* src/extable.txt, src/extable.h: Likewise.

	* src/fileformat.cpp, src/fileformat.h: Likewise:

	* src/Makefile.in: Add povray.cpp.

	* doc/README.dioxin.pov, doc/README.povray, doc/babel31.inc,
	doc/dioxin.mol2, doc/dioxin.pov: New files.

2002-11-25  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/mol.cpp, src/mol.h: Add back << and >> convenience functions
	derived from OBFileFormat methods.
	(PerceiveBondOrders): Fix bug introduced with previous change--do
	not add multiple bonds to an atom already with them. (Causes a
	problem with allene, but this is a rarer situation than the
	reverse.)

	* src/pdb.cpp (WritePDB, WriteDelphiPDB): Fix ATOM records with
	correct specification. Should fix problems with 5-digit atom
	numbers, PR#633719, reported by Kristian Rother.
	
2002-11-20  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* configure: Update with new autoconf-2.56

	* src/Makefile.in, src/math/Makefile.in, test/Makefile.in: Remove
	redundant '/' in path for $(builddir) variable.

	* src/data.cpp: Remove references to depreciated strstream.h
	header.

	* src/data.h: Likewise.

	* src/mol.cpp: Add "flat ring" pass for 5-member and 6-member
	rings. Ensures that these are all marked as potential sp2
	hybrids. Also prevent adding 4 bonds to nitrogen, for now.

	* src/extable.txt, src/extable.h: Update with new Ghemical
	extension .gpr (same format though).

2002-10-11  Michael Banck <mbanck@gmx.net>

	* src/math/Makefile.in: New file.

	* configure.in, src/Makefile.in: Adjust accordingly.

2002-09-19  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/aromatic.h, src/aromatic.txt: Remove selenium aromatic
	pattern for now. Fixes test/smartstest error on
	LABOTEST_LT-2_FE_546.

	* src/c3d.cpp, src/chains.cpp, src/cml.cpp, src/data.cpp,
	src/jaguar.cpp, src/main.cpp, src/mol.cpp, src/oberror.h,
	src/pdb.cpp, src/phmodel.cpp, src/ring.cpp, src/smi.cpp,
	src/xyz.cpp, src/math/matrix3x3.cpp, test/matrixtest.cpp,
	test/ringtest.cpp, test/smartstest.cpp: Warning fixes contributed
	by Erik Kruus.

2002-09-16  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/report.cpp: Fix bug noted by Markku Laukkanen, where the min
	variable would roll past the maximum number of atoms in a
	molecule.

	* src/cml.cpp (debugMolecule): Add default dimension. Contributed
	by Francesco Bresciani.

	* src/obutil.cpp (InternalToCartesian): If vector input is empty,
	quit early. Contributed by Francesco Bresciani.

	* src/oberror.h, src/oberror.cpp: Namespace fixes to handle
	problems on gcc-3.1 as reported by Austin Acton and Andrea Tasso.

2002-09-07  Michael Banck <mbanck@gmx.net>

	* Replaced $DEST with $DESTDIR and $bindir

2002-08-30  Stefan Kebekus <kebekus@users.sourceforge.net>

	* src/math/matrix3x3.cpp, src/math/matrix3x3.h added methods, in
	particular for eigenvector calculations. Modified the API slightly. 
	* src/mol.cpp, src/obutil.cpp removed a couple of undocumented
	mathematical methods that were never used. Make mol use the 
	eigenvector functions from the matrix class now.	
	* test/matrixtest.cpp test suite for matrix algebra and eigenvector 
	calculations
	* src/oberror.cpp, src/oberror.h introduced a generic error class

2002-08-27  Stefan Kebekus <kebekus@users.sourceforge.net>

	* test/main.cpp, test/ringtest.cpp, test/smartstest.cpp resolved
	namespace problems for the tests. The test compile now, but show
	errors.

2002-08-25  Michael Banck <mbanck@gmx.net>

	* Changed charge to int (gaussian.cpp) and divisor/exponent to 
	double (c3d.cpp) in order to compile on g++-3.2

2002-08-19  Stefan Kebekus <kebekus@users.sourceforge.net>

	* src/xyz.cpp more bugfixing. The translator behaves now smarter
	when interpreting the title and the atom types.

2002-08-16  Stefan Kebekus <kebekus@users.sourceforge.net>

	* src/pdb.cpp re-wrote routine to read CONECT records. That
	(hopefully) fixes a number of bugs in the earlier version
	(inability to read 5-digit serial numbers, accidently reading salt
	bridges as connects, etc.)

2002-08-15  Geoff Hutchison <hutchisn@chem.northwestern.edu>

	* src/pdb.cc: Fix ParseConectRecord to make sure we don't read
	beyond the vs.size().


2002-08-14  Stefan Kebekus <kebekus@users.sourceforge.net>

	* src/xyz.cpp Set molecule title to the string that is read from
	the second line of the file. Set atom type strings to strings read
	from the file. The reading method now prints +extensive+ and
	+detailed+ warnings to cerr if is xyz file is bad.

2002-07-26  Stefan Kebekus <kebekus@users.sourceforge.net>

	* src/math/matrix3x3.cpp added a few method useful in symmetry
	detection

2002-07-30  Michael Banck <mbanck@gmx.net>

	* src/gaussian.cpp: Last atom was written twice

2002-07-26  Stefan Kebekus <kebekus@users.sourceforge.net>

	* added documentation

2002-07-24  Michael Banck <mbanck@gmx.net>

	* src/gaussian.cpp: Added some comments/warnings; changed multiplicity
	calculation to use abs(charge) instead of charge
	* src/Makefile.in: splitted CXXFLAGS into CXXFLAGS, DEFS and INCS

2002-07-22  Michael Banck <mbanck@gmx.net>

	* src/gaussian.cpp: Fixed WriteGaussianCart() hopefully:
	- calculate charge and mulitiplicty (the latter being quite shakey
	  ATM)
	- do away with the zero after the atom type
	- don't reference coordinates, write them straight after the atom type
	- add a blank line at the end
	* configure.in: Added --enable-doxygen configure option. If enabled,
	search for doxygen and print an URL if not found
	* configure: Updated for new configure.in

2002-07-18  Stefan Kebekus <kebekus@users.sourceforge.net>

	* isolated the 3-dimensional vector and matrix classes, gave them
	less ambigous names, and put them into a directory "math". Added
	an access function to the vector class.

2002-07-12  Stefan Kebekus <kebekus@users.sourceforge.net>
	
	* added documentation to the Vector class, made some functions
	inline, and commented out some methods of the vector class
	whose names I found misleading, and which were never used anyway.

	* src/Vector.cpp, src/Vector.h: added documentation to the Vector
	class, made some functions inline, and commented out some methods
	of the vector class whose names I found misleading, and which were
	never used anyway.
	
	* configure.in: configure now warns if doxygen is not present

	* configure: Updated with autoconf 2.52d.

2002-07-09  Geoff Hutchison <hutchisn@chem.northwestern.edu>

	* src/fileformat.h, src/mol.h, src/binary.h: Fix missing std:: in
	headers, which causes problems with namespace-compliant
	compilers. Should fix PR #578522.

	* src/residue.cpp: Add missing "using namespace std," which also
	caused problems with namespace-compliant compilers.

	* src/pdb.cpp: For now, disable bond order perception--prevents
	problems with PDB residues observed by Richard Gillian.

2002-07-01  Richard Gillian <reg8@users.sourceforge.net>

	* src/main.cpp: Add -hpH to use the "adjust pH" setting for adding
	hydrogens.

	* src/mol.cpp: When adding hydrogens, if the parent atom is part
	of a residue, copy this information to the new atoms.

	* src/mol2.cpp: Write residue types and original atom names for
	the atoms if they exist.

2002-07-01  Stefan Kebekus <kebekus@users.sourceforge.net>

	* Doxyfile: New file, configuration for Doxygen API documentation.

	* configure: Updated with autoconf 2.52d.

	* Makefile.in: Added automatic API documentation generation. Use
	"make apidoc".

2002-06-29  Vincent Favre-Nicolin <vincefn@users.sourceforge.net>

	* src/mol.h: Give sample documentation from the primer for
	OBResidue, OBAtom and OBMol.

	* src/ring.h: As above for the OBRing class.

2002-06-28  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* configure.in, configure, src/Makefile.in: Adjust for various
	removed (unused) files.

	* src/fileformat.cpp, src/fileformat.h: Use a const char * title
	and declare ReadMolecule and WriteMolecule methods static to allow
	calling w/o needing to instantiate a specific OBFileFormat object.
	
	* src/*.cpp: Adjust all ReadFile methods to use const char* title,
	as above.

	* src/atom.cpp, src/binary.cpp, src/binary.h: Fixes imported from
	OELib.

	* src/generic.cpp, src/generic.h: Add support for torsion and
	angle data and associated handlers for OBGenericData.

	* src/main.cpp: Remove obsolete headers.

	* src/mol.cpp, src/mol.h: Add OBMol::FindTorsions method to fill
	the OBTorsionData in a molecule. Remove the confusing OBPose
	class. (Conformers are conformers and not managed by a separate
	class.) Remove residue code, which has migrated to residue.cpp.

	* src/obutil.cpp: Remove unused code.

	* src/pdb.cpp: Add OBSerialNums, imported from OELib.

2002-06-11  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/Makefile.in: Add proper commandline for new bin2hex for
	PhModelData. Remove smarts.* from file lists.

	* src/Vector.cpp: Import changes contributed by Open
	Eye--improvements in random vector and matrix methods.

	* src/*.h: Improve top comment to refer to Open Babel and update
	copyright dates. Add uniform #ifdef symbols to prevent clashing
	with external header files.

	* src/alchemy.cpp, src/atom.cpp, src/balst.cpp, src/chains.cpp:
	Updates to cut down on buffer overflows, to begin addressing
	PR#565413.

	* src/data.cpp: Make OBExtensionTable types uniform. Make sure
	\0 is inserted in the proper place in the buffer when reading from
	the dataptr in OBGlobalDataBase::Init(). Fixes PR#567420.

	* src/fileformat.cpp: Use updated io_type declarations.

	* src/types.txt: Reformat to remove trailing whitespace and make
	spacing uniform. (Attempt to minimize problems with reading using
	OBTypeTable.)

	* src/obutil.[h,cpp]: Remove SmartsLexReplace, which is now
	migrated to parsmart.[h,cpp].

	* src/parsmart.[h,cpp]: Updated version from code contributed by
	Open Eye. Includes caching and RestrictedMatch method and seems to
	solve memory leak problems in older version. Solves PR#536943 and
	#564096. Includes code from smarts.* files.
	
	* src/phmodel.cpp: Update to work with new SMARTS parser.

2002-06-07  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* Makefile.in: Make sure directories are properly created on "make
	install" which fixes PR#565412.

	* src/Makefile.in: Add support for rebuilding the binary data
	headers when needed. Also generates library first, then just links
	this into main.o for the program. Also as above.

	* src/bin2hex.pl: Import a revised version contributed by Open
	Eye.

	* src/atomtyp.txt, src/extable.txt, src/types.txt: Ditto.

	* src/aromatic.h, src/atomtyp.h, src/element.h, src/extable.h,
	src/phmodel.h, src/resdata.h, src/types.h: Generate using new
	bin2hex.pl.

2002-06-06  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/main.cpp: Add support for -iMIME or -oMIME <mime_type>
	input. This should facilitate use for web applications (e.g. as a
	pass-through filter to other programs).

	* src/data.h, src/data.cpp: Add methods
	IsReadable/IsWritable(io_type).

2002-05-30  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/mol.h, src/atom.cpp (HasNonSingleBond): Add new method to
	simplify and optimize bond typing code.

	* src/mol.cpp (PerceiveBondOrders): Use it.

	* src/cml.cpp: Fix problems with arrays of OBAtom* since these
	cause casting problems in ISO C++. Still has a number of "missing
	return" warnings, though these shouldn't cause problems.

2002-04-23  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/hin.cpp (ReadHIN): Fix bugs in reading HyperChem files as
	contributed by Tommi Hassinen.

2002-04-16  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/qchem.cpp: Fix conversion from internal coordinates at the
	end of a geometry optimization (silly thinko, forgot to convert
	from deg. to rad. for angles).

2002-04-12  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/bin2hex.pl: Set output header files to be "static const
	char" instead of plain char. Should prevent some clobbering,
	e.g. PR #543314.

	* src/aromatic.h, src/atomtyp.h, src/element.h, src/extable.h,
	src/phmodeldata.h, src/resdata.h, src/types.h: Generate as above.

	* src/aromatic.txt, src/atomtyp.txt, src/phmodel.txt: Change
	"OELib" in headers to "Open Babel."

	* src/data.h: Make all _dataptr members "const char" and set
	ParseLine() methods to operate on const char*.

	* src/data.cpp, src/pdb.cpp, src/phmodel.cpp, src/rotor.cpp,
	src/typer.cpp: As above.

	* src/mol.h: Set tokenize() declarations to operate on const
	char*.

	* src/tokenst.cpp: As above.

2002-04-05  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/parsmart.cpp: Added some defensive programming when calling
	delete, also set the pointer to NULL.

	* src/parsmart.h: Add constructor and destructor to
	OBSmartsParser.

	* src/smarts.h: Add full destructors for OBEdge, OBNode.

	* src/smarts.cpp: As above. Should solve memory leaks, PR #508056.
	
2002-03-30  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/*.h: Update headers to use std:: namespace to refer to
	vectors, strings, pairs, streams, etc.

	* src/*.cpp: As a first stage in full namespace support, add a
	full "using namespace std;" (to be followed eventually with only
	importing the necessary namespace portions.
	
2002-03-27  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/cml.cpp: Fix use of strcmp() in code in favor of strncmp()
	to prevent buffer overflow. Fixes problems with 2D/3D export.
	
2002-03-26  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/cml.cpp: New CML support, including built-in XML parser,
	contributed by Peter Murray-Rust.

	* src/Makefile.in: Add cml.cpp to list to support CML compilation.

	* src/data.h, src/data.cpp: Add CML file type.

	* src/extable.txt, src/extable.h: Update to add in .cml files.

	* src.fileformat.cpp, src/fileformat.h: Add in support for CML
	read/write as well as passing options to WriteMolecule methods.

	* src/main.cpp: Add support for passing output options (currently
	only XML/CML) to export methods.
	
	* src/aromatic.txt: Enable support for recognizing Se in aromatic
	rings. (Uncommon, but I've been doing work with selenophene
	recently.)

	* src/pdb.cpp (ParseConectRecord): Fix bugs in CONECT records with
	PDB files. Would only read first 4 columns and wouldn't handle
	multiple bonds. Fixes PR #529744. Use PerceiveBondOrders() when
	using ConnectTheDots() for missing bonds.
	
	* src/ghemical.cpp: Add error checking for reading bonds. Don't
	assume the # of bonds given in the header is correct.

	* src/mol.cpp: Check to see if a bond exists before adding
	it. (Prevents a bond addition (a,b) followed by (b,a).)

2002-02-22  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/Makefile.in: Add viewmol.o.

	* src/viewmol.cpp: New code to read/write ViewMol files as ported
	from the GPL'ed patch sent by J?rg-R?diger Hill.

	* src/data.h, src/data.cpp, src/extable.txt, src/extable.h: Add
	ViewMol format and enable Chem3D and derived formats.

	* src/c3d.cpp: Clean up Read methods to match normal declarations.

	* src/fileformat.h, src/fileformat.cpp: Enable Read/Write for
	Chem3d1, Chem3d2, Mmads and ViewMol formats.

	* src/main.cpp: Add support for using '--' to read from STDIN and
	output to STDOUT. Tackles request #519085.
	
2002-02-17  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* README, configure.in, configure: Update for next snapshot of
	1.99.1.

	* src/Makefile.in: Add -I(builddir) to ensure that
	configure-generated version.h is included.

	* src/alchemy.cpp, src/amber.cpp, src/balst.cpp, src/box.cpp,
	src/cacao.cpp, src/car.cpp, src/ccc.cpp, src/dmol.cpp,
	src/feat.cpp, src/gamess.cpp, src/gaussian.cpp, src/hin.cpp,
	src/jaguar.cpp, src/mopac.cpp, src/mpqc.cpp, src/nwchem.cpp,
	src/qchem.cpp, src/report.cpp, src/unichem.cpp, src/xyz.cpp:
	Make sure that SetTitle() is called with the default title if the
	format doesn't have a title, do not set atom types unless the
	format supports them (the internal atom typer is better), and for
	formats without bond information, call the new
	mol.PerceiveBondOrders() method.

	* src/mol.h, src/mol.cpp: Rename PerceiveBonds in favor of
	more accurate PerceiveBondOrders(). Make sure that
	auto-hybridization is turned off while bond orders are being
	perceived and then let it run afterwards.
	
2002-02-16  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/mol.cpp (PerceiveBonds): Improve accuracy by picking the
	closer atom (shorter bond) when there are two atoms with equal
	electronegativities.

	* src/main.cpp: Add -v flag to give version information but not
	the whole usage report.
	
	* src/base.h, src/binary_io.h, src/bitvec.h, src/commandline.h,
	src/ctransform.h, src/data.h, src/fileformat.h, src/grid.h,
	src/matrix.h, src/mol.h, src/obutil.h, src/smarts.h: Remove
	"using" statements from header files which pollute namespace of
	user code. Fixes PR#493388.

2002-02-14  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/data.h, src/data.cpp: Add MIME support to
	OBExtensionTable. Progress towards feature request #511533.

	* src/extable.txt, src/extable.h: As above (using as much from the
	current Chemical MIME page as possible).

	* src/mol.h, src/mol.cpp (PerceiveBonds): Add preliminary support
	for assigning bond orders based on Roger Sayle's "Cruft to
	Content" algorithm. Still needs functional group and aromatic ring
	recognition passes. Progress towards feature request #514589.
	<http://www.daylight.com/meetings/mug01/Sayle/m4xbondage.html>

	* src/mopac.cc: Fix bugs in read/write methods as pointed out by
	Radek Liboska. Should solve PR #515884.

	* src/atom.cpp (GetType): Fix bug with atom typing by assigning
	the ATN type of an atom if the atom typer didn't pick one before.

	* src/gamess.cpp, src/xyz.cpp, src/unichem.cpp, src/ghemical.cpp:
	Don't try to translate atom types as this causes problems. (Works
	cleanly now that the GetType() bug is fixed.)
	
2002-02-13  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/data.cpp, src/data.h, src/element.txt: Added support for
	electronegativity of elements. (Needed for assigning bond types,
	among other things.)

	* src/element.h: Regenerate using bin2hex.pl.

	* src/ghemical.cpp: Don't try to type atoms, it only causes
	problems with the internal atom typing.

2002-02-07  Geoff Hutchison  <hutchisn@chem.northwestern.edu>
	
	* src/Makefile.in: Add c3d.cpp to the compilation list.

	* src/obutil.cpp: Migrate quat.c code here temporarily.
	
	* src/mol.cpp, src/mol.h: Resolve problems with declarations of
	old quat.c code since they now fall under OpenBabel namespace.

2002-02-01  Michael Banck  <mbanck@gmx.net>

	* test/Makefile.in, src/Makefile.in: Delete Makefile when calling
	distclean target.

	* Makefile.in: Fix install target to create the $(mandir) directory.

2002-02-01  Michael Banck  <mbanck@gmx.net>

	* Makefile.in: Added manpage to install target.

	* configure, configure.in, Makefile.in, src/Makefile.in,
	src/version.h.in: Put version information in configure.in and made
	generated tarballs contain it.
	
2002-01-24  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* doc/FAQ, doc/babel.1: Minor documentation fixes.
	
2002-01-21  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* Makefile.in: Fix problems with the dist target.

	* configure, configure.in: Add 'top_builddir' macro to know where
	objdir is in the case that we're not building in the source
	directory.

	* src/Makefile.in: Run the library through 'ranlib', which fixes
	problems on Mac OS X. Fix problems running 'install' when not
	building in the source directory.

	* src/unichem.cpp: Fix bug in reading files.

	* src/version.h: Set version to 1.99.
	
	* test/Makefile.in: Fix problems with the check target: tests need
	to be copied to the build dir if we're not in the source directory.
	
2002-01-21  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* doc/*: Move documentation to subdirectory--since there will
	eventually be more than these few morsels.

	* doc/Migration: Add a Migration guide for OELib -> Open Babel changes.

2002-01-20  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* FAQ, README: Added FAQ.

2002-01-20  Michael Banck  <mbanck@gmx.net>	

	* babel.1: Initial import of man page, this is only a rough
	draft. Should become babel.1.in perhaps, with the available file
	types determined during the build.

2002-01-18  Michael Banck  <mbanck@gmx.net>	

	* Makefile.in, src/Makefile.in: Fixed install target. Removed
	Primer.html, bin2hex.pl and GNULICENSE from datafiles. Added
	datafiles to install: target.

	* Makefile.in, src/data.cpp, src/Makefile.in, src/patty.cpp,
	src/pdb.cpp, src/phmodel.cpp, src/rotor.cpp, src/typer.cpp:
	Openbabel now first checks for $BABEL_DATADIR, then looks in
	${pkgdatadir} and finally uses the compiled-in values to get the
	relevant data.
	
2002-01-16  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* src/main.cpp: Use the filename as the default filename when
	reading molecules.

	* src/qchem.cpp: Disable reading final z-matrix--currently
	produces very bizzare results for some files.

	* test/Makefile.in: Make sure the main program is compiled too.

	* test/main.cpp, test/ringtest.cpp, test/smartstest.cpp: Update to
	new OpenBabel API. (Took about a minute to change OE->OB and
	change to iterate over OBNodeBase and OBEdgeBase.)

	* test/smartstest.txt: Get rid of DOS line endings--lines are not
	properly read.

2002-01-11  Michael Banck  <mbanck@gmx.net>

	* src/Makefile.in, test/Makefile.in, Makefile.in: Added
	distclean-targets to src and test's Makefiles and removed stamp.h
	again (D'OH).

	* Makefile.in: Added 'stamp.h' to distclean:-target and renamed
	DISTDIR from *openbabel to *openbabel-dist.  ('make dist' fried my
	CVS-repository because it was named 'openbabel') dist:-target
	still does not work properly.
	Removed "\" at the end of first line of distclean:-target

2002-01-09  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* */*: Updated build environment with src/ directory and autoconf
	configure script. Builds cleanly--haven't tried using the "install
	target" The test/ subdir needs fixing to remove OELib references.

2001-12-14  Michael Banck  <mbanck@gmx.net>

	* Makefile.in: Changed oeutil.* to the new name obutil.*.

2001-12-14  Geoff Hutchison  <hutchisn@chem.northwestern.edu>

	* */*: Changed OE prefix to OB -- think I got all the cases, but
	some may crop up.

	Thanks to Michael Banck:
	-Added some GNUish files like AUTHORS, etc.
	-Added configure.in -- will shortly replace configure script with
	an autoconf one. 
	-Need to decide how to work ChangeLog--should it go before
	OELib->OpenBabel change?

	* */*: Many changes to namespace and cleanup of API to work with
	gcc v. 3. Still need to do big OE prefix cleanup.

2001-11-28  ghutchis  <ghutchis@hydra.chem.northwestern.edu>

	* configure, data.cpp, dmol.cpp, main.cpp, mol.cpp, Primer.html,
	Makefile.in: Fixes from Ghemical and a first start at the "Open
	Babel" nomenclature.

2001-11-27  ghutchis  <ghutchis@hydra.chem.northwestern.edu>

	* Import from OE current CVS. (Last GPL'ed version.)

