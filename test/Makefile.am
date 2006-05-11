SUBDIRS = cmltest

noinst_PROGRAMS = matrix rings smarts unitcell aromatic \
		atom bond conversion data format internalcoord \
		mol residue formula

bin_PROGRAMS = roundtrip

if PERL_TESTS
# Use perl's "prove" script and some wrapper .pl scripts to run tests
# (This provides a nicer summary and verbose results on request)
TESTS = wrapper.sh
else
# Standard automake tests
TESTS = matrix smartstest rings unitcell \
	$(top_srcdir)/test/cml.sh $(top_srcdir)/test/test-set.sh
endif

LDADD = ../src/libopenbabel.la @LIBDL@
TESTDATADIR = $(top_srcdir)/test/files/
AM_CPPFLAGS = -I$(top_srcdir)/src -DTESTDATADIR="\"$(TESTDATADIR)\""

EXTRA_DIST = files wrapper.sh cml.sh test-set.sh \
		aromatic.pl

roundtrip_SOURCES = roundtrip.cpp

# Unit Test binaries
aromatic_SOURCES = aromatest.cpp
atom_SOURCES = atom.cpp
bond_SOURCES = bond.cpp
conversion_SOURCES = conversion.cpp
data_SOURCES = data.cpp
format_SOURCES = format.cpp
formula_SOURCES = formula.cpp
internalcoord_SOURCES = internalcoord.cpp
matrix_SOURCES = matrixtest.cpp
mol_SOURCES = mol.cpp
residue_SOURCES = residue.cpp
rings_SOURCES = ringtest.cpp
smarts_SOURCES = smartstest.cpp
unitcell_SOURCES = unitcell.cpp
