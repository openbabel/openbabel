SUBDIRS = cmltest

# Only compile these when running "make check"
check_PROGRAMS = aromatic atom bond cansmi \
	conversion data format formula formalcharge \
	inchiwrite \
	internalcoord invalidsmarts invalidsmiles \
	iterators logp_psa math mol residue rings smarts \
	smilesmatch unitcell

bin_PROGRAMS = roundtrip

if PERL_TESTS
# Use perl's "prove" script and some wrapper .pl scripts to run tests
# (This provides a nicer summary and verbose results on request)
TESTS = wrapper.sh
else
# Standard automake tests (if we don't have prove)
TESTS = cansmi invalidsmarts invalidsmiles smarts rings unitcell \
	$(top_srcdir)/test/cml.sh $(top_srcdir)/test/test-set.sh
endif

LDADD = ../src/libopenbabel.la @LIBDL@
TESTDATADIR = $(top_srcdir)/test/files/
AM_CPPFLAGS = -I$(top_srcdir)/include -DTESTDATADIR="\"$(TESTDATADIR)\""

EXTRA_DIST = files wrapper.sh cml.sh test-set.sh \
		aromatic.pl inchi

roundtrip_SOURCES = roundtrip.cpp

# Unit Test binaries
aromatic_SOURCES = aromatest.cpp
atom_SOURCES = atom.cpp
bond_SOURCES = bond.cpp
cansmi_SOURCES = cansmi.cpp
conversion_SOURCES = conversion.cpp
data_SOURCES = data.cpp
format_SOURCES = format.cpp
formula_SOURCES = formula.cpp
formalcharge_SOURCES = formalcharge.cpp
inchiwrite_SOURCES = inchiwrite.cpp ../src/formats/getinchi.cpp
internalcoord_SOURCES = internalcoord.cpp
invalidsmiles_SOURCES = invalidsmiles.cpp
invalidsmarts_SOURCES = invalidsmarts.cpp
iterators_SOURCES = iterators.cpp
logp_psa_SOURCES = logp_psa.cpp
math_SOURCES = math.cpp
mol_SOURCES = mol.cpp
residue_SOURCES = residue.cpp
rings_SOURCES = ringtest.cpp
smarts_SOURCES = smartstest.cpp
smilesmatch_SOURCES = smilesmatch.cpp
unitcell_SOURCES = unitcell.cpp
