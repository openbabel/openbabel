SUBDIRS = cmltest

noinst_PROGRAMS = matrixtest ringtest smartstest roundtrip unitcell aromatest

if PERL_TESTS
# Use perl's "prove" script and some wrapper .pl scripts to run tests
# (This provides a nicer summary and verbose results on request)
TESTS = wrapper.sh
else
# Standard automake tests
TESTS = matrixtest smartstest ringtest unitcell \
	$(top_srcdir)/test/cml.sh $(top_srcdir)/test/test-set.sh
endif

LDADD = ../src/libopenbabel.la @LIBDL@
AM_CPPFLAGS = -I$(top_srcdir)/src -DTESTDATADIR="\"$(top_srcdir)/test/\""

EXTRA_DIST = smartsresults.txt smartstest.txt \
		unitcell.txt unitcell_results.txt \
		ringresults.txt attype.00.smi \
		cml.sh test-set.sh \
		matrix.pl unitcell.pl rings.pl \
		smarts.pl test-set.pl cml.pl \
		aromatest.pl

matrixtest_SOURCES = matrixtest.cpp
ringtest_SOURCES = ringtest.cpp
smartstest_SOURCES = smartstest.cpp
unitcell_SOURCES = unitcell.cpp
roundtrip_SOURCES = roundtrip.cpp
aromatest_SOURCES = aromatest.cpp