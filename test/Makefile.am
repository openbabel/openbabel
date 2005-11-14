SUBDIRS = cmltest

noinst_PROGRAMS = obtest matrixtest ringtest smartstest roundtrip

LDADD = ../src/libopenbabel.la @LIBDL@
AM_CPPFLAGS = -I$(top_srcdir)/src -DTESTDATADIR="\"$(top_srcdir)/test/\""

EXTRA_DIST = smartsresults.txt smartstest.txt unitcell.txt ringresults.txt attype.00.smi cml.sh test-set.sh

TESTS = matrixtest smartstest ringtest obtest $(top_srcdir)/test/cml.sh $(top_srcdir)/test/test-set.sh

obtest_SOURCES = main.cpp unitcell.cpp
matrixtest_SOURCES = matrixtest.cpp
ringtest_SOURCES = ringtest.cpp
smartstest_SOURCES = smartstest.cpp
roundtrip_SOURCES = roundtrip.cpp
