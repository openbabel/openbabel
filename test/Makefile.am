noinst_PROGRAMS = obtest roundtrip
obtest_SOURCES = matrixtest.cpp ringtest.cpp main.cpp smartstest.cpp unitcell.cpp
roundtrip_SOURCES = roundtrip.cpp
LDADD = ../src/libopenbabel.la
EXTRA_DIST = smartsresults.txt smartstest.txt unitcell.txt ringresults.txt attype.00.smi cml.sh
TESTS = obtest cml.sh
