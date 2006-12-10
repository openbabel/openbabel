
# obdoc is intended to be used for generating file format documentation
# and is not a general end-user program
bin_PROGRAMS = obchiral obgrep obrotate obfit obprop obxyz obgen \
							obrotamer # obdoc

obchiral_SOURCES = obchiral.cpp
obchiral_LDADD = ../src/libopenbabel.la @LIBDL@

obfit_SOURCES = obfit.cpp
obfit_LDADD = ../src/libopenbabel.la @LIBDL@

obgen_SOURCES = obgen.cpp
obgen_LDADD = ../src/libopenbabel.la @LIBDL@

obgrep_SOURCES = obgrep.cpp
obgrep_LDADD = ../src/libopenbabel.la @LIBDL@

obprop_SOURCES = obprop.cpp
obprop_LDADD = ../src/libopenbabel.la @LIBDL@

obrotate_SOURCES = obrotate.cpp
obrotate_LDADD = ../src/libopenbabel.la @LIBDL@

obrotamer_SOURCES = obrotamer.cpp
obrotamer_LDADD = ../src/libopenbabel.la @LIBDL@

obxyz_SOURCES = obxyz.cpp
obxyz_LDADD = ../src/libopenbabel.la @LIBDL@

# obdoc_SOURCES = obdoc.cpp
# obdoc_LDADD = ../src/libopenbabel.la @LIBDL@

AM_CPPFLAGS = -I$(top_srcdir)/include
