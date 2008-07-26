AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_builddir)/include

EXTRA_DIST = CMakeLists.txt obsym.cpp
bin_PROGRAMS = babel \
		obchiral \
		obconformer \
		obenergy \
		obfit \
		obgen \
		obgrep \
		obminimize \
		obprobe \
		obprop \
		obrotamer \
		obrotate 

babel_SOURCES = babel.cpp
babel_LDADD = ../src/libopenbabel.la @LIBDL@

obchiral_SOURCES = obchiral.cpp
obchiral_LDADD = ../src/libopenbabel.la @LIBDL@

obconformer_SOURCES = obconformer.cpp
obconformer_LDADD = ../src/libopenbabel.la @LIBDL@

obenergy_SOURCES = obenergy.cpp
obenergy_LDADD = ../src/libopenbabel.la @LIBDL@

obfit_SOURCES = obfit.cpp
obfit_LDADD = ../src/libopenbabel.la @LIBDL@

obgen_SOURCES = obgen.cpp
obgen_LDADD = ../src/libopenbabel.la @LIBDL@

obgrep_SOURCES = obgrep.cpp
obgrep_LDADD = ../src/libopenbabel.la @LIBDL@

obminimize_SOURCES = obminimize.cpp
obminimize_LDADD = ../src/libopenbabel.la @LIBDL@

obprobe_SOURCES = obprobe.cpp
obprobe_LDADD = ../src/libopenbabel.la @LIBDL@

obprop_SOURCES = obprop.cpp
obprop_LDADD = ../src/libopenbabel.la @LIBDL@

obrotate_SOURCES = obrotate.cpp
obrotate_LDADD = ../src/libopenbabel.la @LIBDL@

obrotamer_SOURCES = obrotamer.cpp
obrotamer_LDADD = ../src/libopenbabel.la @LIBDL@

# obsym_SOURCES = obsym.cpp
# obsym_LDADD = ../src/libopenbabel.la @LIBDL@

# obdoc is intended to be used for generating file format documentation
# and is not a general end-user program

# obdoc_SOURCES = obdoc.cpp
# obdoc_LDADD = ../src/libopenbabel.la @LIBDL@
