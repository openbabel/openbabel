SUBDIRS = math formats

INCLUDES = -I$(srcdir)/../data/

lib_LTLIBRARIES = libopenbabel.la

libopenbabel_la_SOURCES = \
		atom.cpp base.cpp binary.cpp bitvec.cpp bond.cpp \
		bondtyper.cpp chains.cpp chiral.cpp data.cpp \
		obconversion.cpp fingerprint.cpp generic.cpp grid.cpp \
		kekulize.cpp matrix.cpp mol.cpp molchrg.cpp \
		oberror.cpp obutil.cpp parsmart.cpp patty.cpp \
		phmodel.cpp rand.cpp residue.cpp ring.cpp \
		rotor.cpp tokenst.cpp typer.cpp transform.cpp

libopenbabel_la_LIBADD = math/libmath.la formats/libformats.la @LTLIBOBJS@
libopenbabel_la_LDFLAGS = -version-info 0:0:0

libopenbabelincludedir = $(includedir)/openbabel
libopenbabelinclude_HEADERS = \
		babelconfig.h base.h binary.h bitvec.h \
		bondtyper.h chains.h chiral.h data.h \
		fingerprint.h generic.h grid.h \
		matrix.h mol.h molchrg.h \
		oberror.h obifstream.h obutil.h \
		parsmart.h patty.h phmodel.h ring.h \
		rotor.h typer.h 

bin_PROGRAMS = babel
babel_SOURCES = main.cpp
babel_LDADD = libopenbabel.la @LIBOBJS@

EXTRA_DIST = \
	bc32.mak borland.bat oblib.dsp oblib.dsw \
	strncasecmp.c snprintf.c snprintf.h
