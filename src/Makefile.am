SUBDIRS = math formats

INCLUDES = -I$(srcdir)/../data/
AM_CPPFLAGS = -DOB_MODULE_PATH=\"$(pkglibdir)\" -DUSING_DYNAMIC_LIBS

lib_LTLIBRARIES = libopenbabel.la

libopenbabel_la_SOURCES = \
		atom.cpp base.cpp bitvec.cpp bond.cpp \
		bondtyper.cpp chains.cpp chiral.cpp data.cpp \
		dlhandler_unix.cpp obconversion.cpp \
		finger2.cpp fingerprint.cpp \
		fastsearch.cpp generic.cpp grid.cpp \
		kekulize.cpp \
		matrix.cpp mol.cpp molchrg.cpp oberror.cpp \
		obiter.cpp obutil.cpp parsmart.cpp \
		patty.cpp phmodel.cpp rand.cpp residue.cpp ring.cpp \
		rotor.cpp rotamer.cpp tokenst.cpp typer.cpp \
		transform.cpp

libopenbabel_la_LIBADD = math/libmath.la @LTLIBOBJS@
libopenbabel_la_LDFLAGS = -version-info 1:0:0

#libOBPython_la_SOURCES = openbabel_python.cpp
#openbabel_python.cpp: openbabel.i
#	swig -python -c++ -modern -o $@ $<

#libOBPerl_la_SOURCES = openbabel_perl.cpp
#openbabel_perl.cpp: openbabel.i
#       swig -module "Chemistry::OpenBabel" -perl5 -c++ -o $@ $<

pkginclude_HEADERS = \
		base.h bitvec.h \
		bondtyper.h chains.h chiral.h data.h \
		dlhandler.h fastsearch.h fingerprint.h \
		generic.h grid.h matrix.h mol.h molchrg.h \
		obconversion.h oberror.h obifstream.h \
		obiter.h obutil.h \
		parsmart.h patty.h phmodel.h reaction.h \
		ring.h rotor.h rotamer.h typer.h 

bin_PROGRAMS = babel
babel_SOURCES = main.cpp
babel_LDADD = libopenbabel.la  @LIBDL@ @LIBOBJS@ \
	# formats/inchi/libinchi.la

EXTRA_DIST = \
	config.h.in strncasecmp.c snprintf.c snprintf.h \
	openbabel.i openbabel.py openbabel.pm

