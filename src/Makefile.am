EXTRA_DIST = \
	strncasecmp.c snprintf.c snprintf.h zipstream.h zipstreamimpl.h \
	dlhandler_win32.cpp

SUBDIRS = math formats fingerprints
INCLUDES = -I$(top_srcdir)/data/

lib_LTLIBRARIES = libopenbabel.la

if BUILD_SHARED
dlhandler = dlhandler_unix.cpp
else
dlhandler = 
endif

libopenbabel_la_SOURCES = \
		atom.cpp base.cpp bitvec.cpp bond.cpp \
		bondtyper.cpp chains.cpp chiral.cpp data.cpp \
		fingerprint.cpp generic.cpp grid.cpp \
		kekulize.cpp matrix.cpp mol.cpp molchrg.cpp \
		obconversion.cpp oberror.cpp obiter.cpp \
		obutil.cpp parsmart.cpp \
		patty.cpp phmodel.cpp rand.cpp residue.cpp ring.cpp \
		rotor.cpp rotamer.cpp tokenst.cpp typer.cpp \
		transform.cpp $(dlhandler)

libopenbabel_la_LDFLAGS = -no-undefined -version-info 1:2:0

if BUILD_SHARED
AM_CPPFLAGS = -DOB_MODULE_PATH=\"$(pkglibdir)\" -DUSING_DYNAMIC_LIBS

libopenbabel_la_LIBADD = math/libmath.la fingerprints/libfingerprints.la \
			@LTLIBOBJS@ @LIBDL@
else !BUILD_SHARED
libopenbabel_la_LIBADD = math/libmath.la fingerprints/libfingerprints.la \
			formats/libformats.la @LTLIBOBJS@
endif

bin_PROGRAMS = babel
babel_SOURCES = main.cpp
babel_LDADD = libopenbabel.la  @LIBDL@ @LIBOBJS@

libopenbabelincludedir = $(includedir)/openbabel-2.0/openbabel
libopenbabelinclude_HEADERS = \
		babelconfig.h base.h bitvec.h \
		bondtyper.h chains.h chiral.h data.h \
		dlhandler.h fingerprint.h \
		generic.h grid.h matrix.h mol.h molchrg.h \
		obconversion.h oberror.h obiter.h \
		obmolecformat.h obutil.h \
		parsmart.h patty.h phmodel.h reaction.h \
		ring.h rotor.h rotamer.h typer.h 

# Use sed to substitute #include "header.h" to #include <openbabel/header.h>
install-data-hook:
	for header in $(DESTDIR)$(libopenbabelincludedir)/*.h \
	      $(DESTDIR)$(libopenbabelincludedir)/math/*.h; do \
	   sed -e 's!^#include "\([^"]*\)"!#include <openbabel/\1>!' \
	      < $${header} >$${header}.new; \
	   mv $${header}.new $${header}; \
	done

# Remove the configure-generated babelconfig.h
dist-hook:
	rm $(distdir)/babelconfig.h
