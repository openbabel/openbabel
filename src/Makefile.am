EXTRA_DIST = \
	strncasecmp.c snprintf.c snprintf.h zipstream.h zipstreamimpl.h \
	dlhandler_win32.cpp all.h

SUBDIRS = math formats fingerprints
if BUILD_WXGUI
SUBDIRS += GUI
endif

if BUILD_PCH
BUILT_SOURCES = all.h.gch
PCHFLAGS=-Winvalid-pch -x c++-header $(CXXFLAGS) $(INCLUDES)
all.h.gch: all.h Makefile babelconfig.h mol.h base.h bitvec.h generic.h \
	oberror.h obiter.h obconversion.h atom.h bond.h
	rm -f $@
	$(CXX) $(PCHFLAGS) $<
AM_CXXFLAGS= -include $(top_srcdir)/src/all.h

clean-local:
	-rm -f all.h.gch
endif

INCLUDES = -I$(top_srcdir)/data/

bin_PROGRAMS = babel
babel_SOURCES = main.cpp
babel_LDADD = libopenbabel.la  @LIBDL@ @LIBOBJS@

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
		newlinebuf.cpp \
		obconversion.cpp oberror.cpp obiter.cpp \
		obutil.cpp parsmart.cpp \
		patty.cpp phmodel.cpp rand.cpp residue.cpp ring.cpp \
		rotor.cpp rotamer.cpp tokenst.cpp typer.cpp \
		transform.cpp $(dlhandler)

libopenbabel_la_LDFLAGS = -version-info 3:0:1 -no-undefined

if BUILD_SHARED
AM_CPPFLAGS = -DOB_MODULE_PATH=\"$(pkglibdir)\" -DUSING_DYNAMIC_LIBS

libopenbabel_la_LIBADD = math/libmath.la fingerprints/libfingerprints.la \
			@LTLIBOBJS@ @LIBDL@
else !BUILD_SHARED
libopenbabel_la_LIBADD = math/libmath.la fingerprints/libfingerprints.la \
			formats/libformats.la @LTLIBOBJS@
endif

libopenbabelincludedir = $(includedir)/openbabel-2.0/openbabel
libopenbabelinclude_HEADERS = \
		atom.h \
		babelconfig.h base.h bond.h bitvec.h \
		bondtyper.h chains.h chiral.h data.h \
		dlhandler.h fingerprint.h \
		generic.h grid.h internalcoord.h kinetics.h \
		newlinebuf.h \
		matrix.h mol.h molchrg.h \
		obconversion.h oberror.h obiter.h \
		obutil.h \
		parsmart.h patty.h phmodel.h rand.h reaction.h \
		residue.h ring.h rotor.h rotamer.h typer.h 

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
