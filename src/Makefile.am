EXTRA_DIST = \
	strncasecmp.c snprintf.c snprintf.h zipstream.h zipstreamimpl.h \
	CMakeLists.txt RDKitConv.cpp config.h.cmake

if BUILD_SHARED
thisdir = .
endif
SUBDIRS = math fingerprints forcefields descriptors ops $(thisdir) formats 
if BUILD_WXGUI
SUBDIRS += GUI
endif

INCLUDES = -I$(top_builddir)/include -I$(top_srcdir)/data/ -I$(top_srcdir)/include

lib_LTLIBRARIES = libopenbabel.la

if BUILD_SHARED_UNIX
# Unix, e.g., Linux, Mac, BSD, Solaris, etc.
   dlhandler_unix = dlhandler_unix.cpp
else 
	 dlhandler_unix = 
endif

if BUILD_SHARED_WIN
# MinGW
   dlhandler_win32 = dlhandler_win32.cpp
else
	 dlhandler_win32 = 
endif

libopenbabel_la_SOURCES = \
		atom.cpp alias.cpp base.cpp bitvec.cpp bond.cpp \
		bondtyper.cpp canon.cpp builder.cpp \
		chains.cpp chiral.cpp data.cpp descriptor.cpp \
		fingerprint.cpp forcefield.cpp format.cpp generic.cpp \
		griddata.cpp grid.cpp \
		locale.cpp \
		kekulize.cpp matrix.cpp mol.cpp \
		molchrg.cpp \
		obconversion.cpp oberror.cpp obiter.cpp \
		obmolecformat.cpp obutil.cpp \
		op.cpp ops.cpp \
		parsmart.cpp \
		patty.cpp phmodel.cpp plugin.cpp \
		pointgroup.cpp \
		rand.cpp residue.cpp ring.cpp \
		rotor.cpp rotamer.cpp \
		tokenst.cpp typer.cpp \
		transform.cpp $(dlhandler_unix) $(dlhandler_win32)

libopenbabel_la_LDFLAGS = -version-info 3:3:0 -no-undefined

if BUILD_SHARED
AM_CPPFLAGS = -DOB_MODULE_PATH=\"$(pkglibdir)/$(VERSION)\" -DUSING_DYNAMIC_LIBS

libopenbabel_la_LIBADD = math/libmath.la fingerprints/libfingerprints.la \
			forcefields/libforcefields.la ops/libops.la \
			descriptors/libdescriptors.la @LTLIBOBJS@ @LIBDL@
else !BUILD_SHARED
libopenbabel_la_LIBADD = math/libmath.la fingerprints/libfingerprints.la \
			formats/libformats.la forcefields/libforcefields.la \
			ops/libops.la descriptors/libdescriptors.la @LTLIBOBJS@
endif

if BUILD_PCH
BUILT_SOURCES = all.h.gch

PCHFLAGS=-Winvalid-pch -x c++-header $(CPPFLAGS) $(INCLUDES)
all.h.gch: all.h Makefile \
	$(srcdir)/../include/openbabel/babelconfig.h \
  $(srcdir)/../include/openbabel/base.h \
	$(srcdir)/../include/openbabel/generic.h \
	$(srcdir)/../include/openbabel/mol.h \
	$(srcdir)/../include/openbabel/oberror.h
	rm -f $@
	$(CXX) $(PCHFLAGS) $<
AM_CXXFLAGS= -include $(top_srcdir)/src/all.h
clean-local:
	-rm -f all.h.gch
endif
