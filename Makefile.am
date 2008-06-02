# Use our libtool macros
ACLOCAL_AMFLAGS='-I autoconf'

# only build in the scripts subdirectory if possible (perl, python, swig, etc.)
if BUILD_SCRIPTS
  SCRIPT_DIR = scripts
endif

SUBDIRS = data src test tools doc include $(SCRIPT_DIR)
EXTRA_DIST = Doxyfile windows-vc2005
AUTOMAKE_OPTIONS = foreign

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = openbabel-2.0.pc

docs:
	(cd $(top_srcdir); $(DOXYGEN))

dist-hook:
	find $(distdir) -name ".svn" | xargs rm -rf
