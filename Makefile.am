SUBDIRS = src test doc tools data
EXTRA_DIST = Doxyfile bc32.mak windows scripts
AUTOMAKE_OPTIONS = foreign

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = openbabel-2.0.pc

docs:
	$(DOXYGEN)