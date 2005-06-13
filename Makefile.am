SUBDIRS = src test doc tools data
EXTRA_DIST = Doxyfile bc32.mak windows
AUTOMAKE_OPTIONS = foreign

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = openbabel.pc

docs:
	$(DOXYGEN)