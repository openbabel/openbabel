SUBDIRS = openbabel

libopenbabelincludedir = $(includedir)/openbabel-2.0/openbabel
libopenbabelinclude_HEADERS = chemdrawcdx.h

if BUILD_INCHI
inchiincludedir = $(includedir)/inchi102
inchiinclude_HEADERS = inchi_api.h
endif

EXTRA_DIST = inchi102 OB-BGL RDKitConv.h CMakeLists.txt
# .svn directories will be cleaned up by top-level dist-hook
