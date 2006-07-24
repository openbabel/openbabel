EXTRA_DIST = xml.h

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/data
AM_CPPFLAGS =  @XML_CPPFLAGS@
AM_LDFLAGS = @LIBS@ @XML_LIBS@
CXXFLAGS += -O0

if BUILD_SHARED

pkglib_LTLIBRARIES = \
	cmlreactlformat.la pubchem.la xcmlformat.la xmlformat.la

cmlreactlformat_la_SOURCES = cmlreactlformat.cpp xml.cpp
cmlreactlformat_la_LDFLAGS = -module -avoid-version

pubchem_la_SOURCES    = pubchem.cpp xml.cpp
pubchem_la_LDFLAGS    = -module -avoid-version

xcmlformat_la_SOURCES = xcmlformat.cpp xml.cpp
xcmlformat_la_LDFLAGS = -module -avoid-version 

xmlformat_la_SOURCES  = xmlformat.cpp xml.cpp
xmlformat_la_LDFLAGS  = -module -avoid-version

endif # build shared