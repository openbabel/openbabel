INCLUDES = -I$(top_builddir)/include -I$(top_srcdir)/include 
AM_CPPFLAGS =  @XML_CPPFLAGS@
AM_LDFLAGS = @LIBS@ @XML_LIBS@

if BUILD_SHARED

pkglib_LTLIBRARIES = \
  cmlformat.la \
	pubchem.la xmlformat.la cdxmlformat.la

cdxmlformat_la_SOURCES = cdxmlformat.cpp xml.cpp ../obmolecformat.cpp
cdxmlformat_la_LDFLAGS = -module -avoid-version $(AM_LDFLAGS)

cmlformat_la_SOURCES = cmlformat.cpp xml.cpp ../obmolecformat.cpp
cmlformat_la_LDFLAGS = -module -avoid-version  $(AM_LDFLAGS)

pubchem_la_SOURCES    = pubchem.cpp xml.cpp ../obmolecformat.cpp
pubchem_la_LDFLAGS    = -module -avoid-version  $(AM_LDFLAGS)

xmlformat_la_SOURCES  = xmlformat.cpp xml.cpp ../obmolecformat.cpp
xmlformat_la_LDFLAGS  = -module -avoid-version  $(AM_LDFLAGS)

endif # build shared
