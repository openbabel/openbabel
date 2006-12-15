INCLUDES = -I$(top_srcdir)/include
AM_CPPFLAGS =  @XML_CPPFLAGS@
AM_LDFLAGS = @LIBS@ @XML_LIBS@

if BUILD_SHARED

pkglib_LTLIBRARIES = \
  cmlreactlformat.la cmlformat.la \
	pubchem.la xmlformat.la
# cdxmlformat.la

# cdxmlformat_la_SOURCES = cdxmlformat.cpp xml.cpp ../obmolecformat.cpp
# cdxmlformat_la_LDFLAGS = -module -avoid-version 

cmlformat_la_SOURCES = cmlformat.cpp xml.cpp ../obmolecformat.cpp
cmlformat_la_LDFLAGS = -module -avoid-version 

cmlreactlformat_la_SOURCES = cmlreactlformat.cpp xml.cpp ../obmolecformat.cpp
cmlreactlformat_la_LDFLAGS = -module -avoid-version

pubchem_la_SOURCES    = pubchem.cpp xml.cpp ../obmolecformat.cpp
pubchem_la_LDFLAGS    = -module -avoid-version

xmlformat_la_SOURCES  = xmlformat.cpp xml.cpp ../obmolecformat.cpp
xmlformat_la_LDFLAGS  = -module -avoid-version

endif # build shared