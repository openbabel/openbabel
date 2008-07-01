INCLUDES = -I$(top_builddir)/include -I$(top_srcdir)/include 
AM_CPPFLAGS =  @XML_CPPFLAGS@
AM_LDFLAGS = @LIBS@ @XML_LIBS@ ../../libopenbabel.la

if USE_BOOST_FLAGS
AM_CXXFLAGS = $(BOOST_CPPFLAGS)
endif

if BUILD_SHARED
formatsdir = $(pkglibdir)/$(VERSION)
formats_LTLIBRARIES = \
	pubchem.la xmlformat.la cdxmlformat.la cmlformat.la
if HAVE_SHARED_POINTER
    formats_LTLIBRARIES += cmlreactformat.la
endif

cdxmlformat_la_SOURCES = cdxmlformat.cpp xml.cpp 
cdxmlformat_la_LDFLAGS = -module -avoid-version -no-undefined $(AM_LDFLAGS)

cmlformat_la_SOURCES = cmlformat.cpp xml.cpp 
cmlformat_la_LDFLAGS = -module -avoid-version -no-undefined $(AM_LDFLAGS)

cmlreactformat_la_SOURCES = cmlreactformat.cpp xml.cpp 
cmlreactformat_la_LDFLAGS = -module -avoid-version -no-undefined $(AM_LDFLAGS)

pubchem_la_SOURCES    = pubchem.cpp xml.cpp 
pubchem_la_LDFLAGS    = -module -avoid-version -no-undefined $(AM_LDFLAGS)

xmlformat_la_SOURCES  = xmlformat.cpp xml.cpp 
xmlformat_la_LDFLAGS  = -module -avoid-version -no-undefined $(AM_LDFLAGS)

endif # build shared
