EXTRA_DIST = OBGUI.h optswx.h selformats.h stdwx.h \
							OpenBabelGUI.html babel.xpm babel2.ico OBGUI.rc

INCLUDES = -I$(top_builddir)/include -I$(top_srcdir)/data -I$(top_srcdir)/include

AM_CPPFLAGS =  @WX_CPPFLAGS@
AM_CXXFLAGS =  @WX_CXXFLAGS@

bin_PROGRAMS = wxbabel
wxbabel_SOURCES = OBGUI.cpp optswx.cpp stdwx.cpp
wxbabel_LDADD = ../libopenbabel.la  @LIBDL@ @WX_LIBS@
