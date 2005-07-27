EXTRA_DIST = openbabel.i perl python

#man_MANS = babel.1

if BUILD_SWIG
BUILT_SOURCES = perl/openbabel_perl.cpp \
		python/openbabel_python.cpp

perl/openbabel_perl.cpp: openbabel.i $(top_srcdir)/src/mol.h  $(top_srcdir)/src/obconversion.h
	$(SWIG) -module "Chemistry::OpenBabel" -perl5 -c++ -o $@ \
		-I$(top_srcdir)/src $<

python/openbabel_python.cpp: openbabel.i $(top_srcdir)/src/mol.h  $(top_srcdir)/src/obconversion.h
	$(SWIG) -python -modern -c++ -o $@ -I$(top_srcdir)/src $<
endif
