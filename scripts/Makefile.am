EXTRA_DIST = openbabel.i perl python

dist-hook:
	find $(distdir)/perl -name CVS | xargs rm -rf
	(cd $(distdir)/perl; $(MAKE) distclean)
	find $(distdir)/python -name CVS | xargs rm -rf
	(cd $(distdir)/python; $(PYTHON) setup.py clean)
	(cd $(distdir)/python; rm -rf build; rm -rf *.pyc *.so)
	find $(distdir) -name "*~" | xargs rm -rf

# build the Perl and Python wrappers if --enable-maintainer-mode was set
# and SWIG is available -- please use the most recent version of SWIG

if MAINTAINER_MODE
if BUILD_SWIG
BUILT_SOURCES = perl/openbabel_perl.cpp \
		python/openbabel_python.cpp \
	 ruby/openbabel_ruby.cpp

perl/openbabel_perl.cpp: openbabel.i $(top_srcdir)/src/mol.h  $(top_srcdir)/src/obconversion.h
	$(SWIG) -perl5 -c++ -naturalvar -o $@ \
		-I$(top_srcdir)/src $<;
	sed -e "s/^@EXPORT.*/& sub dl_load_flags { 0x01 }/" <perl/OpenBabel.pm >perl/OpenBabel.new;
	mv perl/OpenBabel.new perl/OpenBabel.pm;

python/openbabel_python.cpp: openbabel.i $(top_srcdir)/src/mol.h  $(top_srcdir)/src/obconversion.h
	$(SWIG) -python -modern -naturalvar -fastinit \
		-c++ -o $@ -I$(top_srcdir)/src $<;
	echo "import sys" >python/ob.py;
	echo "if sys.platform.find(\"linux\") != -1:"  >>python/ob.py;
	echo "    import dl" >>python/ob.py;
	echo "    sys.setdlopenflags(sys.getdlopenflags() | dl.RTLD_GLOBAL)" >>python/ob.py;
	echo >>python/ob.py;
	cat python/openbabel.py >>python/ob.py;
	mv python/ob.py python/openbabel.py;

ruby/openbabel_ruby.cpp: openbabel.i $(top_srcdir)/src/mol.h  $(top_srcdir)/src/obconversion.h
	$(SWIG) -ruby -c++ -o $@ -I$(top_srcdir)/src $<;

endif # BUILD_SWIG
endif # MAINTAINER_MODE