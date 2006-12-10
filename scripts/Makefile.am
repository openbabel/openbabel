EXTRA_DIST = openbabel-perl.i perl \
					openbabel-python.i python \
					openbabel-ruby.i ruby \
					openbabel-java.i java

dist-hook:
	find $(distdir) -name .svn | xargs rm -rf
	find $(distdir) -name "*~" | xargs rm -rf
	(cd $(distdir)/perl; $(MAKE) distclean)
	(cd $(distdir)/python; $(PYTHON) setup.py clean)
	(cd $(distdir)/python; rm -rf build; rm -rf *.pyc *.so)
	(cd $(distdir)/ruby; $(MAKE) distclean)

# build the scripting language interfaces if --enable-maintainer-mode was set
# and SWIG is available -- please use the most recent version of SWIG

if MAINTAINER_MODE
if BUILD_SWIG
BUILT_SOURCES = perl/openbabel_perl.cpp python/openbabel_python.cpp \
	java/openbabel_java.cpp

script_includes = $(top_srcdir)/include/openbabel
script_headers = $(script_includes)/mol.h $(script_includes)/obconversion.h \
	$(script_includes)/atom.h $(script_includes)/bond.h \
	$(script_includes)/residue.h

perl/openbabel_perl.cpp: openbabel-perl.i $(script_headers)
	$(SWIG) -perl5 -c++ -naturalvar -o $@ \
		-I$(top_srcdir)/include $<;
	sed -e "s/^@EXPORT.*/& sub dl_load_flags { 0x01 }/" <perl/OpenBabel.pm >perl/OpenBabel.new;
	mv perl/OpenBabel.new perl/OpenBabel.pm;

python/openbabel_python.cpp: openbabel-python.i $(script_headers)
	$(SWIG) -python -modern -naturalvar -fastinit \
		-c++ -o $@ -I$(top_srcdir)/include $<;
	echo "import sys" >python/ob.py;
	echo "if sys.platform.find(\"linux\") != -1:"  >>python/ob.py;
	echo "    import dl" >>python/ob.py;
	echo "    sys.setdlopenflags(sys.getdlopenflags() | dl.RTLD_GLOBAL)" >>python/ob.py;
	echo >>python/ob.py;
	cat python/openbabel.py >>python/ob.py;
	mv python/ob.py python/openbabel.py;

java/openbabel_java.cpp: openbabel-java.i $(script_headers)
	$(SWIG) -java -naturalvar -c++ -o $@ -I$(top_srcdir)/include $<;
	(cd java; \
	  for file in SWIGTYPE_p_*; \
		do tgt=$$(echo $${file} | sed -e "s/SWIGTYPE_p_//" | sed -e "s/std__//g" | sed -e "s/OpenBabel__//g" | sed -e "s/p_/p/g"); \
		mv $$file $$tgt; done)

if BUILD_RUBY
BUILT_SOURCES += ruby/openbabel_ruby.cpp

ruby/openbabel_ruby.cpp: openbabel-ruby.i $(script_headers)
	$(SWIG) -ruby -naturalvar -c++ -o $@ -I$(top_srcdir)/include $<;
	(cd ruby; ruby extconf.rb --with-openbabel-include=$(includedir)/openbabel-2.0/openbabel --with-openbabel-2.0-include=$(includedir)/openbabel-2.0 --with-openbabel-lib=$(libdir) )
endif # BUILD_RUBY

endif # BUILD_SWIG
endif # MAINTAINER_MODE