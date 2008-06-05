EXTRA_DIST = openbabel-perl.i perl \
					openbabel-python.i python \
					openbabel-ruby.i ruby \
					openbabel-java.i java \
		   		openbabel-csharp.i csharp	

# Remove various files before making a distribution archive
dist-hook:
	find $(distdir) -name .svn | xargs rm -rf
	find $(distdir) -name "*~" | xargs rm -rf
	(cd $(distdir)/perl; $(MAKE) distclean)
	(cd $(distdir)/python; $(PYTHON) setup.py clean)
	(cd $(distdir)/python; rm -rf build; rm -f *.pyc *.so)
	(cd $(distdir)/ruby; $(MAKE) distclean)
	(cd $(distdir)/java; rm -f *.so *.o *.class)
	(cd $(distdir)/java; rm -f *.jnilib *.jar)

# build the scripting language interfaces if --enable-maintainer-mode was set
# and SWIG is available -- please use the most recent version of SWIG

if MAINTAINER_MODE
if BUILD_SWIG
BUILT_SOURCES = perl/openbabel_perl.cpp \
	java/openbabel_java.cpp csharp/openbabel_csharp.cpp \
	python/obcore.cpp python/obconversion.cpp python/obtemplate.cpp

script_includes = $(top_srcdir)/include/openbabel
script_headers = $(script_includes)/mol.h $(script_includes)/obconversion.h \
	$(script_includes)/atom.h $(script_includes)/bond.h \
	$(script_includes)/residue.h

perl/openbabel_perl.cpp: openbabel-perl.i $(script_headers)
	$(SWIG) -small -O -perl5 -c++ -templatereduce -naturalvar -o $@ \
		-I$(top_srcdir)/include -I$(top_builddir)/include $<;
	sed -e "s/^@EXPORT.*/& sub dl_load_flags { 0x01 }/" <perl/OpenBabel.pm >perl/OpenBabel.new;
	mv perl/OpenBabel.new perl/OpenBabel.pm;

python/obconversion.cpp: openbabel-python2.i $(script_headers)
	$(SWIG) -small -O -python -templatereduce -naturalvar \
		-c++ -o $@ -I$(top_srcdir)/include -I$(top_builddir)/include $<;

python/obtemplate.cpp: openbabel-python3.i $(script_headers)
	$(SWIG) -small -O -python -templatereduce -naturalvar \
		-c++ -o $@ -I$(top_srcdir)/include -I$(top_builddir)/include $<;

python/obcore.cpp: openbabel-python.i $(script_headers)
	$(SWIG) -small -O -python -templatereduce -naturalvar \
		-c++ -o $@ -I$(top_srcdir)/include -I$(top_builddir)/include $<;

java/openbabel_java.cpp: openbabel-java.i $(script_headers)
	$(SWIG) -small -O -java -naturalvar -templatereduce -c++ -o $@ \
	-I$(top_srcdir)/include -I$(top_builddir)/include $<;

csharp/openbabel_csharp.cpp: openbabel-csharp.i $(script_headers)
	$(SWIG) -small -O -csharp -naturalvar -templatereduce -c++ -o $@ \
	-I$(top_srcdir)/include -I$(top_builddir)/include $<;

if BUILD_RUBY
BUILT_SOURCES += ruby/openbabel_ruby.cpp

ruby/openbabel_ruby.cpp: openbabel-ruby.i $(script_headers)
	$(SWIG) -small -O -ruby -naturalvar -autorename -c++ -o $@ \
	-templatereduce -I$(top_srcdir)/include -I$(top_builddir)/include $<;
	(cd ruby; ruby extconf.rb \
  --with-openbabel-include=$(includedir)/openbabel-2.0 \
  --with-openbabel-lib=$(libdir) )
	sed -e 's/rb_define_module(\"Openbabel\")/rb_define_module(\"OpenBabel\")/' -e 's/void Init_OpenBabel/void Init_openbabel/' <ruby/openbabel_ruby.cpp >ruby/openbabel_ruby.new 
	mv ruby/openbabel_ruby.new ruby/openbabel_ruby.cpp
endif # BUILD_RUBY

endif # BUILD_SWIG
endif # MAINTAINER_MODE
