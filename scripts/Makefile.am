EXTRA_DIST = openbabel-perl.i perl \
					openbabel-python.i python \
					openbabel-ruby.i ruby \
					openbabel-java.i java \
		   		openbabel-mono.i csharp	

# Remove various files before making a distribution archive
dist-hook:
	find $(distdir) -name .svn | xargs rm -rf
	find $(distdir) -name "*~" | xargs rm -rf
	-(cd $(distdir)/perl; $(MAKE) distclean)
	(cd $(distdir)/python; $(PYTHON) setup.py clean)
	(cd $(distdir)/python; rm -rf build; rm -f *.pyc *.so)
	-(cd $(distdir)/ruby; $(MAKE) distclean)
	(cd $(distdir)/java; rm -f *.so *.o *.class)
	(cd $(distdir)/java; rm -f *.jnilib)
	(cd $(distdir)/java; rm -rf META_INF)
	(cd $(distdir)/java; find . -name "*.java" -print | grep -v OBTest | xargs rm -f)

# build the scripting language interfaces if --enable-maintainer-mode was set
# and SWIG is available -- please use the most recent version of SWIG

if MAINTAINER_MODE
if BUILD_SWIG
BUILT_SOURCES = perl/openbabel_perl.cpp python/openbabel_python.cpp

script_includes = $(top_srcdir)/include/openbabel
script_headers = $(script_includes)/mol.h $(script_includes)/obconversion.h \
	$(script_includes)/atom.h $(script_includes)/bond.h \
	$(script_includes)/residue.h

perl/openbabel_perl.cpp: openbabel-perl.i $(script_headers)
	$(SWIG) -small -O -perl5 -c++ -templatereduce -naturalvar -o $@ \
		-I$(top_srcdir)/include -I$(top_builddir)/include -outdir perl $<;
	sed -e "s/^@EXPORT.*/& sub dl_load_flags { 0x01 }/" <perl/OpenBabel.pm >perl/OpenBabel.new;
	mv perl/OpenBabel.new perl/OpenBabel.pm;

python/openbabel_python.cpp: openbabel-python.i $(script_headers)
	$(SWIG) -small -O -python -templatereduce -naturalvar -outdir python \
		-c++ -o $@ -I$(top_srcdir)/include -I$(top_builddir)/include $<;
	echo "import sys" >python/ob.py;
	echo "if sys.platform.find(\"linux\") != -1:"  >>python/ob.py;
	echo "    try:" >>python/ob.py;
	echo "        import dl" >>python/ob.py;
	echo "    except ImportError:" >>python/ob.py;
	echo "        import DLFCN as dl" >>python/ob.py;
	echo "    sys.setdlopenflags(sys.getdlopenflags() | dl.RTLD_GLOBAL)" >>python/ob.py;
	echo >>python/ob.py;
	cat python/openbabel.py >>python/ob.py;
	mv python/ob.py python/openbabel.py;

if BUILD_JAVA
BUILT_SOURCES += java/openbabel_java.cpp

java/openbabel_java.cpp: openbabel-java.i $(script_headers)
	(cd java; rm -rf org; mkdir org; mkdir org/openbabel)
	$(SWIG) -small -O -java -package org.openbabel \
	  -outdir java/org/openbabel -naturalvar -templatereduce -c++ -o $@ \
	  -I$(top_srcdir)/include -I$(top_builddir)/include $<;
	$(JAVAC) java/org/openbabel/*.java
	(cd java; ls org/openbabel/*.class > classes.list)
	(cd java; $(JAR) cf openbabel.jar @classes.list)
	(cd java; rm -rf org classes.list)
endif BUILD_JAVA

if BUILD_CSHARP
BUILT_SOURCES += csharp/openbabel_csharp.cpp

csharp/openbabel_csharp.cpp: openbabel-mono.i $(script_headers)
	$(SWIG) -small -O -csharp -templatereduce -c++ -o $@ \
	-namespace OpenBabel -outdir csharp \
	-I$(top_srcdir)/include -I$(top_builddir)/include $<;
	(cd csharp; mv example.cs example.tmp);
	(cd csharp; $(CSHARP) /optimize /target:library /out:OBDotNet.dll *.cs);
	(cd csharp; rm *.cs; mv example.tmp example.cs);
endif # BUILD_CSHARP

if BUILD_RUBY
BUILT_SOURCES += ruby/openbabel_ruby.cpp

ruby/openbabel_ruby.cpp: openbabel-ruby.i $(script_headers)
	$(SWIG) -small -O -ruby -naturalvar -autorename -c++ -o $@ \
	-templatereduce -I$(top_srcdir)/include -I$(top_builddir)/include $<;
	sed -e 's/rb_define_module(\"Openbabel\")/rb_define_module(\"OpenBabel\")/' -e 's/void Init_OpenBabel/void Init_openbabel/' <ruby/openbabel_ruby.cpp >ruby/openbabel_ruby.new 
	mv ruby/openbabel_ruby.new ruby/openbabel_ruby.cpp
endif # BUILD_RUBY

endif # BUILD_SWIG
endif # MAINTAINER_MODE
