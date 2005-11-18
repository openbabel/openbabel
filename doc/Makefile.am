EXTRA_DIST = $(man_MANS) README.dioxin.pov README.povray babel31.inc \
		dioxin.mol2 dioxin.inc dioxin.png dioxin.pov \
		Doxyfile-man \
		FAQ.html OBTwoMigration.html \
		babel.html obchiral.html obfit.html obgrep.html \
		obprop.html obrotate.html 

man_MANS = babel.1 obfit.1 obgrep.1 obrotate.1 obprop.1 obchiral.1

# build the man pages only if --enable-maintainer-mode is used
# and Doxygen is available. Please install the most recent version of Doxygen.

# if MAINTAINER_MODE
# if BUILD_DOCS
# babel.1:	$(top_srcdir)/src/main.cpp $(srcdir)/Doxyfile-man
# 		$(DOXYGEN) $(srcdir)/Doxyfile-man
# 		mv man/man1/babel.1 .
# 		rm -rf man

# obfit.1:	$(top_srcdir)/tools/obfit.cpp $(srcdir)/Doxyfile-man
# 		sed -e 's!/src/main.cpp!/tools/obfit.cpp!' \
# 		   <$(srcdir)/Doxyfile-man | $(DOXYGEN) -
# 		mv man/man1/obfit.1 .
# 		rm -rf man

# obgrep.1:	$(top_srcdir)/tools/obgrep.cpp $(srcdir)/Doxyfile-man
# 		sed -e 's!/src/main.cpp!/tools/obgrep.cpp!' \
# 		   <$(srcdir)/Doxyfile-man | $(DOXYGEN) -
# 		mv man/man1/obgrep.1 .
# 		rm -rf man

# obprop.1:	$(top_srcdir)/tools/obprop.cpp $(srcdir)/Doxyfile-man
# 		sed -e 's!/src/main.cpp!/tools/obprop.cpp!' \
# 		   <$(srcdir)/Doxyfile-man | $(DOXYGEN) -
# 		mv man/man1/obprop.1 .
# 		rm -rf man

# obrotate.1:	$(top_srcdir)/tools/obrotate.cpp $(srcdir)/Doxyfile-man
# 		sed -e 's!/src/main.cpp!/tools/obrotate.cpp!' \
# 		   <$(srcdir)/Doxyfile-man | $(DOXYGEN) -
# 		mv man/man1/obrotate.1 .
# 		rm -rf man

# obchiral.1:	$(top_srcdir)/tools/obchiral.cpp $(srcdir)/Doxyfile-man
# 		sed -e 's!/src/main.cpp!/tools/obchiral.cpp!' \
# 		   <$(srcdir)/Doxyfile-man | $(DOXYGEN) -
# 		mv man/man1/obchiral.1 .
# 		rm -rf man

# endif
# endif

