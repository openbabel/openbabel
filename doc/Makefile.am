EXTRA_DIST = $(man_MANS) Migration.html README.dioxin.pov FAQ.html README.povray babel31.inc dioxin.mol2 dioxin.inc dioxin.png dioxin.pov

man_MANS = babel.1

if BUILD_DOCS
babel.1:	$(srcdir)/../src/main.cpp
		$(DOXYGEN) $(srcdir)/Doxyfile-man
endif

