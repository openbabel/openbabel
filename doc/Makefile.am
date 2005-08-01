EXTRA_DIST = $(man_MANS) Migration.html README.dioxin.pov FAQ.html README.povray babel31.inc dioxin.mol2 dioxin.inc dioxin.png dioxin.pov

man_MANS = babel.1

if BUILD_DOCS
babel.1:	$(top_srcdir)/src/main.cpp $(srcdir)/Doxyfile-man
		$(DOXYGEN) $(srcdir)/Doxyfile-man
		mv man/man1/openbabel.1 babel.1
		rm -rf man
endif

