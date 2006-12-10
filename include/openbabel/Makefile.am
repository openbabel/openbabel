libopenbabelincludedir = $(includedir)/openbabel-2.0/openbabel
libopenbabelinclude_HEADERS = \
		atom.h \
		babelconfig.h base.h bond.h bitvec.h \
		bondtyper.h chains.h chiral.h data.h \
		dlhandler.h fingerprint.h forcefield.h \
		generic.h grid.h internalcoord.h kinetics.h \
		newlinebuf.h \
		matrix.h mol.h molchrg.h \
		obconversion.h oberror.h obiter.h \
		obutil.h \
		parsmart.h patty.h phmodel.h rand.h reaction.h \
		residue.h ring.h rotor.h rotamer.h typer.h 

# Use sed to substitute #include "header.h" to #include <openbabel/header.h>
install-data-hook:
	for header in $(DESTDIR)$(libopenbabelincludedir)/*.h \
	      $(DESTDIR)$(libopenbabelincludedir)/math/*.h; do \
	   sed -e 's!^#include "\([^"]*\)"!#include <openbabel/\1>!' \
	      < $${header} >$${header}.new; \
	   mv $${header}.new $${header}; \
	done

# Remove the configure-generated babelconfig.h
dist-hook:
	rm $(distdir)/babelconfig.h
