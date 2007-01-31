SUBDIRS = math

libopenbabelincludedir = $(includedir)/openbabel-2.0/openbabel
libopenbabelinclude_HEADERS = \
		atom.h \
		babelconfig.h base.h bitvec.h bond.h \
		bondtyper.h canon.h \
		chains.h chiral.h data.h \
		dlhandler.h fingerprint.h forcefield.h \
		generic.h grid.h internalcoord.h kinetics.h \
		matrix.h mol.h molchrg.h \
		newlinebuf.h \
		obconversion.h oberror.h obiter.h \
		obmolecformat.h obutil.h \
		parsmart.h patty.h phmodel.h pluginiter.h \
		rand.h reaction.h \
		residue.h ring.h rotamer.h rotor.h typer.h \
		xml.h

# Remove the configure-generated babelconfig.h
dist-hook:
	rm $(distdir)/babelconfig.h
