SUBDIRS = math

libopenbabelincludedir = $(includedir)/openbabel-2.0/openbabel
libopenbabelinclude_HEADERS = \
		atom.h \
		babelconfig.h base.h bitvec.h bond.h \
		bondtyper.h canon.h builder.h \
		chains.h chiral.h data.h descriptor.h \
		dlhandler.h fingerprint.h forcefield.h format.h \
		generic.h grid.h griddata.h groupcontrib.h \
		internalcoord.h kinetics.h \
		lineend.h  matrix.h mol.h molchrg.h \
		obconversion.h oberror.h obiter.h \
		obmolecformat.h obutil.h \
		op.h dock.h \
		parsmart.h patty.h phmodel.h plugin.h \
		rand.h reaction.h \
		residue.h ring.h rotamer.h rotor.h typer.h \
		xml.h

# Remove the configure-generated babelconfig.h
dist-hook:
	rm $(distdir)/babelconfig.h
