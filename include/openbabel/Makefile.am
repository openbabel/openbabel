SUBDIRS = math

# will not be automatically removed otherwise
DISTCLEANFILES = babelconfig.h

libopenbabelincludedir = $(includedir)/openbabel-2.0/openbabel
libopenbabelinclude_HEADERS = \
	alias.h \
	atom.h \
	atomclass.h \
	babelconfig.h \
	base.h \
	bitvec.h \
	bond.h \
	bondtyper.h \
	builder.h \
	canon.h \
	chains.h \
	chiral.h \
	data.h \
	descriptor.h \
	dlhandler.h \
	fingerprint.h \
	forcefield.h \
	format.h \
	generic.h \
	grid.h \
	griddata.h \
	groupcontrib.h \
	internalcoord.h \
	kinetics.h \
	locale.h \
	lineend.h \
	matrix.h \
	mol.h \
	molchrg.h \
	obconversion.h \
	oberror.h \
	obiter.h \
	obmolecformat.h \
	obutil.h \
	op.h \
	optransform.h \
	parsmart.h \
	patty.h \
	phmodel.h \
	plugin.h \
	pointgroup.h \
	rand.h \
	reaction.h \
	residue.h \
	ring.h \
	rotamer.h \
	rotor.h \
	tokenst.h \
	typer.h \
	xml.h

# Remove the configure-generated babelconfig.h
dist-hook:
	rm $(distdir)/babelconfig.h
