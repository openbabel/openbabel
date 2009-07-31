# Install these as version-specific data. Newer should replace older files.
tabledir     = $(pkgdatadir)/$(VERSION)
table_DATA = \
		aromatic.txt bondtyp.txt isotope.txt resdata.txt \
		atomtyp.txt element.txt isotope-small.txt phmodel.txt \
		patterns.txt torlib.txt types.txt space-groups.txt \
		SMARTS_InteLigand.txt fragments.txt \
		ghemical.prm UFF.prm ringtyp.txt \
		mmffoop.par mmffchg.par mmfftor.par mmffang.par \
		mmffpbci.par mmffdef.par \
		mmffvdw.par mmffstbn.par mmffbndk.par mmffprop.par mmffbond.par \
		mmffdfsb.par mmffs_oop.par mmffs_tor.par mmff94.ff mmff94s.ff\
		logp.txt psa.txt mr.txt plugindefines.txt MACCS.txt

BUILT_SOURCES = aromatic.h atomtyp.h bondtyp.h element.h \
		isotope.h phmodeldata.h resdata.h \
		torlib.h types.h

EXTRA_DIST = \
	bin2hex.pl CMakeLists.txt $(table_DATA) $(BUILT_SOURCES)

if BIN2HEX
aromatic.h:	bin2hex.pl aromatic.txt
		(cd $(srcdir); ./bin2hex.pl aromatic.txt AromaticData >aromatic.h)

atomtyp.h:	bin2hex.pl atomtyp.txt
		(cd $(srcdir); ./bin2hex.pl atomtyp.txt AtomTypeData >atomtyp.h)

bondtyp.h:	bin2hex.pl bondtyp.txt
		(cd $(srcdir); ./bin2hex.pl bondtyp.txt BondTypeData >bondtyp.h)

element.h:	bin2hex.pl element.txt
		(cd $(srcdir); ./bin2hex.pl element.txt ElementData >element.h)

isotope.h:	bin2hex.pl isotope-small.txt
		(cd $(srcdir); ./bin2hex.pl isotope-small.txt IsotopeData ISOTOPE >isotope.h)

phmodeldata.h:	bin2hex.pl phmodel.txt
		(cd $(srcdir); ./bin2hex.pl phmodel.txt PhModelData PHMODELDATA >phmodeldata.h)

resdata.h:	bin2hex.pl resdata.txt
		(cd $(srcdir); ./bin2hex.pl resdata.txt ResidueData >resdata.h)

torlib.h:	bin2hex.pl torlib.txt
		(cd $(srcdir); ./bin2hex.pl torlib.txt TorsionDefaults >torlib.h)

types.h:	bin2hex.pl types.txt
		(cd $(srcdir); ./bin2hex.pl types.txt TypesData >types.h)
endif