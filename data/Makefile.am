pkgdata_DATA = \
		aromatic.txt bondtyp.txt isotope.txt resdata.txt \
		atomtyp.txt element.txt isotope-small.txt phmodel.txt types.txt

EXTRA_DIST = \
	bin2hex.pl $(pkgdata_DATA)

BUILT_SOURCES = aromatic.h atomtyp.h bondtyp.h element.h extable.h \
		isotope.h phmodeldata.h resdata.h types.h

aromatic.h:	bin2hex.pl aromatic.txt
		(cd $(srcdir); ./bin2hex.pl aromatic.txt AromaticData >aromatic.h)

atomtyp.h:	bin2hex.pl atomtyp.txt
		(cd $(srcdir); ./bin2hex.pl atomtyp.txt AtomTypeData >atomtyp.h)

bondtyp.h:	bin2hex.pl bondtyp.txt
		(cd $(srcdir); ./bin2hex.pl bondtyp.txt BondTypeData >bondtyp.h)

element.h:	bin2hex.pl element.txt
		(cd $(srcdir); ./bin2hex.pl element.txt ElementData >element.h)

extable.h:	bin2hex.pl extable.txt
		(cd $(srcdir); ./bin2hex.pl extable.txt ExtensionTableData >extable.h)

isotope.h:	bin2hex.pl isotope-small.txt
		(cd $(srcdir); ./bin2hex.pl isotope-small.txt IsotopeData ISOTOPE >isotope.h)

phmodeldata.h:	bin2hex.pl phmodel.txt
		(cd $(srcdir); ./bin2hex.pl phmodel.txt PhModelData PHMODELDATA >phmodeldata.h)

resdata.h:	bin2hex.pl resdata.txt
		(cd $(srcdir); ./bin2hex.pl resdata.txt ResidueData >resdata.h)

types.h:	bin2hex.pl types.txt
		(cd $(srcdir); ./bin2hex.pl types.txt TypesData >types.h)
