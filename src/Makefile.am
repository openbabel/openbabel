SUBDIRS = math windows

lib_LTLIBRARIES = libopenbabel.la

libopenbabel_la_SOURCES = \
		alchemy.cpp amber.cpp atom.cpp balst.cpp base.cpp bgf.cpp \
		binary.cpp bitvec.cpp bond.cpp box.cpp c3d.cpp cacao.cpp \
		cache.cpp car.cpp ccc.cpp chains.cpp chdrw.cpp chemtool.cpp chiral.cpp \
		cml.cpp csr.cpp cssr.cpp data.cpp dmol.cpp feat.cpp fh.cpp \
		fileformat.cpp gamess.cpp gaussian.cpp generic.cpp \
		ghemical.cpp grid.cpp gromos96.cpp hin.cpp jaguar.cpp \
		matrix.cpp mdl.cpp mm3.cpp mmod.cpp mol.cpp mol2.cpp \
		molchrg.cpp molvector.cpp mopac.cpp mpqc.cpp nwchem.cpp \
		oberror.cpp obutil.cpp parsmart.cpp parsmi.cpp patty.cpp \
		pdb.cpp phmodel.cpp povray.cpp qchem.cpp rand.cpp report.cpp \
		residue.cpp ring.cpp rotor.cpp shelx.cpp smi.cpp tinker.cpp \
		tokenst.cpp typer.cpp unichem.cpp viewmol.cpp xed.cpp xyz.cpp \
		zindo.cpp

libopenbabel_la_LIBADD = math/libmath.la @LTLIBOBJS@

libopenbabelincludedir = $(includedir)/openbabel
libopenbabelinclude_HEADERS = \
		aromatic.h atomtyp.h babelconfig.h base.h binary.h bitvec.h \
		bondtyp.h bondtyper.h chains.h chiral.h data.h element.h \
		extable.h fileformat.h generic.h grid.h isotope.h matrix.h \
		mol.h molchrg.h molvector.h oberror.h obifstream.h obutil.h \
		parsmart.h patty.h phmodel.h phmodeldata.h resdata.h ring.h \
		rotor.h smi.h typer.h types.h

bin_PROGRAMS = babel
babel_SOURCES = main.cpp
babel_LDADD = libopenbabel.la @LIBOBJS@

pkgdata_DATA = \
		aromatic.txt bondtyp.txt extable.txt isotope.txt resdata.txt \
		atomtyp.txt element.txt isotope-small.txt phmodel.txt types.txt

EXTRA_DIST = \
	bc32.mak borland.bat oblib.dsp oblib.dsw \
	bin2hex.pl $(pkgdata_DATA) \
	strncasecmp.c snprintf.c snprintf.h

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
