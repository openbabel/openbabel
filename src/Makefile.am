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
		aromatic.h atomtyp.h base.h binary.h bitvec.h bondtyp.h \
		bondtyper.h chains.h chiral.h data.h element.h extable.h \
		fileformat.h generic.h grid.h isotope.h matrix.h mol.h \
		molchrg.h molvector.h oberror.h obifstream.h obutil.h \
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
	newmain.cpp bondtyper.cpp bondtyper.h \
	strncasecmp.c snprintf.c snprintf.h
