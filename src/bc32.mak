###################################################################
#
#       Borland specific directives ---
#
.SWAP
.AUTODEPEND

openbabel: libopenbabel.lib main.obj
	ilink32 -Tpe -ap c0x32 main.obj  ,babel,, import32.lib cw32mt.lib libopenbabel.lib

#left out: main.obj  newmain.obj povray.obj
libopenbabel.lib: alchemy.obj amber.obj atom.obj balst.obj base.obj bgf.obj binary.obj bitgrid.obj bitvec.obj bond.obj box.obj c3d.obj cacao.obj cache.obj car.obj ccc.obj chains.obj chdrw.obj  chiral.obj cml.obj csr.obj cssr.obj cwrap.obj data.obj dmol.obj feat.obj fh.obj fileformat.obj gamess.obj gaussian.obj generic.obj ghemical.obj grid.obj gromos96.obj hin.obj jaguar.obj matrix.obj mdl.obj mm3.obj mmod.obj mol.obj mol2.obj molchrg.obj molvector.obj mopac.obj mpqc.obj nwchem.obj oberror.obj obutil.obj parsmart.obj parsmi.obj patty.obj pdb.obj phmodel.obj qchem.obj report.obj residue.obj ring.obj rotor.obj smi.obj tinker.obj tokenst.obj typer.obj unichem.obj viewmol.obj xed.obj xyz.obj zindo.obj math\matrix3x3.obj math\vector3.obj rand.obj
	tlib "libopenbabel.lib" -+alchemy.obj -+amber.obj -+atom.obj -+balst.obj -+base.obj -+bgf.obj -+binary.obj -+bitgrid.obj -+bitvec.obj -+bond.obj -+box.obj -+c3d.obj -+cacao.obj -+cache.obj -+car.obj -+ccc.obj -+chains.obj -+chdrw.obj  -+chiral.obj -+cml.obj -+csr.obj -+cssr.obj -+cwrap.obj -+data.obj -+dmol.obj -+feat.obj -+fh.obj -+fileformat.obj -+gamess.obj -+gaussian.obj -+generic.obj -+ghemical.obj -+grid.obj -+gromos96.obj -+hin.obj -+jaguar.obj -+matrix.obj -+mdl.obj -+mm3.obj -+mmod.obj -+mol.obj -+mol2.obj -+molchrg.obj -+molvector.obj -+mopac.obj -+mpqc.obj -+nwchem.obj -+oberror.obj -+obutil.obj -+parsmart.obj -+parsmi.obj -+patty.obj -+pdb.obj -+phmodel.obj -+qchem.obj -+report.obj -+residue.obj -+ring.obj -+rotor.obj -+smi.obj -+tinker.obj -+tokenst.obj -+typer.obj -+unichem.obj -+viewmol.obj -+xed.obj -+xyz.obj -+zindo.obj -+matrix3x3.obj -+vector3.obj -+rand.obj


lib: libopenbabel.lib

.c.obj:
	bcc32 -c $*.c

.cpp.obj:
	bcc32 -tWC -DWIN32 -DDATADIR=\"C:\dev\openbabel\" -c $*.cpp

tidy:
	del /F *.obj math\*.obj
   
clean: tidy
	del  /F *.exe *.lib

default: openbabel
