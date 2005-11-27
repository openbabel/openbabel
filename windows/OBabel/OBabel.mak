# Microsoft Developer Studio Generated NMAKE File, Based on OBabel.dsp
!IF "$(CFG)" == ""
CFG=OBabel - Win32 Debug
!MESSAGE No configuration specified. Defaulting to OBabel - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "OBabel - Win32 Release" && "$(CFG)" != "OBabel - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "OBabel.mak" CFG="OBabel - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "OBabel - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "OBabel - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "OBabel - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release

ALL : ".\babel.exe"


CLEAN :
	-@erase "$(INTDIR)\alchemyformat.obj"
	-@erase "$(INTDIR)\amberformat.obj"
	-@erase "$(INTDIR)\APIInterface.obj"
	-@erase "$(INTDIR)\atom.obj"
	-@erase "$(INTDIR)\balstformat.obj"
	-@erase "$(INTDIR)\base.obj"
	-@erase "$(INTDIR)\bgfformat.obj"
	-@erase "$(INTDIR)\bitvec.obj"
	-@erase "$(INTDIR)\bond.obj"
	-@erase "$(INTDIR)\bondtyper.obj"
	-@erase "$(INTDIR)\boxformat.obj"
	-@erase "$(INTDIR)\cacaoformat.obj"
	-@erase "$(INTDIR)\cacheformat.obj"
	-@erase "$(INTDIR)\carformat.obj"
	-@erase "$(INTDIR)\cccformat.obj"
	-@erase "$(INTDIR)\chains.obj"
	-@erase "$(INTDIR)\chem3dformat.obj"
	-@erase "$(INTDIR)\chemdrawformat.obj"
	-@erase "$(INTDIR)\chemtoolformat.obj"
	-@erase "$(INTDIR)\chiral.obj"
	-@erase "$(INTDIR)\cmlreactlformat.obj"
	-@erase "$(INTDIR)\copyformat.obj"
	-@erase "$(INTDIR)\crkformat.obj"
	-@erase "$(INTDIR)\CSRformat.obj"
	-@erase "$(INTDIR)\cssrformat.obj"
	-@erase "$(INTDIR)\data.obj"
	-@erase "$(INTDIR)\dlhandler_win32.obj"
	-@erase "$(INTDIR)\dmolformat.obj"
	-@erase "$(INTDIR)\fastsearchformat.obj"
	-@erase "$(INTDIR)\featformat.obj"
	-@erase "$(INTDIR)\fhformat.obj"
	-@erase "$(INTDIR)\finger2.obj"
	-@erase "$(INTDIR)\finger3.obj"
	-@erase "$(INTDIR)\fingerprint.obj"
	-@erase "$(INTDIR)\fingerprintformat.obj"
	-@erase "$(INTDIR)\freefracformat.obj"
	-@erase "$(INTDIR)\gamessformat.obj"
	-@erase "$(INTDIR)\gaussformat.obj"
	-@erase "$(INTDIR)\generic.obj"
	-@erase "$(INTDIR)\ghemicalformat.obj"
	-@erase "$(INTDIR)\grid.obj"
	-@erase "$(INTDIR)\gromos96format.obj"
	-@erase "$(INTDIR)\hinformat.obj"
	-@erase "$(INTDIR)\inchiformat.obj"
	-@erase "$(INTDIR)\jaguarformat.obj"
	-@erase "$(INTDIR)\kekulize.obj"
	-@erase "$(INTDIR)\main.obj"
	-@erase "$(INTDIR)\matrix.obj"
	-@erase "$(INTDIR)\matrix3x3.obj"
	-@erase "$(INTDIR)\mdlformat.obj"
	-@erase "$(INTDIR)\mmodformat.obj"
	-@erase "$(INTDIR)\mol.obj"
	-@erase "$(INTDIR)\mol2format.obj"
	-@erase "$(INTDIR)\molchrg.obj"
	-@erase "$(INTDIR)\mopacformat.obj"
	-@erase "$(INTDIR)\mpdformat.obj"
	-@erase "$(INTDIR)\mpqcformat.obj"
	-@erase "$(INTDIR)\nwchemformat.obj"
	-@erase "$(INTDIR)\obconversion.obj"
	-@erase "$(INTDIR)\oberror.obj"
	-@erase "$(INTDIR)\obiter.obj"
	-@erase "$(INTDIR)\obutil.obj"
	-@erase "$(INTDIR)\parsmart.obj"
	-@erase "$(INTDIR)\patty.obj"
	-@erase "$(INTDIR)\pcmodelformat.obj"
	-@erase "$(INTDIR)\pdbformat.obj"
	-@erase "$(INTDIR)\phmodel.obj"
	-@erase "$(INTDIR)\povrayformat.obj"
	-@erase "$(INTDIR)\PQSformat.obj"
	-@erase "$(INTDIR)\pubchem.obj"
	-@erase "$(INTDIR)\qchemformat.obj"
	-@erase "$(INTDIR)\rand.obj"
	-@erase "$(INTDIR)\reportformat.obj"
	-@erase "$(INTDIR)\residue.obj"
	-@erase "$(INTDIR)\ring.obj"
	-@erase "$(INTDIR)\rotamer.obj"
	-@erase "$(INTDIR)\rotor.obj"
	-@erase "$(INTDIR)\rxnformat.obj"
	-@erase "$(INTDIR)\shelxformat.obj"
	-@erase "$(INTDIR)\smilesformat.obj"
	-@erase "$(INTDIR)\tinkerformat.obj"
	-@erase "$(INTDIR)\tokenst.obj"
	-@erase "$(INTDIR)\transform.obj"
	-@erase "$(INTDIR)\turbomoleformat.obj"
	-@erase "$(INTDIR)\typer.obj"
	-@erase "$(INTDIR)\unichemformat.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vector3.obj"
	-@erase "$(INTDIR)\viewmolformat.obj"
	-@erase "$(INTDIR)\xcmlformat.obj"
	-@erase "$(INTDIR)\xedformat.obj"
	-@erase "$(INTDIR)\xml.obj"
	-@erase "$(INTDIR)\xmlformat.obj"
	-@erase "$(INTDIR)\xyzformat.obj"
	-@erase "$(INTDIR)\yasaraformat.obj"
	-@erase "$(INTDIR)\zindoformat.obj"
	-@erase ".\babel.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP_PROJ=/nologo /MT /W3 /GR /GX /O1 /I "..\..\src" /I ".." /I "../../data" /I "..\..\src\formats" /I "..\..\src\formats\xml" /D "NDEBUG" /D "WIN32" /D "_CONSOLE" /D "_MBCS" /D "INCHI_LINK_AS_DLL" /D "HAVE_CONFIG_H" /Fp"$(INTDIR)\OBabel.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\OBabel.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=libinchi.lib /nologo /subsystem:console /incremental:no /pdb:"$(OUTDIR)\babel.pdb" /machine:I386 /out:"babel.exe" /libpath:".." 
LINK32_OBJS= \
	"$(INTDIR)\alchemyformat.obj" \
	"$(INTDIR)\amberformat.obj" \
	"$(INTDIR)\APIInterface.obj" \
	"$(INTDIR)\atom.obj" \
	"$(INTDIR)\balstformat.obj" \
	"$(INTDIR)\base.obj" \
	"$(INTDIR)\bgfformat.obj" \
	"$(INTDIR)\bitvec.obj" \
	"$(INTDIR)\bond.obj" \
	"$(INTDIR)\bondtyper.obj" \
	"$(INTDIR)\boxformat.obj" \
	"$(INTDIR)\cacaoformat.obj" \
	"$(INTDIR)\cacheformat.obj" \
	"$(INTDIR)\carformat.obj" \
	"$(INTDIR)\cccformat.obj" \
	"$(INTDIR)\chains.obj" \
	"$(INTDIR)\chem3dformat.obj" \
	"$(INTDIR)\chemdrawformat.obj" \
	"$(INTDIR)\chemtoolformat.obj" \
	"$(INTDIR)\chiral.obj" \
	"$(INTDIR)\cmlreactlformat.obj" \
	"$(INTDIR)\copyformat.obj" \
	"$(INTDIR)\crkformat.obj" \
	"$(INTDIR)\CSRformat.obj" \
	"$(INTDIR)\cssrformat.obj" \
	"$(INTDIR)\data.obj" \
	"$(INTDIR)\dlhandler_win32.obj" \
	"$(INTDIR)\dmolformat.obj" \
	"$(INTDIR)\fastsearchformat.obj" \
	"$(INTDIR)\featformat.obj" \
	"$(INTDIR)\fhformat.obj" \
	"$(INTDIR)\finger2.obj" \
	"$(INTDIR)\finger3.obj" \
	"$(INTDIR)\fingerprint.obj" \
	"$(INTDIR)\fingerprintformat.obj" \
	"$(INTDIR)\freefracformat.obj" \
	"$(INTDIR)\gamessformat.obj" \
	"$(INTDIR)\gaussformat.obj" \
	"$(INTDIR)\generic.obj" \
	"$(INTDIR)\ghemicalformat.obj" \
	"$(INTDIR)\grid.obj" \
	"$(INTDIR)\gromos96format.obj" \
	"$(INTDIR)\hinformat.obj" \
	"$(INTDIR)\inchiformat.obj" \
	"$(INTDIR)\jaguarformat.obj" \
	"$(INTDIR)\kekulize.obj" \
	"$(INTDIR)\main.obj" \
	"$(INTDIR)\matrix.obj" \
	"$(INTDIR)\matrix3x3.obj" \
	"$(INTDIR)\mdlformat.obj" \
	"$(INTDIR)\mmodformat.obj" \
	"$(INTDIR)\mol.obj" \
	"$(INTDIR)\mol2format.obj" \
	"$(INTDIR)\molchrg.obj" \
	"$(INTDIR)\mopacformat.obj" \
	"$(INTDIR)\mpdformat.obj" \
	"$(INTDIR)\mpqcformat.obj" \
	"$(INTDIR)\nwchemformat.obj" \
	"$(INTDIR)\obconversion.obj" \
	"$(INTDIR)\oberror.obj" \
	"$(INTDIR)\obiter.obj" \
	"$(INTDIR)\obutil.obj" \
	"$(INTDIR)\parsmart.obj" \
	"$(INTDIR)\patty.obj" \
	"$(INTDIR)\pcmodelformat.obj" \
	"$(INTDIR)\pdbformat.obj" \
	"$(INTDIR)\phmodel.obj" \
	"$(INTDIR)\povrayformat.obj" \
	"$(INTDIR)\PQSformat.obj" \
	"$(INTDIR)\pubchem.obj" \
	"$(INTDIR)\qchemformat.obj" \
	"$(INTDIR)\rand.obj" \
	"$(INTDIR)\reportformat.obj" \
	"$(INTDIR)\residue.obj" \
	"$(INTDIR)\ring.obj" \
	"$(INTDIR)\rotamer.obj" \
	"$(INTDIR)\rotor.obj" \
	"$(INTDIR)\rxnformat.obj" \
	"$(INTDIR)\shelxformat.obj" \
	"$(INTDIR)\smilesformat.obj" \
	"$(INTDIR)\tinkerformat.obj" \
	"$(INTDIR)\tokenst.obj" \
	"$(INTDIR)\transform.obj" \
	"$(INTDIR)\turbomoleformat.obj" \
	"$(INTDIR)\typer.obj" \
	"$(INTDIR)\unichemformat.obj" \
	"$(INTDIR)\vector3.obj" \
	"$(INTDIR)\viewmolformat.obj" \
	"$(INTDIR)\xcmlformat.obj" \
	"$(INTDIR)\xedformat.obj" \
	"$(INTDIR)\xml.obj" \
	"$(INTDIR)\xmlformat.obj" \
	"$(INTDIR)\xyzformat.obj" \
	"$(INTDIR)\yasaraformat.obj" \
	"$(INTDIR)\zindoformat.obj" \
	"..\libxml2.lib" \
	"..\zdll.lib"

".\babel.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"

OUTDIR=.\Debug
INTDIR=.\Debug
# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

ALL : "$(OUTDIR)\babel.exe" "$(OUTDIR)\OBabel.bsc"


CLEAN :
	-@erase "$(INTDIR)\alchemyformat.obj"
	-@erase "$(INTDIR)\alchemyformat.sbr"
	-@erase "$(INTDIR)\amberformat.obj"
	-@erase "$(INTDIR)\amberformat.sbr"
	-@erase "$(INTDIR)\APIInterface.obj"
	-@erase "$(INTDIR)\APIInterface.sbr"
	-@erase "$(INTDIR)\atom.obj"
	-@erase "$(INTDIR)\atom.sbr"
	-@erase "$(INTDIR)\balstformat.obj"
	-@erase "$(INTDIR)\balstformat.sbr"
	-@erase "$(INTDIR)\base.obj"
	-@erase "$(INTDIR)\base.sbr"
	-@erase "$(INTDIR)\bgfformat.obj"
	-@erase "$(INTDIR)\bgfformat.sbr"
	-@erase "$(INTDIR)\bitvec.obj"
	-@erase "$(INTDIR)\bitvec.sbr"
	-@erase "$(INTDIR)\bond.obj"
	-@erase "$(INTDIR)\bond.sbr"
	-@erase "$(INTDIR)\bondtyper.obj"
	-@erase "$(INTDIR)\bondtyper.sbr"
	-@erase "$(INTDIR)\boxformat.obj"
	-@erase "$(INTDIR)\boxformat.sbr"
	-@erase "$(INTDIR)\cacaoformat.obj"
	-@erase "$(INTDIR)\cacaoformat.sbr"
	-@erase "$(INTDIR)\cacheformat.obj"
	-@erase "$(INTDIR)\cacheformat.sbr"
	-@erase "$(INTDIR)\carformat.obj"
	-@erase "$(INTDIR)\carformat.sbr"
	-@erase "$(INTDIR)\cccformat.obj"
	-@erase "$(INTDIR)\cccformat.sbr"
	-@erase "$(INTDIR)\chains.obj"
	-@erase "$(INTDIR)\chains.sbr"
	-@erase "$(INTDIR)\chem3dformat.obj"
	-@erase "$(INTDIR)\chem3dformat.sbr"
	-@erase "$(INTDIR)\chemdrawformat.obj"
	-@erase "$(INTDIR)\chemdrawformat.sbr"
	-@erase "$(INTDIR)\chemtoolformat.obj"
	-@erase "$(INTDIR)\chemtoolformat.sbr"
	-@erase "$(INTDIR)\chiral.obj"
	-@erase "$(INTDIR)\chiral.sbr"
	-@erase "$(INTDIR)\cmlreactlformat.obj"
	-@erase "$(INTDIR)\cmlreactlformat.sbr"
	-@erase "$(INTDIR)\copyformat.obj"
	-@erase "$(INTDIR)\copyformat.sbr"
	-@erase "$(INTDIR)\crkformat.obj"
	-@erase "$(INTDIR)\crkformat.sbr"
	-@erase "$(INTDIR)\CSRformat.obj"
	-@erase "$(INTDIR)\CSRformat.sbr"
	-@erase "$(INTDIR)\cssrformat.obj"
	-@erase "$(INTDIR)\cssrformat.sbr"
	-@erase "$(INTDIR)\data.obj"
	-@erase "$(INTDIR)\data.sbr"
	-@erase "$(INTDIR)\dlhandler_win32.obj"
	-@erase "$(INTDIR)\dlhandler_win32.sbr"
	-@erase "$(INTDIR)\dmolformat.obj"
	-@erase "$(INTDIR)\dmolformat.sbr"
	-@erase "$(INTDIR)\fastsearchformat.obj"
	-@erase "$(INTDIR)\fastsearchformat.sbr"
	-@erase "$(INTDIR)\featformat.obj"
	-@erase "$(INTDIR)\featformat.sbr"
	-@erase "$(INTDIR)\fhformat.obj"
	-@erase "$(INTDIR)\fhformat.sbr"
	-@erase "$(INTDIR)\finger2.obj"
	-@erase "$(INTDIR)\finger2.sbr"
	-@erase "$(INTDIR)\finger3.obj"
	-@erase "$(INTDIR)\finger3.sbr"
	-@erase "$(INTDIR)\fingerprint.obj"
	-@erase "$(INTDIR)\fingerprint.sbr"
	-@erase "$(INTDIR)\fingerprintformat.obj"
	-@erase "$(INTDIR)\fingerprintformat.sbr"
	-@erase "$(INTDIR)\freefracformat.obj"
	-@erase "$(INTDIR)\freefracformat.sbr"
	-@erase "$(INTDIR)\gamessformat.obj"
	-@erase "$(INTDIR)\gamessformat.sbr"
	-@erase "$(INTDIR)\gaussformat.obj"
	-@erase "$(INTDIR)\gaussformat.sbr"
	-@erase "$(INTDIR)\generic.obj"
	-@erase "$(INTDIR)\generic.sbr"
	-@erase "$(INTDIR)\ghemicalformat.obj"
	-@erase "$(INTDIR)\ghemicalformat.sbr"
	-@erase "$(INTDIR)\grid.obj"
	-@erase "$(INTDIR)\grid.sbr"
	-@erase "$(INTDIR)\gromos96format.obj"
	-@erase "$(INTDIR)\gromos96format.sbr"
	-@erase "$(INTDIR)\hinformat.obj"
	-@erase "$(INTDIR)\hinformat.sbr"
	-@erase "$(INTDIR)\inchiformat.obj"
	-@erase "$(INTDIR)\inchiformat.sbr"
	-@erase "$(INTDIR)\jaguarformat.obj"
	-@erase "$(INTDIR)\jaguarformat.sbr"
	-@erase "$(INTDIR)\kekulize.obj"
	-@erase "$(INTDIR)\kekulize.sbr"
	-@erase "$(INTDIR)\main.obj"
	-@erase "$(INTDIR)\main.sbr"
	-@erase "$(INTDIR)\matrix.obj"
	-@erase "$(INTDIR)\matrix.sbr"
	-@erase "$(INTDIR)\matrix3x3.obj"
	-@erase "$(INTDIR)\matrix3x3.sbr"
	-@erase "$(INTDIR)\mdlformat.obj"
	-@erase "$(INTDIR)\mdlformat.sbr"
	-@erase "$(INTDIR)\mmodformat.obj"
	-@erase "$(INTDIR)\mmodformat.sbr"
	-@erase "$(INTDIR)\mol.obj"
	-@erase "$(INTDIR)\mol.sbr"
	-@erase "$(INTDIR)\mol2format.obj"
	-@erase "$(INTDIR)\mol2format.sbr"
	-@erase "$(INTDIR)\molchrg.obj"
	-@erase "$(INTDIR)\molchrg.sbr"
	-@erase "$(INTDIR)\mopacformat.obj"
	-@erase "$(INTDIR)\mopacformat.sbr"
	-@erase "$(INTDIR)\mpdformat.obj"
	-@erase "$(INTDIR)\mpdformat.sbr"
	-@erase "$(INTDIR)\mpqcformat.obj"
	-@erase "$(INTDIR)\mpqcformat.sbr"
	-@erase "$(INTDIR)\nwchemformat.obj"
	-@erase "$(INTDIR)\nwchemformat.sbr"
	-@erase "$(INTDIR)\obconversion.obj"
	-@erase "$(INTDIR)\obconversion.sbr"
	-@erase "$(INTDIR)\oberror.obj"
	-@erase "$(INTDIR)\oberror.sbr"
	-@erase "$(INTDIR)\obiter.obj"
	-@erase "$(INTDIR)\obiter.sbr"
	-@erase "$(INTDIR)\obutil.obj"
	-@erase "$(INTDIR)\obutil.sbr"
	-@erase "$(INTDIR)\parsmart.obj"
	-@erase "$(INTDIR)\parsmart.sbr"
	-@erase "$(INTDIR)\patty.obj"
	-@erase "$(INTDIR)\patty.sbr"
	-@erase "$(INTDIR)\pcmodelformat.obj"
	-@erase "$(INTDIR)\pcmodelformat.sbr"
	-@erase "$(INTDIR)\pdbformat.obj"
	-@erase "$(INTDIR)\pdbformat.sbr"
	-@erase "$(INTDIR)\phmodel.obj"
	-@erase "$(INTDIR)\phmodel.sbr"
	-@erase "$(INTDIR)\povrayformat.obj"
	-@erase "$(INTDIR)\povrayformat.sbr"
	-@erase "$(INTDIR)\PQSformat.obj"
	-@erase "$(INTDIR)\PQSformat.sbr"
	-@erase "$(INTDIR)\pubchem.obj"
	-@erase "$(INTDIR)\pubchem.sbr"
	-@erase "$(INTDIR)\qchemformat.obj"
	-@erase "$(INTDIR)\qchemformat.sbr"
	-@erase "$(INTDIR)\rand.obj"
	-@erase "$(INTDIR)\rand.sbr"
	-@erase "$(INTDIR)\reportformat.obj"
	-@erase "$(INTDIR)\reportformat.sbr"
	-@erase "$(INTDIR)\residue.obj"
	-@erase "$(INTDIR)\residue.sbr"
	-@erase "$(INTDIR)\ring.obj"
	-@erase "$(INTDIR)\ring.sbr"
	-@erase "$(INTDIR)\rotamer.obj"
	-@erase "$(INTDIR)\rotamer.sbr"
	-@erase "$(INTDIR)\rotor.obj"
	-@erase "$(INTDIR)\rotor.sbr"
	-@erase "$(INTDIR)\rxnformat.obj"
	-@erase "$(INTDIR)\rxnformat.sbr"
	-@erase "$(INTDIR)\shelxformat.obj"
	-@erase "$(INTDIR)\shelxformat.sbr"
	-@erase "$(INTDIR)\smilesformat.obj"
	-@erase "$(INTDIR)\smilesformat.sbr"
	-@erase "$(INTDIR)\tinkerformat.obj"
	-@erase "$(INTDIR)\tinkerformat.sbr"
	-@erase "$(INTDIR)\tokenst.obj"
	-@erase "$(INTDIR)\tokenst.sbr"
	-@erase "$(INTDIR)\transform.obj"
	-@erase "$(INTDIR)\transform.sbr"
	-@erase "$(INTDIR)\turbomoleformat.obj"
	-@erase "$(INTDIR)\turbomoleformat.sbr"
	-@erase "$(INTDIR)\typer.obj"
	-@erase "$(INTDIR)\typer.sbr"
	-@erase "$(INTDIR)\unichemformat.obj"
	-@erase "$(INTDIR)\unichemformat.sbr"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vc60.pdb"
	-@erase "$(INTDIR)\vector3.obj"
	-@erase "$(INTDIR)\vector3.sbr"
	-@erase "$(INTDIR)\viewmolformat.obj"
	-@erase "$(INTDIR)\viewmolformat.sbr"
	-@erase "$(INTDIR)\xcmlformat.obj"
	-@erase "$(INTDIR)\xcmlformat.sbr"
	-@erase "$(INTDIR)\xedformat.obj"
	-@erase "$(INTDIR)\xedformat.sbr"
	-@erase "$(INTDIR)\xml.obj"
	-@erase "$(INTDIR)\xml.sbr"
	-@erase "$(INTDIR)\xmlformat.obj"
	-@erase "$(INTDIR)\xmlformat.sbr"
	-@erase "$(INTDIR)\xyzformat.obj"
	-@erase "$(INTDIR)\xyzformat.sbr"
	-@erase "$(INTDIR)\yasaraformat.obj"
	-@erase "$(INTDIR)\yasaraformat.sbr"
	-@erase "$(INTDIR)\zindoformat.obj"
	-@erase "$(INTDIR)\zindoformat.sbr"
	-@erase "$(OUTDIR)\babel.exe"
	-@erase "$(OUTDIR)\babel.ilk"
	-@erase "$(OUTDIR)\babel.pdb"
	-@erase "$(OUTDIR)\OBabel.bsc"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP_PROJ=/nologo /MTd /W3 /Gm /GR /GX /ZI /Od /I "..\..\src" /I ".." /I "../../data" /I "..\..\src\formats" /I "..\..\src\formats\xml" /D "_DEBUG" /D "WIN32" /D "_CONSOLE" /D "_MBCS" /D "INCHI_LINK_AS_DLL" /D "HAVE_CONFIG_H" /FR"$(INTDIR)\\" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\OBabel.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\alchemyformat.sbr" \
	"$(INTDIR)\amberformat.sbr" \
	"$(INTDIR)\APIInterface.sbr" \
	"$(INTDIR)\atom.sbr" \
	"$(INTDIR)\balstformat.sbr" \
	"$(INTDIR)\base.sbr" \
	"$(INTDIR)\bgfformat.sbr" \
	"$(INTDIR)\bitvec.sbr" \
	"$(INTDIR)\bond.sbr" \
	"$(INTDIR)\bondtyper.sbr" \
	"$(INTDIR)\boxformat.sbr" \
	"$(INTDIR)\cacaoformat.sbr" \
	"$(INTDIR)\cacheformat.sbr" \
	"$(INTDIR)\carformat.sbr" \
	"$(INTDIR)\cccformat.sbr" \
	"$(INTDIR)\chains.sbr" \
	"$(INTDIR)\chem3dformat.sbr" \
	"$(INTDIR)\chemdrawformat.sbr" \
	"$(INTDIR)\chemtoolformat.sbr" \
	"$(INTDIR)\chiral.sbr" \
	"$(INTDIR)\cmlreactlformat.sbr" \
	"$(INTDIR)\copyformat.sbr" \
	"$(INTDIR)\crkformat.sbr" \
	"$(INTDIR)\CSRformat.sbr" \
	"$(INTDIR)\cssrformat.sbr" \
	"$(INTDIR)\data.sbr" \
	"$(INTDIR)\dlhandler_win32.sbr" \
	"$(INTDIR)\dmolformat.sbr" \
	"$(INTDIR)\fastsearchformat.sbr" \
	"$(INTDIR)\featformat.sbr" \
	"$(INTDIR)\fhformat.sbr" \
	"$(INTDIR)\finger2.sbr" \
	"$(INTDIR)\finger3.sbr" \
	"$(INTDIR)\fingerprint.sbr" \
	"$(INTDIR)\fingerprintformat.sbr" \
	"$(INTDIR)\freefracformat.sbr" \
	"$(INTDIR)\gamessformat.sbr" \
	"$(INTDIR)\gaussformat.sbr" \
	"$(INTDIR)\generic.sbr" \
	"$(INTDIR)\ghemicalformat.sbr" \
	"$(INTDIR)\grid.sbr" \
	"$(INTDIR)\gromos96format.sbr" \
	"$(INTDIR)\hinformat.sbr" \
	"$(INTDIR)\inchiformat.sbr" \
	"$(INTDIR)\jaguarformat.sbr" \
	"$(INTDIR)\kekulize.sbr" \
	"$(INTDIR)\main.sbr" \
	"$(INTDIR)\matrix.sbr" \
	"$(INTDIR)\matrix3x3.sbr" \
	"$(INTDIR)\mdlformat.sbr" \
	"$(INTDIR)\mmodformat.sbr" \
	"$(INTDIR)\mol.sbr" \
	"$(INTDIR)\mol2format.sbr" \
	"$(INTDIR)\molchrg.sbr" \
	"$(INTDIR)\mopacformat.sbr" \
	"$(INTDIR)\mpdformat.sbr" \
	"$(INTDIR)\mpqcformat.sbr" \
	"$(INTDIR)\nwchemformat.sbr" \
	"$(INTDIR)\obconversion.sbr" \
	"$(INTDIR)\oberror.sbr" \
	"$(INTDIR)\obiter.sbr" \
	"$(INTDIR)\obutil.sbr" \
	"$(INTDIR)\parsmart.sbr" \
	"$(INTDIR)\patty.sbr" \
	"$(INTDIR)\pcmodelformat.sbr" \
	"$(INTDIR)\pdbformat.sbr" \
	"$(INTDIR)\phmodel.sbr" \
	"$(INTDIR)\povrayformat.sbr" \
	"$(INTDIR)\PQSformat.sbr" \
	"$(INTDIR)\pubchem.sbr" \
	"$(INTDIR)\qchemformat.sbr" \
	"$(INTDIR)\rand.sbr" \
	"$(INTDIR)\reportformat.sbr" \
	"$(INTDIR)\residue.sbr" \
	"$(INTDIR)\ring.sbr" \
	"$(INTDIR)\rotamer.sbr" \
	"$(INTDIR)\rotor.sbr" \
	"$(INTDIR)\rxnformat.sbr" \
	"$(INTDIR)\shelxformat.sbr" \
	"$(INTDIR)\smilesformat.sbr" \
	"$(INTDIR)\tinkerformat.sbr" \
	"$(INTDIR)\tokenst.sbr" \
	"$(INTDIR)\transform.sbr" \
	"$(INTDIR)\turbomoleformat.sbr" \
	"$(INTDIR)\typer.sbr" \
	"$(INTDIR)\unichemformat.sbr" \
	"$(INTDIR)\vector3.sbr" \
	"$(INTDIR)\viewmolformat.sbr" \
	"$(INTDIR)\xcmlformat.sbr" \
	"$(INTDIR)\xedformat.sbr" \
	"$(INTDIR)\xml.sbr" \
	"$(INTDIR)\xmlformat.sbr" \
	"$(INTDIR)\xyzformat.sbr" \
	"$(INTDIR)\yasaraformat.sbr" \
	"$(INTDIR)\zindoformat.sbr"

"$(OUTDIR)\OBabel.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib libinchi.lib /nologo /subsystem:console /incremental:yes /pdb:"$(OUTDIR)\babel.pdb" /debug /machine:I386 /out:"$(OUTDIR)\babel.exe" /pdbtype:sept /libpath:".." 
LINK32_OBJS= \
	"$(INTDIR)\alchemyformat.obj" \
	"$(INTDIR)\amberformat.obj" \
	"$(INTDIR)\APIInterface.obj" \
	"$(INTDIR)\atom.obj" \
	"$(INTDIR)\balstformat.obj" \
	"$(INTDIR)\base.obj" \
	"$(INTDIR)\bgfformat.obj" \
	"$(INTDIR)\bitvec.obj" \
	"$(INTDIR)\bond.obj" \
	"$(INTDIR)\bondtyper.obj" \
	"$(INTDIR)\boxformat.obj" \
	"$(INTDIR)\cacaoformat.obj" \
	"$(INTDIR)\cacheformat.obj" \
	"$(INTDIR)\carformat.obj" \
	"$(INTDIR)\cccformat.obj" \
	"$(INTDIR)\chains.obj" \
	"$(INTDIR)\chem3dformat.obj" \
	"$(INTDIR)\chemdrawformat.obj" \
	"$(INTDIR)\chemtoolformat.obj" \
	"$(INTDIR)\chiral.obj" \
	"$(INTDIR)\cmlreactlformat.obj" \
	"$(INTDIR)\copyformat.obj" \
	"$(INTDIR)\crkformat.obj" \
	"$(INTDIR)\CSRformat.obj" \
	"$(INTDIR)\cssrformat.obj" \
	"$(INTDIR)\data.obj" \
	"$(INTDIR)\dlhandler_win32.obj" \
	"$(INTDIR)\dmolformat.obj" \
	"$(INTDIR)\fastsearchformat.obj" \
	"$(INTDIR)\featformat.obj" \
	"$(INTDIR)\fhformat.obj" \
	"$(INTDIR)\finger2.obj" \
	"$(INTDIR)\finger3.obj" \
	"$(INTDIR)\fingerprint.obj" \
	"$(INTDIR)\fingerprintformat.obj" \
	"$(INTDIR)\freefracformat.obj" \
	"$(INTDIR)\gamessformat.obj" \
	"$(INTDIR)\gaussformat.obj" \
	"$(INTDIR)\generic.obj" \
	"$(INTDIR)\ghemicalformat.obj" \
	"$(INTDIR)\grid.obj" \
	"$(INTDIR)\gromos96format.obj" \
	"$(INTDIR)\hinformat.obj" \
	"$(INTDIR)\inchiformat.obj" \
	"$(INTDIR)\jaguarformat.obj" \
	"$(INTDIR)\kekulize.obj" \
	"$(INTDIR)\main.obj" \
	"$(INTDIR)\matrix.obj" \
	"$(INTDIR)\matrix3x3.obj" \
	"$(INTDIR)\mdlformat.obj" \
	"$(INTDIR)\mmodformat.obj" \
	"$(INTDIR)\mol.obj" \
	"$(INTDIR)\mol2format.obj" \
	"$(INTDIR)\molchrg.obj" \
	"$(INTDIR)\mopacformat.obj" \
	"$(INTDIR)\mpdformat.obj" \
	"$(INTDIR)\mpqcformat.obj" \
	"$(INTDIR)\nwchemformat.obj" \
	"$(INTDIR)\obconversion.obj" \
	"$(INTDIR)\oberror.obj" \
	"$(INTDIR)\obiter.obj" \
	"$(INTDIR)\obutil.obj" \
	"$(INTDIR)\parsmart.obj" \
	"$(INTDIR)\patty.obj" \
	"$(INTDIR)\pcmodelformat.obj" \
	"$(INTDIR)\pdbformat.obj" \
	"$(INTDIR)\phmodel.obj" \
	"$(INTDIR)\povrayformat.obj" \
	"$(INTDIR)\PQSformat.obj" \
	"$(INTDIR)\pubchem.obj" \
	"$(INTDIR)\qchemformat.obj" \
	"$(INTDIR)\rand.obj" \
	"$(INTDIR)\reportformat.obj" \
	"$(INTDIR)\residue.obj" \
	"$(INTDIR)\ring.obj" \
	"$(INTDIR)\rotamer.obj" \
	"$(INTDIR)\rotor.obj" \
	"$(INTDIR)\rxnformat.obj" \
	"$(INTDIR)\shelxformat.obj" \
	"$(INTDIR)\smilesformat.obj" \
	"$(INTDIR)\tinkerformat.obj" \
	"$(INTDIR)\tokenst.obj" \
	"$(INTDIR)\transform.obj" \
	"$(INTDIR)\turbomoleformat.obj" \
	"$(INTDIR)\typer.obj" \
	"$(INTDIR)\unichemformat.obj" \
	"$(INTDIR)\vector3.obj" \
	"$(INTDIR)\viewmolformat.obj" \
	"$(INTDIR)\xcmlformat.obj" \
	"$(INTDIR)\xedformat.obj" \
	"$(INTDIR)\xml.obj" \
	"$(INTDIR)\xmlformat.obj" \
	"$(INTDIR)\xyzformat.obj" \
	"$(INTDIR)\yasaraformat.obj" \
	"$(INTDIR)\zindoformat.obj" \
	"..\libxml2.lib" \
	"..\zdll.lib"

"$(OUTDIR)\babel.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 

.c{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("OBabel.dep")
!INCLUDE "OBabel.dep"
!ELSE 
!MESSAGE Warning: cannot find "OBabel.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "OBabel - Win32 Release" || "$(CFG)" == "OBabel - Win32 Debug"
SOURCE=..\..\src\formats\alchemyformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\alchemyformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\alchemyformat.obj"	"$(INTDIR)\alchemyformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\amberformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\amberformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\amberformat.obj"	"$(INTDIR)\amberformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\APIInterface.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\APIInterface.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\APIInterface.obj"	"$(INTDIR)\APIInterface.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\atom.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\atom.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\atom.obj"	"$(INTDIR)\atom.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\balstformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\balstformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\balstformat.obj"	"$(INTDIR)\balstformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\base.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\base.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\base.obj"	"$(INTDIR)\base.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\bgfformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\bgfformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\bgfformat.obj"	"$(INTDIR)\bgfformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\bitvec.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\bitvec.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\bitvec.obj"	"$(INTDIR)\bitvec.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\bond.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\bond.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\bond.obj"	"$(INTDIR)\bond.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\bondtyper.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\bondtyper.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\bondtyper.obj"	"$(INTDIR)\bondtyper.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\boxformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\boxformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\boxformat.obj"	"$(INTDIR)\boxformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\cacaoformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\cacaoformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\cacaoformat.obj"	"$(INTDIR)\cacaoformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\cacheformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\cacheformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\cacheformat.obj"	"$(INTDIR)\cacheformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\carformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\carformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\carformat.obj"	"$(INTDIR)\carformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\cccformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\cccformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\cccformat.obj"	"$(INTDIR)\cccformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\chains.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\chains.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\chains.obj"	"$(INTDIR)\chains.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\chem3dformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\chem3dformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\chem3dformat.obj"	"$(INTDIR)\chem3dformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\chemdrawformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\chemdrawformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\chemdrawformat.obj"	"$(INTDIR)\chemdrawformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\chemtoolformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\chemtoolformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\chemtoolformat.obj"	"$(INTDIR)\chemtoolformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\chiral.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\chiral.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\chiral.obj"	"$(INTDIR)\chiral.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\xml\cmlreactlformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\cmlreactlformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\cmlreactlformat.obj"	"$(INTDIR)\cmlreactlformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\copyformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\copyformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\copyformat.obj"	"$(INTDIR)\copyformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\crkformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\crkformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\crkformat.obj"	"$(INTDIR)\crkformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\CSRformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\CSRformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\CSRformat.obj"	"$(INTDIR)\CSRformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\cssrformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\cssrformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\cssrformat.obj"	"$(INTDIR)\cssrformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\data.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\data.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\data.obj"	"$(INTDIR)\data.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\dlhandler_win32.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\dlhandler_win32.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\dlhandler_win32.obj"	"$(INTDIR)\dlhandler_win32.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\dmolformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\dmolformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\dmolformat.obj"	"$(INTDIR)\dmolformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\fastsearchformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\fastsearchformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\fastsearchformat.obj"	"$(INTDIR)\fastsearchformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\featformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\featformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\featformat.obj"	"$(INTDIR)\featformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\fhformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\fhformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\fhformat.obj"	"$(INTDIR)\fhformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\fingerprints\finger2.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\finger2.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\finger2.obj"	"$(INTDIR)\finger2.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\fingerprints\finger3.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\finger3.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\finger3.obj"	"$(INTDIR)\finger3.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\fingerprint.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\fingerprint.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\fingerprint.obj"	"$(INTDIR)\fingerprint.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\fingerprintformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\fingerprintformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\fingerprintformat.obj"	"$(INTDIR)\fingerprintformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\freefracformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\freefracformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\freefracformat.obj"	"$(INTDIR)\freefracformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\gamessformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\gamessformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\gamessformat.obj"	"$(INTDIR)\gamessformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\gaussformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\gaussformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\gaussformat.obj"	"$(INTDIR)\gaussformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\generic.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\generic.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\generic.obj"	"$(INTDIR)\generic.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\ghemicalformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\ghemicalformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\ghemicalformat.obj"	"$(INTDIR)\ghemicalformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\grid.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\grid.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\grid.obj"	"$(INTDIR)\grid.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\gromos96format.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\gromos96format.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\gromos96format.obj"	"$(INTDIR)\gromos96format.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\hinformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\hinformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\hinformat.obj"	"$(INTDIR)\hinformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\inchiformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\inchiformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\inchiformat.obj"	"$(INTDIR)\inchiformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\jaguarformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\jaguarformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\jaguarformat.obj"	"$(INTDIR)\jaguarformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\kekulize.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\kekulize.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\kekulize.obj"	"$(INTDIR)\kekulize.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\main.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\main.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\main.obj"	"$(INTDIR)\main.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\matrix.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\matrix.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\matrix.obj"	"$(INTDIR)\matrix.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\math\matrix3x3.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\matrix3x3.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\matrix3x3.obj"	"$(INTDIR)\matrix3x3.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\mdlformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\mdlformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\mdlformat.obj"	"$(INTDIR)\mdlformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\mmodformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\mmodformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\mmodformat.obj"	"$(INTDIR)\mmodformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\mol.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\mol.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\mol.obj"	"$(INTDIR)\mol.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\mol2format.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\mol2format.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\mol2format.obj"	"$(INTDIR)\mol2format.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\molchrg.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\molchrg.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\molchrg.obj"	"$(INTDIR)\molchrg.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\mopacformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\mopacformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\mopacformat.obj"	"$(INTDIR)\mopacformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\mpdformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\mpdformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\mpdformat.obj"	"$(INTDIR)\mpdformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\mpqcformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\mpqcformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\mpqcformat.obj"	"$(INTDIR)\mpqcformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\nwchemformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\nwchemformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\nwchemformat.obj"	"$(INTDIR)\nwchemformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\obconversion.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\obconversion.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\obconversion.obj"	"$(INTDIR)\obconversion.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\oberror.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\oberror.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\oberror.obj"	"$(INTDIR)\oberror.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\obiter.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\obiter.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\obiter.obj"	"$(INTDIR)\obiter.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\obutil.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\obutil.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\obutil.obj"	"$(INTDIR)\obutil.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\parsmart.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\parsmart.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\parsmart.obj"	"$(INTDIR)\parsmart.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\patty.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\patty.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\patty.obj"	"$(INTDIR)\patty.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\pcmodelformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\pcmodelformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\pcmodelformat.obj"	"$(INTDIR)\pcmodelformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\pdbformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\pdbformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\pdbformat.obj"	"$(INTDIR)\pdbformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\phmodel.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\phmodel.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\phmodel.obj"	"$(INTDIR)\phmodel.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\povrayformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\povrayformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\povrayformat.obj"	"$(INTDIR)\povrayformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\PQSformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\PQSformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\PQSformat.obj"	"$(INTDIR)\PQSformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\xml\pubchem.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\pubchem.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\pubchem.obj"	"$(INTDIR)\pubchem.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\qchemformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\qchemformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\qchemformat.obj"	"$(INTDIR)\qchemformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\rand.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\rand.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\rand.obj"	"$(INTDIR)\rand.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\reportformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\reportformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\reportformat.obj"	"$(INTDIR)\reportformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\residue.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\residue.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\residue.obj"	"$(INTDIR)\residue.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\ring.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\ring.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\ring.obj"	"$(INTDIR)\ring.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\rotamer.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\rotamer.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\rotamer.obj"	"$(INTDIR)\rotamer.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\rotor.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\rotor.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\rotor.obj"	"$(INTDIR)\rotor.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\rxnformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\rxnformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\rxnformat.obj"	"$(INTDIR)\rxnformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\shelxformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\shelxformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\shelxformat.obj"	"$(INTDIR)\shelxformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\smilesformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\smilesformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\smilesformat.obj"	"$(INTDIR)\smilesformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\tinkerformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\tinkerformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\tinkerformat.obj"	"$(INTDIR)\tinkerformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\tokenst.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\tokenst.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\tokenst.obj"	"$(INTDIR)\tokenst.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\transform.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\transform.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\transform.obj"	"$(INTDIR)\transform.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\turbomoleformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\turbomoleformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\turbomoleformat.obj"	"$(INTDIR)\turbomoleformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\typer.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\typer.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\typer.obj"	"$(INTDIR)\typer.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\unichemformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\unichemformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\unichemformat.obj"	"$(INTDIR)\unichemformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\math\vector3.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\vector3.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\vector3.obj"	"$(INTDIR)\vector3.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\viewmolformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\viewmolformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\viewmolformat.obj"	"$(INTDIR)\viewmolformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\xml\xcmlformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\xcmlformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\xcmlformat.obj"	"$(INTDIR)\xcmlformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\xedformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\xedformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\xedformat.obj"	"$(INTDIR)\xedformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\xml\xml.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\xml.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\xml.obj"	"$(INTDIR)\xml.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\xml\xmlformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\xmlformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\xmlformat.obj"	"$(INTDIR)\xmlformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\xyzformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\xyzformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\xyzformat.obj"	"$(INTDIR)\xyzformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\yasaraformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\yasaraformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\yasaraformat.obj"	"$(INTDIR)\yasaraformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\zindoformat.cpp

!IF  "$(CFG)" == "OBabel - Win32 Release"


"$(INTDIR)\zindoformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"


"$(INTDIR)\zindoformat.obj"	"$(INTDIR)\zindoformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 


!ENDIF 

