# Microsoft Developer Studio Generated NMAKE File, Based on OBGUIs.dsp
!IF "$(CFG)" == ""
CFG=OBGUIs - Win32 Debug
!MESSAGE No configuration specified. Defaulting to OBGUIs - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "OBGUIs - Win32 Release" && "$(CFG)" != "OBGUIs - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "OBGUIs.mak" CFG="OBGUIs - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "OBGUIs - Win32 Release" (based on "Win32 (x86) Application")
!MESSAGE "OBGUIs - Win32 Debug" (based on "Win32 (x86) Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

CPP=cl.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "OBGUIs - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release

ALL : ".\OBGUIs.exe"


CLEAN :
	-@erase "$(INTDIR)\APIInterface.obj"
	-@erase "$(INTDIR)\atom.obj"
	-@erase "$(INTDIR)\base.obj"
	-@erase "$(INTDIR)\bitvec.obj"
	-@erase "$(INTDIR)\bond.obj"
	-@erase "$(INTDIR)\bondtyper.obj"
	-@erase "$(INTDIR)\chains.obj"
	-@erase "$(INTDIR)\chiral.obj"
	-@erase "$(INTDIR)\cmlreactlformat.obj"
	-@erase "$(INTDIR)\copyformat.obj"
	-@erase "$(INTDIR)\data.obj"
	-@erase "$(INTDIR)\dlhandler_win32.obj"
	-@erase "$(INTDIR)\dmolformat.obj"
	-@erase "$(INTDIR)\DynamicOptions.obj"
	-@erase "$(INTDIR)\finger2.obj"
	-@erase "$(INTDIR)\finger3.obj"
	-@erase "$(INTDIR)\fingerprint.obj"
	-@erase "$(INTDIR)\fingerprintformat.obj"
	-@erase "$(INTDIR)\freefracformat.obj"
	-@erase "$(INTDIR)\generic.obj"
	-@erase "$(INTDIR)\grid.obj"
	-@erase "$(INTDIR)\inchiformat.obj"
	-@erase "$(INTDIR)\kekulize.obj"
	-@erase "$(INTDIR)\matrix.obj"
	-@erase "$(INTDIR)\matrix3x3.obj"
	-@erase "$(INTDIR)\mdlformat.obj"
	-@erase "$(INTDIR)\mol.obj"
	-@erase "$(INTDIR)\mol2format.obj"
	-@erase "$(INTDIR)\molchrg.obj"
	-@erase "$(INTDIR)\mpdformat.obj"
	-@erase "$(INTDIR)\obconversion.obj"
	-@erase "$(INTDIR)\oberror.obj"
	-@erase "$(INTDIR)\OBGUI.obj"
	-@erase "$(INTDIR)\OBGUI.res"
	-@erase "$(INTDIR)\OBGUIDlg.obj"
	-@erase "$(INTDIR)\obiter.obj"
	-@erase "$(INTDIR)\obutil.obj"
	-@erase "$(INTDIR)\parsmart.obj"
	-@erase "$(INTDIR)\patty.obj"
	-@erase "$(INTDIR)\phmodel.obj"
	-@erase "$(INTDIR)\pubchem.obj"
	-@erase "$(INTDIR)\qchemformat.obj"
	-@erase "$(INTDIR)\rand.obj"
	-@erase "$(INTDIR)\residue.obj"
	-@erase "$(INTDIR)\ring.obj"
	-@erase "$(INTDIR)\rotamer.obj"
	-@erase "$(INTDIR)\rotor.obj"
	-@erase "$(INTDIR)\rxnformat.obj"
	-@erase "$(INTDIR)\smilesformat.obj"
	-@erase "$(INTDIR)\StdAfx.obj"
	-@erase "$(INTDIR)\tokenst.obj"
	-@erase "$(INTDIR)\transform.obj"
	-@erase "$(INTDIR)\typer.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vector3.obj"
	-@erase "$(INTDIR)\xcmlformat.obj"
	-@erase "$(INTDIR)\xml.obj"
	-@erase "$(INTDIR)\xmlformat.obj"
	-@erase ".\OBGUIs.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP_PROJ=/nologo /MT /W3 /GR /GX /I "..\..\src" /I ".." /I "../../data" /I "..\..\src\formats" /I "..\..\src\formats\xml" /I "..\OBGUI" /D "NDEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "GUI" /D "INCHI_LINK_AS_DLL" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 
MTL_PROJ=/nologo /D "NDEBUG" /mktyplib203 /win32 
RSC_PROJ=/l 0x809 /fo"$(INTDIR)\OBGUI.res" /d "NDEBUG" 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\OBGUIs.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib nafxcw.lib libcmt.lib libinchi.lib Shlwapi.lib /nologo /subsystem:windows /incremental:no /pdb:"$(OUTDIR)\OBGUIs.pdb" /machine:I386 /nodefaultlib:"nafxcw.lib libcmt.lib" /out:"OBGUIs.exe" /libpath:".." /libpath:"..\..\src\formats\cmlpp\\builds\windows\vc6\cmlpplib\release" 
LINK32_OBJS= \
	"$(INTDIR)\APIInterface.obj" \
	"$(INTDIR)\atom.obj" \
	"$(INTDIR)\base.obj" \
	"$(INTDIR)\bitvec.obj" \
	"$(INTDIR)\bond.obj" \
	"$(INTDIR)\bondtyper.obj" \
	"$(INTDIR)\chains.obj" \
	"$(INTDIR)\chiral.obj" \
	"$(INTDIR)\cmlreactlformat.obj" \
	"$(INTDIR)\copyformat.obj" \
	"$(INTDIR)\data.obj" \
	"$(INTDIR)\dlhandler_win32.obj" \
	"$(INTDIR)\dmolformat.obj" \
	"$(INTDIR)\DynamicOptions.obj" \
	"$(INTDIR)\finger2.obj" \
	"$(INTDIR)\finger3.obj" \
	"$(INTDIR)\fingerprint.obj" \
	"$(INTDIR)\fingerprintformat.obj" \
	"$(INTDIR)\freefracformat.obj" \
	"$(INTDIR)\generic.obj" \
	"$(INTDIR)\grid.obj" \
	"$(INTDIR)\inchiformat.obj" \
	"$(INTDIR)\kekulize.obj" \
	"$(INTDIR)\matrix.obj" \
	"$(INTDIR)\matrix3x3.obj" \
	"$(INTDIR)\mdlformat.obj" \
	"$(INTDIR)\mol.obj" \
	"$(INTDIR)\mol2format.obj" \
	"$(INTDIR)\molchrg.obj" \
	"$(INTDIR)\mpdformat.obj" \
	"$(INTDIR)\obconversion.obj" \
	"$(INTDIR)\oberror.obj" \
	"$(INTDIR)\OBGUI.obj" \
	"$(INTDIR)\OBGUIDlg.obj" \
	"$(INTDIR)\obiter.obj" \
	"$(INTDIR)\obutil.obj" \
	"$(INTDIR)\parsmart.obj" \
	"$(INTDIR)\patty.obj" \
	"$(INTDIR)\phmodel.obj" \
	"$(INTDIR)\pubchem.obj" \
	"$(INTDIR)\qchemformat.obj" \
	"$(INTDIR)\rand.obj" \
	"$(INTDIR)\residue.obj" \
	"$(INTDIR)\ring.obj" \
	"$(INTDIR)\rotamer.obj" \
	"$(INTDIR)\rotor.obj" \
	"$(INTDIR)\rxnformat.obj" \
	"$(INTDIR)\smilesformat.obj" \
	"$(INTDIR)\StdAfx.obj" \
	"$(INTDIR)\tokenst.obj" \
	"$(INTDIR)\transform.obj" \
	"$(INTDIR)\typer.obj" \
	"$(INTDIR)\vector3.obj" \
	"$(INTDIR)\xcmlformat.obj" \
	"$(INTDIR)\xml.obj" \
	"$(INTDIR)\xmlformat.obj" \
	"$(INTDIR)\OBGUI.res" \
	"..\libxml2.lib" \
	"..\zdll.lib"

".\OBGUIs.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"

OUTDIR=.\Debug
INTDIR=.\Debug
# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

ALL : "$(OUTDIR)\OBGUIs.exe" "$(OUTDIR)\OBGUIs.bsc"


CLEAN :
	-@erase "$(INTDIR)\APIInterface.obj"
	-@erase "$(INTDIR)\APIInterface.sbr"
	-@erase "$(INTDIR)\atom.obj"
	-@erase "$(INTDIR)\atom.sbr"
	-@erase "$(INTDIR)\base.obj"
	-@erase "$(INTDIR)\base.sbr"
	-@erase "$(INTDIR)\bitvec.obj"
	-@erase "$(INTDIR)\bitvec.sbr"
	-@erase "$(INTDIR)\bond.obj"
	-@erase "$(INTDIR)\bond.sbr"
	-@erase "$(INTDIR)\bondtyper.obj"
	-@erase "$(INTDIR)\bondtyper.sbr"
	-@erase "$(INTDIR)\chains.obj"
	-@erase "$(INTDIR)\chains.sbr"
	-@erase "$(INTDIR)\chiral.obj"
	-@erase "$(INTDIR)\chiral.sbr"
	-@erase "$(INTDIR)\cmlreactlformat.obj"
	-@erase "$(INTDIR)\cmlreactlformat.sbr"
	-@erase "$(INTDIR)\copyformat.obj"
	-@erase "$(INTDIR)\copyformat.sbr"
	-@erase "$(INTDIR)\data.obj"
	-@erase "$(INTDIR)\data.sbr"
	-@erase "$(INTDIR)\dlhandler_win32.obj"
	-@erase "$(INTDIR)\dlhandler_win32.sbr"
	-@erase "$(INTDIR)\dmolformat.obj"
	-@erase "$(INTDIR)\dmolformat.sbr"
	-@erase "$(INTDIR)\DynamicOptions.obj"
	-@erase "$(INTDIR)\DynamicOptions.sbr"
	-@erase "$(INTDIR)\fastsearchformat.obj"
	-@erase "$(INTDIR)\fastsearchformat.sbr"
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
	-@erase "$(INTDIR)\generic.obj"
	-@erase "$(INTDIR)\generic.sbr"
	-@erase "$(INTDIR)\grid.obj"
	-@erase "$(INTDIR)\grid.sbr"
	-@erase "$(INTDIR)\inchiformat.obj"
	-@erase "$(INTDIR)\inchiformat.sbr"
	-@erase "$(INTDIR)\kekulize.obj"
	-@erase "$(INTDIR)\kekulize.sbr"
	-@erase "$(INTDIR)\matrix.obj"
	-@erase "$(INTDIR)\matrix.sbr"
	-@erase "$(INTDIR)\matrix3x3.obj"
	-@erase "$(INTDIR)\matrix3x3.sbr"
	-@erase "$(INTDIR)\mdlformat.obj"
	-@erase "$(INTDIR)\mdlformat.sbr"
	-@erase "$(INTDIR)\mol.obj"
	-@erase "$(INTDIR)\mol.sbr"
	-@erase "$(INTDIR)\mol2format.obj"
	-@erase "$(INTDIR)\mol2format.sbr"
	-@erase "$(INTDIR)\molchrg.obj"
	-@erase "$(INTDIR)\molchrg.sbr"
	-@erase "$(INTDIR)\mpdformat.obj"
	-@erase "$(INTDIR)\mpdformat.sbr"
	-@erase "$(INTDIR)\obconversion.obj"
	-@erase "$(INTDIR)\obconversion.sbr"
	-@erase "$(INTDIR)\oberror.obj"
	-@erase "$(INTDIR)\oberror.sbr"
	-@erase "$(INTDIR)\OBGUI.obj"
	-@erase "$(INTDIR)\OBGUI.res"
	-@erase "$(INTDIR)\OBGUI.sbr"
	-@erase "$(INTDIR)\OBGUIDlg.obj"
	-@erase "$(INTDIR)\OBGUIDlg.sbr"
	-@erase "$(INTDIR)\obiter.obj"
	-@erase "$(INTDIR)\obiter.sbr"
	-@erase "$(INTDIR)\obutil.obj"
	-@erase "$(INTDIR)\obutil.sbr"
	-@erase "$(INTDIR)\parsmart.obj"
	-@erase "$(INTDIR)\parsmart.sbr"
	-@erase "$(INTDIR)\patty.obj"
	-@erase "$(INTDIR)\patty.sbr"
	-@erase "$(INTDIR)\phmodel.obj"
	-@erase "$(INTDIR)\phmodel.sbr"
	-@erase "$(INTDIR)\pubchem.obj"
	-@erase "$(INTDIR)\pubchem.sbr"
	-@erase "$(INTDIR)\qchemformat.obj"
	-@erase "$(INTDIR)\qchemformat.sbr"
	-@erase "$(INTDIR)\rand.obj"
	-@erase "$(INTDIR)\rand.sbr"
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
	-@erase "$(INTDIR)\smilesformat.obj"
	-@erase "$(INTDIR)\smilesformat.sbr"
	-@erase "$(INTDIR)\StdAfx.obj"
	-@erase "$(INTDIR)\StdAfx.sbr"
	-@erase "$(INTDIR)\tokenst.obj"
	-@erase "$(INTDIR)\tokenst.sbr"
	-@erase "$(INTDIR)\transform.obj"
	-@erase "$(INTDIR)\transform.sbr"
	-@erase "$(INTDIR)\typer.obj"
	-@erase "$(INTDIR)\typer.sbr"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vc60.pdb"
	-@erase "$(INTDIR)\vector3.obj"
	-@erase "$(INTDIR)\vector3.sbr"
	-@erase "$(INTDIR)\xcmlformat.obj"
	-@erase "$(INTDIR)\xcmlformat.sbr"
	-@erase "$(INTDIR)\xml.obj"
	-@erase "$(INTDIR)\xml.sbr"
	-@erase "$(INTDIR)\xmlformat.obj"
	-@erase "$(INTDIR)\xmlformat.sbr"
	-@erase "$(OUTDIR)\OBGUIs.bsc"
	-@erase "$(OUTDIR)\OBGUIs.exe"
	-@erase "$(OUTDIR)\OBGUIs.ilk"
	-@erase "$(OUTDIR)\OBGUIs.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP_PROJ=/nologo /MTd /W3 /Gm /GR /GX /ZI /Od /I "..\..\src\formats\xml" /I "..\..\src" /I ".." /I "../../data" /I "..\..\src\formats" /I "..\OBGUI" /D "_DEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "GUI" /D "INCHI_LINK_AS_DLL" /FR"$(INTDIR)\\" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 
MTL_PROJ=/nologo /D "_DEBUG" /mktyplib203 /win32 
RSC_PROJ=/l 0x809 /fo"$(INTDIR)\OBGUI.res" /d "_DEBUG" 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\OBGUIs.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\APIInterface.sbr" \
	"$(INTDIR)\atom.sbr" \
	"$(INTDIR)\base.sbr" \
	"$(INTDIR)\bitvec.sbr" \
	"$(INTDIR)\bond.sbr" \
	"$(INTDIR)\bondtyper.sbr" \
	"$(INTDIR)\chains.sbr" \
	"$(INTDIR)\chiral.sbr" \
	"$(INTDIR)\cmlreactlformat.sbr" \
	"$(INTDIR)\copyformat.sbr" \
	"$(INTDIR)\data.sbr" \
	"$(INTDIR)\dlhandler_win32.sbr" \
	"$(INTDIR)\dmolformat.sbr" \
	"$(INTDIR)\DynamicOptions.sbr" \
	"$(INTDIR)\fastsearchformat.sbr" \
	"$(INTDIR)\finger2.sbr" \
	"$(INTDIR)\finger3.sbr" \
	"$(INTDIR)\fingerprint.sbr" \
	"$(INTDIR)\fingerprintformat.sbr" \
	"$(INTDIR)\freefracformat.sbr" \
	"$(INTDIR)\generic.sbr" \
	"$(INTDIR)\grid.sbr" \
	"$(INTDIR)\inchiformat.sbr" \
	"$(INTDIR)\kekulize.sbr" \
	"$(INTDIR)\matrix.sbr" \
	"$(INTDIR)\matrix3x3.sbr" \
	"$(INTDIR)\mdlformat.sbr" \
	"$(INTDIR)\mol.sbr" \
	"$(INTDIR)\mol2format.sbr" \
	"$(INTDIR)\molchrg.sbr" \
	"$(INTDIR)\mpdformat.sbr" \
	"$(INTDIR)\obconversion.sbr" \
	"$(INTDIR)\oberror.sbr" \
	"$(INTDIR)\OBGUI.sbr" \
	"$(INTDIR)\OBGUIDlg.sbr" \
	"$(INTDIR)\obiter.sbr" \
	"$(INTDIR)\obutil.sbr" \
	"$(INTDIR)\parsmart.sbr" \
	"$(INTDIR)\patty.sbr" \
	"$(INTDIR)\phmodel.sbr" \
	"$(INTDIR)\pubchem.sbr" \
	"$(INTDIR)\qchemformat.sbr" \
	"$(INTDIR)\rand.sbr" \
	"$(INTDIR)\residue.sbr" \
	"$(INTDIR)\ring.sbr" \
	"$(INTDIR)\rotamer.sbr" \
	"$(INTDIR)\rotor.sbr" \
	"$(INTDIR)\rxnformat.sbr" \
	"$(INTDIR)\smilesformat.sbr" \
	"$(INTDIR)\StdAfx.sbr" \
	"$(INTDIR)\tokenst.sbr" \
	"$(INTDIR)\transform.sbr" \
	"$(INTDIR)\typer.sbr" \
	"$(INTDIR)\vector3.sbr" \
	"$(INTDIR)\xcmlformat.sbr" \
	"$(INTDIR)\xml.sbr" \
	"$(INTDIR)\xmlformat.sbr"

"$(OUTDIR)\OBGUIs.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib nafxcwd.lib libcmtd.lib libinchi.lib Shlwapi.lib /nologo /subsystem:windows /incremental:yes /pdb:"$(OUTDIR)\OBGUIs.pdb" /debug /machine:I386 /nodefaultlib:"nafxcwd.lib libcmtd.lib" /out:"$(OUTDIR)\OBGUIs.exe" /libpath:".." /libpath:"..\..\src\formats\cmlpp\\builds\windows\vc6\cmlpplib\debug" 
LINK32_OBJS= \
	"$(INTDIR)\APIInterface.obj" \
	"$(INTDIR)\atom.obj" \
	"$(INTDIR)\base.obj" \
	"$(INTDIR)\bitvec.obj" \
	"$(INTDIR)\bond.obj" \
	"$(INTDIR)\bondtyper.obj" \
	"$(INTDIR)\chains.obj" \
	"$(INTDIR)\chiral.obj" \
	"$(INTDIR)\cmlreactlformat.obj" \
	"$(INTDIR)\copyformat.obj" \
	"$(INTDIR)\data.obj" \
	"$(INTDIR)\dlhandler_win32.obj" \
	"$(INTDIR)\dmolformat.obj" \
	"$(INTDIR)\DynamicOptions.obj" \
	"$(INTDIR)\fastsearchformat.obj" \
	"$(INTDIR)\finger2.obj" \
	"$(INTDIR)\finger3.obj" \
	"$(INTDIR)\fingerprint.obj" \
	"$(INTDIR)\fingerprintformat.obj" \
	"$(INTDIR)\freefracformat.obj" \
	"$(INTDIR)\generic.obj" \
	"$(INTDIR)\grid.obj" \
	"$(INTDIR)\inchiformat.obj" \
	"$(INTDIR)\kekulize.obj" \
	"$(INTDIR)\matrix.obj" \
	"$(INTDIR)\matrix3x3.obj" \
	"$(INTDIR)\mdlformat.obj" \
	"$(INTDIR)\mol.obj" \
	"$(INTDIR)\mol2format.obj" \
	"$(INTDIR)\molchrg.obj" \
	"$(INTDIR)\mpdformat.obj" \
	"$(INTDIR)\obconversion.obj" \
	"$(INTDIR)\oberror.obj" \
	"$(INTDIR)\OBGUI.obj" \
	"$(INTDIR)\OBGUIDlg.obj" \
	"$(INTDIR)\obiter.obj" \
	"$(INTDIR)\obutil.obj" \
	"$(INTDIR)\parsmart.obj" \
	"$(INTDIR)\patty.obj" \
	"$(INTDIR)\phmodel.obj" \
	"$(INTDIR)\pubchem.obj" \
	"$(INTDIR)\qchemformat.obj" \
	"$(INTDIR)\rand.obj" \
	"$(INTDIR)\residue.obj" \
	"$(INTDIR)\ring.obj" \
	"$(INTDIR)\rotamer.obj" \
	"$(INTDIR)\rotor.obj" \
	"$(INTDIR)\rxnformat.obj" \
	"$(INTDIR)\smilesformat.obj" \
	"$(INTDIR)\StdAfx.obj" \
	"$(INTDIR)\tokenst.obj" \
	"$(INTDIR)\transform.obj" \
	"$(INTDIR)\typer.obj" \
	"$(INTDIR)\vector3.obj" \
	"$(INTDIR)\xcmlformat.obj" \
	"$(INTDIR)\xml.obj" \
	"$(INTDIR)\xmlformat.obj" \
	"$(INTDIR)\OBGUI.res" \
	"..\libxml2.lib" \
	"..\zdll.lib"

"$(OUTDIR)\OBGUIs.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
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
!IF EXISTS("OBGUIs.dep")
!INCLUDE "OBGUIs.dep"
!ELSE 
!MESSAGE Warning: cannot find "OBGUIs.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "OBGUIs - Win32 Release" || "$(CFG)" == "OBGUIs - Win32 Debug"
SOURCE=..\..\src\formats\APIInterface.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\APIInterface.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\APIInterface.obj"	"$(INTDIR)\APIInterface.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\assignbonds.cpp
SOURCE=..\..\src\atom.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\atom.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\atom.obj"	"$(INTDIR)\atom.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\base.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\base.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\base.obj"	"$(INTDIR)\base.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\bitvec.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\bitvec.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\bitvec.obj"	"$(INTDIR)\bitvec.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\bond.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\bond.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\bond.obj"	"$(INTDIR)\bond.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\bondtyper.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\bondtyper.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\bondtyper.obj"	"$(INTDIR)\bondtyper.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\chains.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\chains.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\chains.obj"	"$(INTDIR)\chains.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\chiral.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\chiral.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\chiral.obj"	"$(INTDIR)\chiral.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\xml\cmlreactlformat.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\cmlreactlformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\cmlreactlformat.obj"	"$(INTDIR)\cmlreactlformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\copyformat.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\copyformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\copyformat.obj"	"$(INTDIR)\copyformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\data.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\data.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\data.obj"	"$(INTDIR)\data.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\dlhandler_win32.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\dlhandler_win32.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\dlhandler_win32.obj"	"$(INTDIR)\dlhandler_win32.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\dmolformat.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\dmolformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\dmolformat.obj"	"$(INTDIR)\dmolformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\OBGUI\DynamicOptions.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\DynamicOptions.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\DynamicOptions.obj"	"$(INTDIR)\DynamicOptions.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\fastsearchformat.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"

!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\fastsearchformat.obj"	"$(INTDIR)\fastsearchformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\fingerprints\finger2.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\finger2.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\finger2.obj"	"$(INTDIR)\finger2.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\fingerprints\finger3.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\finger3.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\finger3.obj"	"$(INTDIR)\finger3.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\fingerprint.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\fingerprint.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\fingerprint.obj"	"$(INTDIR)\fingerprint.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\fingerprintformat.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\fingerprintformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\fingerprintformat.obj"	"$(INTDIR)\fingerprintformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\freefracformat.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\freefracformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\freefracformat.obj"	"$(INTDIR)\freefracformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\generic.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\generic.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\generic.obj"	"$(INTDIR)\generic.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\grid.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\grid.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\grid.obj"	"$(INTDIR)\grid.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\inchiformat.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\inchiformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\inchiformat.obj"	"$(INTDIR)\inchiformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\kekulize.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\kekulize.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\kekulize.obj"	"$(INTDIR)\kekulize.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\matrix.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\matrix.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\matrix.obj"	"$(INTDIR)\matrix.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\math\matrix3x3.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\matrix3x3.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\matrix3x3.obj"	"$(INTDIR)\matrix3x3.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\mdlformat.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\mdlformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\mdlformat.obj"	"$(INTDIR)\mdlformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\mol.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\mol.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\mol.obj"	"$(INTDIR)\mol.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\mol2format.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\mol2format.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\mol2format.obj"	"$(INTDIR)\mol2format.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\molchrg.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\molchrg.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\molchrg.obj"	"$(INTDIR)\molchrg.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\mpdformat.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\mpdformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\mpdformat.obj"	"$(INTDIR)\mpdformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\obconversion.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\obconversion.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\obconversion.obj"	"$(INTDIR)\obconversion.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\oberror.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\oberror.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\oberror.obj"	"$(INTDIR)\oberror.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\OBGUI\OBGUI.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\OBGUI.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\OBGUI.obj"	"$(INTDIR)\OBGUI.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\OBGUI\OBGUI.rc

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\OBGUI.res" : $(SOURCE) "$(INTDIR)"
	$(RSC) /l 0x809 /fo"$(INTDIR)\OBGUI.res" /i "\My Documents\MSVC\OpenBabel Ultimate\CVSforOB2\openbabel\windows\OBGUI" /d "NDEBUG" $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\OBGUI.res" : $(SOURCE) "$(INTDIR)"
	$(RSC) /l 0x809 /fo"$(INTDIR)\OBGUI.res" /i "\My Documents\MSVC\OpenBabel Ultimate\CVSforOB2\openbabel\windows\OBGUI" /d "_DEBUG" $(SOURCE)


!ENDIF 

SOURCE=..\OBGUI\OBGUIDlg.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\OBGUIDlg.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\OBGUIDlg.obj"	"$(INTDIR)\OBGUIDlg.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\obiter.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\obiter.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\obiter.obj"	"$(INTDIR)\obiter.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\obutil.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\obutil.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\obutil.obj"	"$(INTDIR)\obutil.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\parsmart.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\parsmart.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\parsmart.obj"	"$(INTDIR)\parsmart.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\patty.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\patty.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\patty.obj"	"$(INTDIR)\patty.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\phmodel.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\phmodel.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\phmodel.obj"	"$(INTDIR)\phmodel.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\xml\pubchem.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\pubchem.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\pubchem.obj"	"$(INTDIR)\pubchem.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\qchemformat.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\qchemformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\qchemformat.obj"	"$(INTDIR)\qchemformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\rand.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\rand.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\rand.obj"	"$(INTDIR)\rand.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\residue.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\residue.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\residue.obj"	"$(INTDIR)\residue.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\ring.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\ring.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\ring.obj"	"$(INTDIR)\ring.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\rotamer.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\rotamer.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\rotamer.obj"	"$(INTDIR)\rotamer.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\rotor.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\rotor.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\rotor.obj"	"$(INTDIR)\rotor.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\rxnformat.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\rxnformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\rxnformat.obj"	"$(INTDIR)\rxnformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\smilesformat.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\smilesformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\smilesformat.obj"	"$(INTDIR)\smilesformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\OBGUI\StdAfx.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\StdAfx.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\StdAfx.obj"	"$(INTDIR)\StdAfx.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\tokenst.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\tokenst.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\tokenst.obj"	"$(INTDIR)\tokenst.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\transform.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\transform.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\transform.obj"	"$(INTDIR)\transform.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\typer.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\typer.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\typer.obj"	"$(INTDIR)\typer.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\math\vector3.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\vector3.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\vector3.obj"	"$(INTDIR)\vector3.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\xml\xcmlformat.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\xcmlformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\xcmlformat.obj"	"$(INTDIR)\xcmlformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\xml\xml.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\xml.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\xml.obj"	"$(INTDIR)\xml.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\xml\xmlformat.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\xmlformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\xmlformat.obj"	"$(INTDIR)\xmlformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 


!ENDIF 

