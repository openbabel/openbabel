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

!IF  "$(CFG)" == "OBGUIs - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release

ALL : ".\OBGUIs.exe"


CLEAN :
	-@erase "$(INTDIR)\atom.obj"
	-@erase "$(INTDIR)\base.obj"
	-@erase "$(INTDIR)\bgf.obj"
	-@erase "$(INTDIR)\binary.obj"
	-@erase "$(INTDIR)\bitvec.obj"
	-@erase "$(INTDIR)\bond.obj"
	-@erase "$(INTDIR)\box.obj"
	-@erase "$(INTDIR)\c3d.obj"
	-@erase "$(INTDIR)\cache.obj"
	-@erase "$(INTDIR)\car.obj"
	-@erase "$(INTDIR)\ccc.obj"
	-@erase "$(INTDIR)\chains.obj"
	-@erase "$(INTDIR)\chdrw.obj"
	-@erase "$(INTDIR)\chiral.obj"
	-@erase "$(INTDIR)\cml.obj"
	-@erase "$(INTDIR)\cmlformat.obj"
	-@erase "$(INTDIR)\csr.obj"
	-@erase "$(INTDIR)\cssr.obj"
	-@erase "$(INTDIR)\data.obj"
	-@erase "$(INTDIR)\dlhandler_win32..obj"
	-@erase "$(INTDIR)\dmol.obj"
	-@erase "$(INTDIR)\DynamicOptions.obj"
	-@erase "$(INTDIR)\generic.obj"
	-@erase "$(INTDIR)\grid.obj"
	-@erase "$(INTDIR)\hin.obj"
	-@erase "$(INTDIR)\matrix.obj"
	-@erase "$(INTDIR)\matrix3x3.obj"
	-@erase "$(INTDIR)\mdlformat.obj"
	-@erase "$(INTDIR)\mol.obj"
	-@erase "$(INTDIR)\mol2format.obj"
	-@erase "$(INTDIR)\molchrg.obj"
	-@erase "$(INTDIR)\obconversion.obj"
	-@erase "$(INTDIR)\oberror.obj"
	-@erase "$(INTDIR)\OBGUI.obj"
	-@erase "$(INTDIR)\OBGUI.res"
	-@erase "$(INTDIR)\OBGUIDlg.obj"
	-@erase "$(INTDIR)\obutil.obj"
	-@erase "$(INTDIR)\parsmart.obj"
	-@erase "$(INTDIR)\patty.obj"
	-@erase "$(INTDIR)\phmodel.obj"
	-@erase "$(INTDIR)\rand.obj"
	-@erase "$(INTDIR)\report.obj"
	-@erase "$(INTDIR)\residue.obj"
	-@erase "$(INTDIR)\ring.obj"
	-@erase "$(INTDIR)\rotor.obj"
	-@erase "$(INTDIR)\rxnformat.obj"
	-@erase "$(INTDIR)\smilesformat.obj"
	-@erase "$(INTDIR)\StdAfx.obj"
	-@erase "$(INTDIR)\tinker.obj"
	-@erase "$(INTDIR)\tokenst.obj"
	-@erase "$(INTDIR)\transform.obj"
	-@erase "$(INTDIR)\turbomoleformat.obj"
	-@erase "$(INTDIR)\typer.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vector3.obj"
	-@erase ".\OBGUIs.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MT /W3 /GR /GX /I "..\obgui" /I "..\..\src" /D "NDEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "GUI" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

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

MTL=midl.exe
MTL_PROJ=/nologo /D "NDEBUG" /mktyplib203 /win32 
RSC=rc.exe
RSC_PROJ=/l 0x809 /fo"$(INTDIR)\OBGUI.res" /d "NDEBUG" 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\OBGUIs.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=nafxcw.lib libcmt.lib /nologo /subsystem:windows /incremental:no /pdb:"$(OUTDIR)\OBGUIs.pdb" /machine:I386 /nodefaultlib:"nafxcw.lib libcmt.lib" /out:"OBGUIs.exe" 
LINK32_OBJS= \
	"$(INTDIR)\atom.obj" \
	"$(INTDIR)\base.obj" \
	"$(INTDIR)\bgf.obj" \
	"$(INTDIR)\binary.obj" \
	"$(INTDIR)\bitvec.obj" \
	"$(INTDIR)\bond.obj" \
	"$(INTDIR)\box.obj" \
	"$(INTDIR)\c3d.obj" \
	"$(INTDIR)\cache.obj" \
	"$(INTDIR)\car.obj" \
	"$(INTDIR)\ccc.obj" \
	"$(INTDIR)\chains.obj" \
	"$(INTDIR)\chdrw.obj" \
	"$(INTDIR)\chiral.obj" \
	"$(INTDIR)\cml.obj" \
	"$(INTDIR)\cmlformat.obj" \
	"$(INTDIR)\csr.obj" \
	"$(INTDIR)\cssr.obj" \
	"$(INTDIR)\data.obj" \
	"$(INTDIR)\dlhandler_win32..obj" \
	"$(INTDIR)\dmol.obj" \
	"$(INTDIR)\DynamicOptions.obj" \
	"$(INTDIR)\generic.obj" \
	"$(INTDIR)\grid.obj" \
	"$(INTDIR)\hin.obj" \
	"$(INTDIR)\matrix.obj" \
	"$(INTDIR)\matrix3x3.obj" \
	"$(INTDIR)\mdlformat.obj" \
	"$(INTDIR)\mol.obj" \
	"$(INTDIR)\mol2format.obj" \
	"$(INTDIR)\molchrg.obj" \
	"$(INTDIR)\obconversion.obj" \
	"$(INTDIR)\oberror.obj" \
	"$(INTDIR)\OBGUI.obj" \
	"$(INTDIR)\OBGUIDlg.obj" \
	"$(INTDIR)\obutil.obj" \
	"$(INTDIR)\parsmart.obj" \
	"$(INTDIR)\patty.obj" \
	"$(INTDIR)\phmodel.obj" \
	"$(INTDIR)\rand.obj" \
	"$(INTDIR)\report.obj" \
	"$(INTDIR)\residue.obj" \
	"$(INTDIR)\ring.obj" \
	"$(INTDIR)\rotor.obj" \
	"$(INTDIR)\rxnformat.obj" \
	"$(INTDIR)\smilesformat.obj" \
	"$(INTDIR)\StdAfx.obj" \
	"$(INTDIR)\tinker.obj" \
	"$(INTDIR)\tokenst.obj" \
	"$(INTDIR)\transform.obj" \
	"$(INTDIR)\turbomoleformat.obj" \
	"$(INTDIR)\typer.obj" \
	"$(INTDIR)\vector3.obj" \
	"$(INTDIR)\OBGUI.res"

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
	-@erase "$(INTDIR)\atom.obj"
	-@erase "$(INTDIR)\atom.sbr"
	-@erase "$(INTDIR)\base.obj"
	-@erase "$(INTDIR)\base.sbr"
	-@erase "$(INTDIR)\bgf.obj"
	-@erase "$(INTDIR)\bgf.sbr"
	-@erase "$(INTDIR)\binary.obj"
	-@erase "$(INTDIR)\binary.sbr"
	-@erase "$(INTDIR)\bitvec.obj"
	-@erase "$(INTDIR)\bitvec.sbr"
	-@erase "$(INTDIR)\bond.obj"
	-@erase "$(INTDIR)\bond.sbr"
	-@erase "$(INTDIR)\box.obj"
	-@erase "$(INTDIR)\box.sbr"
	-@erase "$(INTDIR)\c3d.obj"
	-@erase "$(INTDIR)\c3d.sbr"
	-@erase "$(INTDIR)\cache.obj"
	-@erase "$(INTDIR)\cache.sbr"
	-@erase "$(INTDIR)\car.obj"
	-@erase "$(INTDIR)\car.sbr"
	-@erase "$(INTDIR)\ccc.obj"
	-@erase "$(INTDIR)\ccc.sbr"
	-@erase "$(INTDIR)\chains.obj"
	-@erase "$(INTDIR)\chains.sbr"
	-@erase "$(INTDIR)\chdrw.obj"
	-@erase "$(INTDIR)\chdrw.sbr"
	-@erase "$(INTDIR)\chiral.obj"
	-@erase "$(INTDIR)\chiral.sbr"
	-@erase "$(INTDIR)\cml.obj"
	-@erase "$(INTDIR)\cml.sbr"
	-@erase "$(INTDIR)\cmlformat.obj"
	-@erase "$(INTDIR)\cmlformat.sbr"
	-@erase "$(INTDIR)\csr.obj"
	-@erase "$(INTDIR)\csr.sbr"
	-@erase "$(INTDIR)\cssr.obj"
	-@erase "$(INTDIR)\cssr.sbr"
	-@erase "$(INTDIR)\data.obj"
	-@erase "$(INTDIR)\data.sbr"
	-@erase "$(INTDIR)\dlhandler_win32..obj"
	-@erase "$(INTDIR)\dlhandler_win32..sbr"
	-@erase "$(INTDIR)\dmol.obj"
	-@erase "$(INTDIR)\dmol.sbr"
	-@erase "$(INTDIR)\DynamicOptions.obj"
	-@erase "$(INTDIR)\DynamicOptions.sbr"
	-@erase "$(INTDIR)\generic.obj"
	-@erase "$(INTDIR)\generic.sbr"
	-@erase "$(INTDIR)\grid.obj"
	-@erase "$(INTDIR)\grid.sbr"
	-@erase "$(INTDIR)\hin.obj"
	-@erase "$(INTDIR)\hin.sbr"
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
	-@erase "$(INTDIR)\obconversion.obj"
	-@erase "$(INTDIR)\obconversion.sbr"
	-@erase "$(INTDIR)\oberror.obj"
	-@erase "$(INTDIR)\oberror.sbr"
	-@erase "$(INTDIR)\OBGUI.obj"
	-@erase "$(INTDIR)\OBGUI.res"
	-@erase "$(INTDIR)\OBGUI.sbr"
	-@erase "$(INTDIR)\OBGUIDlg.obj"
	-@erase "$(INTDIR)\OBGUIDlg.sbr"
	-@erase "$(INTDIR)\obutil.obj"
	-@erase "$(INTDIR)\obutil.sbr"
	-@erase "$(INTDIR)\parsmart.obj"
	-@erase "$(INTDIR)\parsmart.sbr"
	-@erase "$(INTDIR)\patty.obj"
	-@erase "$(INTDIR)\patty.sbr"
	-@erase "$(INTDIR)\phmodel.obj"
	-@erase "$(INTDIR)\phmodel.sbr"
	-@erase "$(INTDIR)\rand.obj"
	-@erase "$(INTDIR)\rand.sbr"
	-@erase "$(INTDIR)\report.obj"
	-@erase "$(INTDIR)\report.sbr"
	-@erase "$(INTDIR)\residue.obj"
	-@erase "$(INTDIR)\residue.sbr"
	-@erase "$(INTDIR)\ring.obj"
	-@erase "$(INTDIR)\ring.sbr"
	-@erase "$(INTDIR)\rotor.obj"
	-@erase "$(INTDIR)\rotor.sbr"
	-@erase "$(INTDIR)\rxnformat.obj"
	-@erase "$(INTDIR)\rxnformat.sbr"
	-@erase "$(INTDIR)\smilesformat.obj"
	-@erase "$(INTDIR)\smilesformat.sbr"
	-@erase "$(INTDIR)\StdAfx.obj"
	-@erase "$(INTDIR)\StdAfx.sbr"
	-@erase "$(INTDIR)\tinker.obj"
	-@erase "$(INTDIR)\tinker.sbr"
	-@erase "$(INTDIR)\tokenst.obj"
	-@erase "$(INTDIR)\tokenst.sbr"
	-@erase "$(INTDIR)\transform.obj"
	-@erase "$(INTDIR)\transform.sbr"
	-@erase "$(INTDIR)\turbomoleformat.obj"
	-@erase "$(INTDIR)\turbomoleformat.sbr"
	-@erase "$(INTDIR)\typer.obj"
	-@erase "$(INTDIR)\typer.sbr"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vc60.pdb"
	-@erase "$(INTDIR)\vector3.obj"
	-@erase "$(INTDIR)\vector3.sbr"
	-@erase "$(OUTDIR)\OBGUIs.bsc"
	-@erase "$(OUTDIR)\OBGUIs.exe"
	-@erase "$(OUTDIR)\OBGUIs.ilk"
	-@erase "$(OUTDIR)\OBGUIs.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MTd /W3 /Gm /GR /GX /ZI /Od /I "..\obgui" /I "..\..\src" /D "_DEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "GUI" /FR"$(INTDIR)\\" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

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

MTL=midl.exe
MTL_PROJ=/nologo /D "_DEBUG" /mktyplib203 /win32 
RSC=rc.exe
RSC_PROJ=/l 0x809 /fo"$(INTDIR)\OBGUI.res" /d "_DEBUG" 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\OBGUIs.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\atom.sbr" \
	"$(INTDIR)\base.sbr" \
	"$(INTDIR)\bgf.sbr" \
	"$(INTDIR)\binary.sbr" \
	"$(INTDIR)\bitvec.sbr" \
	"$(INTDIR)\bond.sbr" \
	"$(INTDIR)\box.sbr" \
	"$(INTDIR)\c3d.sbr" \
	"$(INTDIR)\cache.sbr" \
	"$(INTDIR)\car.sbr" \
	"$(INTDIR)\ccc.sbr" \
	"$(INTDIR)\chains.sbr" \
	"$(INTDIR)\chdrw.sbr" \
	"$(INTDIR)\chiral.sbr" \
	"$(INTDIR)\cml.sbr" \
	"$(INTDIR)\cmlformat.sbr" \
	"$(INTDIR)\csr.sbr" \
	"$(INTDIR)\cssr.sbr" \
	"$(INTDIR)\data.sbr" \
	"$(INTDIR)\dlhandler_win32..sbr" \
	"$(INTDIR)\dmol.sbr" \
	"$(INTDIR)\DynamicOptions.sbr" \
	"$(INTDIR)\generic.sbr" \
	"$(INTDIR)\grid.sbr" \
	"$(INTDIR)\hin.sbr" \
	"$(INTDIR)\matrix.sbr" \
	"$(INTDIR)\matrix3x3.sbr" \
	"$(INTDIR)\mdlformat.sbr" \
	"$(INTDIR)\mol.sbr" \
	"$(INTDIR)\mol2format.sbr" \
	"$(INTDIR)\molchrg.sbr" \
	"$(INTDIR)\obconversion.sbr" \
	"$(INTDIR)\oberror.sbr" \
	"$(INTDIR)\OBGUI.sbr" \
	"$(INTDIR)\OBGUIDlg.sbr" \
	"$(INTDIR)\obutil.sbr" \
	"$(INTDIR)\parsmart.sbr" \
	"$(INTDIR)\patty.sbr" \
	"$(INTDIR)\phmodel.sbr" \
	"$(INTDIR)\rand.sbr" \
	"$(INTDIR)\report.sbr" \
	"$(INTDIR)\residue.sbr" \
	"$(INTDIR)\ring.sbr" \
	"$(INTDIR)\rotor.sbr" \
	"$(INTDIR)\rxnformat.sbr" \
	"$(INTDIR)\smilesformat.sbr" \
	"$(INTDIR)\StdAfx.sbr" \
	"$(INTDIR)\tinker.sbr" \
	"$(INTDIR)\tokenst.sbr" \
	"$(INTDIR)\transform.sbr" \
	"$(INTDIR)\turbomoleformat.sbr" \
	"$(INTDIR)\typer.sbr" \
	"$(INTDIR)\vector3.sbr"

"$(OUTDIR)\OBGUIs.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib nafxcwd.lib libcmtd.lib /nologo /subsystem:windows /incremental:yes /pdb:"$(OUTDIR)\OBGUIs.pdb" /debug /machine:I386 /nodefaultlib:"nafxcwd.lib libcmtd.lib" /out:"$(OUTDIR)\OBGUIs.exe" /pdbtype:sept 
LINK32_OBJS= \
	"$(INTDIR)\atom.obj" \
	"$(INTDIR)\base.obj" \
	"$(INTDIR)\bgf.obj" \
	"$(INTDIR)\binary.obj" \
	"$(INTDIR)\bitvec.obj" \
	"$(INTDIR)\bond.obj" \
	"$(INTDIR)\box.obj" \
	"$(INTDIR)\c3d.obj" \
	"$(INTDIR)\cache.obj" \
	"$(INTDIR)\car.obj" \
	"$(INTDIR)\ccc.obj" \
	"$(INTDIR)\chains.obj" \
	"$(INTDIR)\chdrw.obj" \
	"$(INTDIR)\chiral.obj" \
	"$(INTDIR)\cml.obj" \
	"$(INTDIR)\cmlformat.obj" \
	"$(INTDIR)\csr.obj" \
	"$(INTDIR)\cssr.obj" \
	"$(INTDIR)\data.obj" \
	"$(INTDIR)\dlhandler_win32..obj" \
	"$(INTDIR)\dmol.obj" \
	"$(INTDIR)\DynamicOptions.obj" \
	"$(INTDIR)\generic.obj" \
	"$(INTDIR)\grid.obj" \
	"$(INTDIR)\hin.obj" \
	"$(INTDIR)\matrix.obj" \
	"$(INTDIR)\matrix3x3.obj" \
	"$(INTDIR)\mdlformat.obj" \
	"$(INTDIR)\mol.obj" \
	"$(INTDIR)\mol2format.obj" \
	"$(INTDIR)\molchrg.obj" \
	"$(INTDIR)\obconversion.obj" \
	"$(INTDIR)\oberror.obj" \
	"$(INTDIR)\OBGUI.obj" \
	"$(INTDIR)\OBGUIDlg.obj" \
	"$(INTDIR)\obutil.obj" \
	"$(INTDIR)\parsmart.obj" \
	"$(INTDIR)\patty.obj" \
	"$(INTDIR)\phmodel.obj" \
	"$(INTDIR)\rand.obj" \
	"$(INTDIR)\report.obj" \
	"$(INTDIR)\residue.obj" \
	"$(INTDIR)\ring.obj" \
	"$(INTDIR)\rotor.obj" \
	"$(INTDIR)\rxnformat.obj" \
	"$(INTDIR)\smilesformat.obj" \
	"$(INTDIR)\StdAfx.obj" \
	"$(INTDIR)\tinker.obj" \
	"$(INTDIR)\tokenst.obj" \
	"$(INTDIR)\transform.obj" \
	"$(INTDIR)\turbomoleformat.obj" \
	"$(INTDIR)\typer.obj" \
	"$(INTDIR)\vector3.obj" \
	"$(INTDIR)\OBGUI.res"

"$(OUTDIR)\OBGUIs.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("OBGUIs.dep")
!INCLUDE "OBGUIs.dep"
!ELSE 
!MESSAGE Warning: cannot find "OBGUIs.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "OBGUIs - Win32 Release" || "$(CFG)" == "OBGUIs - Win32 Debug"
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

SOURCE=..\..\src\bgf.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\bgf.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\bgf.obj"	"$(INTDIR)\bgf.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\binary.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\binary.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\binary.obj"	"$(INTDIR)\binary.sbr" : $(SOURCE) "$(INTDIR)"
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

SOURCE=..\..\src\box.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\box.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\box.obj"	"$(INTDIR)\box.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\c3d.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\c3d.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\c3d.obj"	"$(INTDIR)\c3d.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\cache.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\cache.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\cache.obj"	"$(INTDIR)\cache.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\car.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\car.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\car.obj"	"$(INTDIR)\car.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\ccc.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\ccc.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\ccc.obj"	"$(INTDIR)\ccc.sbr" : $(SOURCE) "$(INTDIR)"
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

SOURCE=..\..\src\chdrw.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\chdrw.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\chdrw.obj"	"$(INTDIR)\chdrw.sbr" : $(SOURCE) "$(INTDIR)"
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

SOURCE=..\..\src\formats\cml.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\cml.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\cml.obj"	"$(INTDIR)\cml.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\cmlformat.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\cmlformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\cmlformat.obj"	"$(INTDIR)\cmlformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\csr.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\csr.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\csr.obj"	"$(INTDIR)\csr.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\cssr.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\cssr.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\cssr.obj"	"$(INTDIR)\cssr.sbr" : $(SOURCE) "$(INTDIR)"
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

SOURCE=..\..\src\dlhandler_win32..cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\dlhandler_win32..obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\dlhandler_win32..obj"	"$(INTDIR)\dlhandler_win32..sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\dmol.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\dmol.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\dmol.obj"	"$(INTDIR)\dmol.sbr" : $(SOURCE) "$(INTDIR)"
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

SOURCE=..\..\src\hin.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\hin.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\hin.obj"	"$(INTDIR)\hin.sbr" : $(SOURCE) "$(INTDIR)"
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
	$(RSC) /l 0x809 /fo"$(INTDIR)\OBGUI.res" /i "\My Documents\MSVC\OpenBabel Ultimate\windows\OBGUI" /d "NDEBUG" $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\OBGUI.res" : $(SOURCE) "$(INTDIR)"
	$(RSC) /l 0x809 /fo"$(INTDIR)\OBGUI.res" /i "\My Documents\MSVC\OpenBabel Ultimate\windows\OBGUI" /d "_DEBUG" $(SOURCE)


!ENDIF 

SOURCE=..\OBGUI\OBGUIDlg.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\OBGUIDlg.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\OBGUIDlg.obj"	"$(INTDIR)\OBGUIDlg.sbr" : $(SOURCE) "$(INTDIR)"
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

SOURCE=..\..\src\rand.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\rand.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\rand.obj"	"$(INTDIR)\rand.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\report.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\report.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\report.obj"	"$(INTDIR)\report.sbr" : $(SOURCE) "$(INTDIR)"
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

SOURCE=..\..\src\tinker.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\tinker.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\tinker.obj"	"$(INTDIR)\tinker.sbr" : $(SOURCE) "$(INTDIR)"
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

SOURCE=..\..\src\formats\turbomoleformat.cpp

!IF  "$(CFG)" == "OBGUIs - Win32 Release"


"$(INTDIR)\turbomoleformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBGUIs - Win32 Debug"


"$(INTDIR)\turbomoleformat.obj"	"$(INTDIR)\turbomoleformat.sbr" : $(SOURCE) "$(INTDIR)"
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


!ENDIF 

