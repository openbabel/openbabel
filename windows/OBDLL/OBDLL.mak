# Microsoft Developer Studio Generated NMAKE File, Based on OBDLL.dsp
!IF "$(CFG)" == ""
CFG=OBDLL - Win32 Debug
!MESSAGE No configuration specified. Defaulting to OBDLL - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "OBDLL - Win32 Release" && "$(CFG)" != "OBDLL - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "OBDLL.mak" CFG="OBDLL - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "OBDLL - Win32 Release" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE "OBDLL - Win32 Debug" (based on "Win32 (x86) Dynamic-Link Library")
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

!IF  "$(CFG)" == "OBDLL - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release

ALL : ".\OBDLL.dll"


CLEAN :
	-@erase "$(INTDIR)\atom.obj"
	-@erase "$(INTDIR)\base.obj"
	-@erase "$(INTDIR)\bitvec.obj"
	-@erase "$(INTDIR)\bond.obj"
	-@erase "$(INTDIR)\bondtyper.obj"
	-@erase "$(INTDIR)\chains.obj"
	-@erase "$(INTDIR)\chiral.obj"
	-@erase "$(INTDIR)\data.obj"
	-@erase "$(INTDIR)\fingerprint.obj"
	-@erase "$(INTDIR)\generic.obj"
	-@erase "$(INTDIR)\grid.obj"
	-@erase "$(INTDIR)\kekulize.obj"
	-@erase "$(INTDIR)\matrix.obj"
	-@erase "$(INTDIR)\matrix3x3.obj"
	-@erase "$(INTDIR)\mol.obj"
	-@erase "$(INTDIR)\molchrg.obj"
	-@erase "$(INTDIR)\oberror.obj"
	-@erase "$(INTDIR)\obiter.obj"
	-@erase "$(INTDIR)\obutil.obj"
	-@erase "$(INTDIR)\parsmart.obj"
	-@erase "$(INTDIR)\patty.obj"
	-@erase "$(INTDIR)\phmodel.obj"
	-@erase "$(INTDIR)\rand.obj"
	-@erase "$(INTDIR)\residue.obj"
	-@erase "$(INTDIR)\ring.obj"
	-@erase "$(INTDIR)\rotamer.obj"
	-@erase "$(INTDIR)\rotor.obj"
	-@erase "$(INTDIR)\tokenst.obj"
	-@erase "$(INTDIR)\transform.obj"
	-@erase "$(INTDIR)\typer.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vector3.obj"
	-@erase "$(OUTDIR)\OBDLL.exp"
	-@erase "$(OUTDIR)\OBDLL.lib"
	-@erase ".\OBDLL.dll"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP_PROJ=/nologo /MD /W3 /GR /GX /I "..\math ..\src" /I "..\..\src" /I ".." /I "../../data" /D "NDEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "__KCC" /D "USING_DYNAMIC_LIBS" /D "OBDLL_EXPORTS" /D "HAVE_CONFIG_H" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 
MTL_PROJ=/nologo /D "NDEBUG" /mktyplib203 /win32 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\OBDLL.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib /nologo /dll /incremental:no /pdb:"$(OUTDIR)\OBDLL.pdb" /machine:I386 /out:"OBDLL.dll" /implib:"$(OUTDIR)\OBDLL.lib" 
LINK32_OBJS= \
	"$(INTDIR)\atom.obj" \
	"$(INTDIR)\base.obj" \
	"$(INTDIR)\bitvec.obj" \
	"$(INTDIR)\bond.obj" \
	"$(INTDIR)\bondtyper.obj" \
	"$(INTDIR)\chains.obj" \
	"$(INTDIR)\chiral.obj" \
	"$(INTDIR)\data.obj" \
	"$(INTDIR)\fingerprint.obj" \
	"$(INTDIR)\generic.obj" \
	"$(INTDIR)\grid.obj" \
	"$(INTDIR)\kekulize.obj" \
	"$(INTDIR)\matrix.obj" \
	"$(INTDIR)\matrix3x3.obj" \
	"$(INTDIR)\mol.obj" \
	"$(INTDIR)\molchrg.obj" \
	"$(INTDIR)\oberror.obj" \
	"$(INTDIR)\obiter.obj" \
	"$(INTDIR)\obutil.obj" \
	"$(INTDIR)\parsmart.obj" \
	"$(INTDIR)\patty.obj" \
	"$(INTDIR)\phmodel.obj" \
	"$(INTDIR)\rand.obj" \
	"$(INTDIR)\residue.obj" \
	"$(INTDIR)\ring.obj" \
	"$(INTDIR)\rotamer.obj" \
	"$(INTDIR)\rotor.obj" \
	"$(INTDIR)\tokenst.obj" \
	"$(INTDIR)\transform.obj" \
	"$(INTDIR)\typer.obj" \
	"$(INTDIR)\vector3.obj"

".\OBDLL.dll" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"

OUTDIR=.\Debug
INTDIR=.\Debug
# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

ALL : "$(OUTDIR)\OBDLL.dll" "$(OUTDIR)\OBDLL.bsc"


CLEAN :
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
	-@erase "$(INTDIR)\data.obj"
	-@erase "$(INTDIR)\data.sbr"
	-@erase "$(INTDIR)\fingerprint.obj"
	-@erase "$(INTDIR)\fingerprint.sbr"
	-@erase "$(INTDIR)\generic.obj"
	-@erase "$(INTDIR)\generic.sbr"
	-@erase "$(INTDIR)\grid.obj"
	-@erase "$(INTDIR)\grid.sbr"
	-@erase "$(INTDIR)\kekulize.obj"
	-@erase "$(INTDIR)\kekulize.sbr"
	-@erase "$(INTDIR)\matrix.obj"
	-@erase "$(INTDIR)\matrix.sbr"
	-@erase "$(INTDIR)\matrix3x3.obj"
	-@erase "$(INTDIR)\matrix3x3.sbr"
	-@erase "$(INTDIR)\mol.obj"
	-@erase "$(INTDIR)\mol.sbr"
	-@erase "$(INTDIR)\molchrg.obj"
	-@erase "$(INTDIR)\molchrg.sbr"
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
	-@erase "$(INTDIR)\phmodel.obj"
	-@erase "$(INTDIR)\phmodel.sbr"
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
	-@erase "$(OUTDIR)\OBDLL.bsc"
	-@erase "$(OUTDIR)\OBDLL.dll"
	-@erase "$(OUTDIR)\OBDLL.exp"
	-@erase "$(OUTDIR)\OBDLL.ilk"
	-@erase "$(OUTDIR)\OBDLL.lib"
	-@erase "$(OUTDIR)\OBDLL.pdb"
	-@erase ".\OBDLL.map"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP_PROJ=/nologo /MDd /Gm /GR /GX /ZI /Od /I "..\math ..\src" /I "..\..\src" /I ".." /I "../../data" /D "_DEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "__KCC" /D "USING_DYNAMIC_LIBS" /D "OBDLL_EXPORTS" /D "HAVE_CONFIG_H" /FR"$(INTDIR)\\" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 
MTL_PROJ=/nologo /D "_DEBUG" /mktyplib203 /win32 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\OBDLL.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\atom.sbr" \
	"$(INTDIR)\base.sbr" \
	"$(INTDIR)\bitvec.sbr" \
	"$(INTDIR)\bond.sbr" \
	"$(INTDIR)\bondtyper.sbr" \
	"$(INTDIR)\chains.sbr" \
	"$(INTDIR)\chiral.sbr" \
	"$(INTDIR)\data.sbr" \
	"$(INTDIR)\fingerprint.sbr" \
	"$(INTDIR)\generic.sbr" \
	"$(INTDIR)\grid.sbr" \
	"$(INTDIR)\kekulize.sbr" \
	"$(INTDIR)\matrix.sbr" \
	"$(INTDIR)\matrix3x3.sbr" \
	"$(INTDIR)\mol.sbr" \
	"$(INTDIR)\molchrg.sbr" \
	"$(INTDIR)\oberror.sbr" \
	"$(INTDIR)\obiter.sbr" \
	"$(INTDIR)\obutil.sbr" \
	"$(INTDIR)\parsmart.sbr" \
	"$(INTDIR)\patty.sbr" \
	"$(INTDIR)\phmodel.sbr" \
	"$(INTDIR)\rand.sbr" \
	"$(INTDIR)\residue.sbr" \
	"$(INTDIR)\ring.sbr" \
	"$(INTDIR)\rotamer.sbr" \
	"$(INTDIR)\rotor.sbr" \
	"$(INTDIR)\tokenst.sbr" \
	"$(INTDIR)\transform.sbr" \
	"$(INTDIR)\typer.sbr" \
	"$(INTDIR)\vector3.sbr"

"$(OUTDIR)\OBDLL.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib /nologo /dll /incremental:yes /pdb:"$(OUTDIR)\OBDLL.pdb" /map:"OBDLL.map" /debug /machine:I386 /out:"$(OUTDIR)\OBDLL.dll" /implib:"$(OUTDIR)\OBDLL.lib" /pdbtype:sept 
LINK32_OBJS= \
	"$(INTDIR)\atom.obj" \
	"$(INTDIR)\base.obj" \
	"$(INTDIR)\bitvec.obj" \
	"$(INTDIR)\bond.obj" \
	"$(INTDIR)\bondtyper.obj" \
	"$(INTDIR)\chains.obj" \
	"$(INTDIR)\chiral.obj" \
	"$(INTDIR)\data.obj" \
	"$(INTDIR)\fingerprint.obj" \
	"$(INTDIR)\generic.obj" \
	"$(INTDIR)\grid.obj" \
	"$(INTDIR)\kekulize.obj" \
	"$(INTDIR)\matrix.obj" \
	"$(INTDIR)\matrix3x3.obj" \
	"$(INTDIR)\mol.obj" \
	"$(INTDIR)\molchrg.obj" \
	"$(INTDIR)\oberror.obj" \
	"$(INTDIR)\obiter.obj" \
	"$(INTDIR)\obutil.obj" \
	"$(INTDIR)\parsmart.obj" \
	"$(INTDIR)\patty.obj" \
	"$(INTDIR)\phmodel.obj" \
	"$(INTDIR)\rand.obj" \
	"$(INTDIR)\residue.obj" \
	"$(INTDIR)\ring.obj" \
	"$(INTDIR)\rotamer.obj" \
	"$(INTDIR)\rotor.obj" \
	"$(INTDIR)\tokenst.obj" \
	"$(INTDIR)\transform.obj" \
	"$(INTDIR)\typer.obj" \
	"$(INTDIR)\vector3.obj"

"$(OUTDIR)\OBDLL.dll" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
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
!IF EXISTS("OBDLL.dep")
!INCLUDE "OBDLL.dep"
!ELSE 
!MESSAGE Warning: cannot find "OBDLL.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "OBDLL - Win32 Release" || "$(CFG)" == "OBDLL - Win32 Debug"
SOURCE=..\..\src\atom.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\atom.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\atom.obj"	"$(INTDIR)\atom.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\base.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\base.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\base.obj"	"$(INTDIR)\base.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\bitvec.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\bitvec.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\bitvec.obj"	"$(INTDIR)\bitvec.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\bond.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\bond.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\bond.obj"	"$(INTDIR)\bond.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\bondtyper.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\bondtyper.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\bondtyper.obj"	"$(INTDIR)\bondtyper.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\chains.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\chains.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\chains.obj"	"$(INTDIR)\chains.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\chiral.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\chiral.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\chiral.obj"	"$(INTDIR)\chiral.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\data.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\data.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\data.obj"	"$(INTDIR)\data.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\fingerprint.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\fingerprint.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\fingerprint.obj"	"$(INTDIR)\fingerprint.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\generic.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\generic.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\generic.obj"	"$(INTDIR)\generic.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\grid.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\grid.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\grid.obj"	"$(INTDIR)\grid.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\kekulize.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\kekulize.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\kekulize.obj"	"$(INTDIR)\kekulize.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\matrix.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\matrix.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\matrix.obj"	"$(INTDIR)\matrix.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\math\matrix3x3.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\matrix3x3.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\matrix3x3.obj"	"$(INTDIR)\matrix3x3.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\mol.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\mol.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\mol.obj"	"$(INTDIR)\mol.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\molchrg.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\molchrg.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\molchrg.obj"	"$(INTDIR)\molchrg.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\oberror.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\oberror.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\oberror.obj"	"$(INTDIR)\oberror.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\obiter.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\obiter.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\obiter.obj"	"$(INTDIR)\obiter.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\obutil.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\obutil.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\obutil.obj"	"$(INTDIR)\obutil.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\parsmart.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\parsmart.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\parsmart.obj"	"$(INTDIR)\parsmart.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\patty.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\patty.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\patty.obj"	"$(INTDIR)\patty.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\phmodel.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\phmodel.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\phmodel.obj"	"$(INTDIR)\phmodel.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\rand.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\rand.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\rand.obj"	"$(INTDIR)\rand.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\residue.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\residue.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\residue.obj"	"$(INTDIR)\residue.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\ring.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\ring.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\ring.obj"	"$(INTDIR)\ring.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\rotamer.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\rotamer.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\rotamer.obj"	"$(INTDIR)\rotamer.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\rotor.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\rotor.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\rotor.obj"	"$(INTDIR)\rotor.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\tokenst.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\tokenst.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\tokenst.obj"	"$(INTDIR)\tokenst.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\transform.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\transform.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\transform.obj"	"$(INTDIR)\transform.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\typer.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\typer.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\typer.obj"	"$(INTDIR)\typer.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\math\vector3.cpp

!IF  "$(CFG)" == "OBDLL - Win32 Release"


"$(INTDIR)\vector3.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"


"$(INTDIR)\vector3.obj"	"$(INTDIR)\vector3.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 


!ENDIF 

