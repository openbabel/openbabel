# Microsoft Developer Studio Generated NMAKE File, Based on ProgramTest2.dsp
!IF "$(CFG)" == ""
CFG=ProgramTest2 - Win32 Debug
!MESSAGE No configuration specified. Defaulting to ProgramTest2 - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "ProgramTest2 - Win32 Release" && "$(CFG)" != "ProgramTest2 - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "ProgramTest2.mak" CFG="ProgramTest2 - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "ProgramTest2 - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "ProgramTest2 - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "ProgramTest2 - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release

ALL : ".\ProgramTest2.exe"


CLEAN :
	-@erase "$(INTDIR)\prog2.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase ".\ProgramTest2.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MD /W3 /GX /O2 /I "..\..\src" /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /Fp"$(INTDIR)\ProgramTest2.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

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

RSC=rc.exe
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\ProgramTest2.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib obconv.lib /nologo /subsystem:console /incremental:no /pdb:"$(OUTDIR)\ProgramTest2.pdb" /machine:I386 /out:"ProgramTest2.exe" /libpath:"..\obconv\release" 
LINK32_OBJS= \
	"$(INTDIR)\prog2.obj"

".\ProgramTest2.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "ProgramTest2 - Win32 Debug"

OUTDIR=.\Debug
INTDIR=.\Debug
# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

ALL : "$(OUTDIR)\ProgramTest2.exe" "$(OUTDIR)\ProgramTest2.bsc"


CLEAN :
	-@erase "$(INTDIR)\prog2.obj"
	-@erase "$(INTDIR)\prog2.sbr"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vc60.pdb"
	-@erase "$(OUTDIR)\ProgramTest2.bsc"
	-@erase "$(OUTDIR)\ProgramTest2.exe"
	-@erase "$(OUTDIR)\ProgramTest2.ilk"
	-@erase "$(OUTDIR)\ProgramTest2.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MDd /W3 /Gm /GR /GX /ZI /Od /I "..\..\src" /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /FR"$(INTDIR)\\" /Fp"$(INTDIR)\ProgramTest2.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

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

RSC=rc.exe
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\ProgramTest2.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\prog2.sbr"

"$(OUTDIR)\ProgramTest2.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib obconv.lib /nologo /subsystem:console /incremental:yes /pdb:"$(OUTDIR)\ProgramTest2.pdb" /debug /machine:I386 /out:"$(OUTDIR)\ProgramTest2.exe" /pdbtype:sept /libpath:"..\obconv\debug" 
LINK32_OBJS= \
	"$(INTDIR)\prog2.obj"

"$(OUTDIR)\ProgramTest2.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

SOURCE="$(InputPath)"
PostBuild_Desc=Copy obdll.dll, obconv.dll, obformats.obf
DS_POSTBUILD_DEP=$(INTDIR)\postbld.dep

ALL : $(DS_POSTBUILD_DEP)

# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

$(DS_POSTBUILD_DEP) : "$(OUTDIR)\ProgramTest2.exe" "$(OUTDIR)\ProgramTest2.bsc"
   copy  ..\obdll\debug\obdll.dll  .\debug  /Y
	copy  ..\obconv\debug\obconv.dll  .\debug  /Y
	copy  ..\obformats\debug\obformats.obf  .  /Y
	echo Helper for Post-build step > "$(DS_POSTBUILD_DEP)"

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("ProgramTest2.dep")
!INCLUDE "ProgramTest2.dep"
!ELSE 
!MESSAGE Warning: cannot find "ProgramTest2.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "ProgramTest2 - Win32 Release" || "$(CFG)" == "ProgramTest2 - Win32 Debug"
SOURCE=.\prog2.cpp

!IF  "$(CFG)" == "ProgramTest2 - Win32 Release"


"$(INTDIR)\prog2.obj" : $(SOURCE) "$(INTDIR)"


!ELSEIF  "$(CFG)" == "ProgramTest2 - Win32 Debug"


"$(INTDIR)\prog2.obj"	"$(INTDIR)\prog2.sbr" : $(SOURCE) "$(INTDIR)"


!ENDIF 


!ENDIF 

