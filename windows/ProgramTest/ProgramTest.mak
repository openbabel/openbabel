# Microsoft Developer Studio Generated NMAKE File, Based on ProgramTest.dsp
!IF "$(CFG)" == ""
CFG=ProgramTest - Win32 Debug
!MESSAGE No configuration specified. Defaulting to ProgramTest - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "ProgramTest - Win32 Release" && "$(CFG)" != "ProgramTest - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "ProgramTest.mak" CFG="ProgramTest - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "ProgramTest - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "ProgramTest - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "ProgramTest - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release

!IF "$(RECURSE)" == "0" 

ALL : ".\ProgramTest.exe"

!ELSE 

ALL : "OBFormats - Win32 Release" "OBDLL - Win32 Release" "OBConv - Win32 Release" ".\ProgramTest.exe"

!ENDIF 

!IF "$(RECURSE)" == "1" 
CLEAN :"OBConv - Win32 ReleaseCLEAN" "OBDLL - Win32 ReleaseCLEAN" "OBFormats - Win32 ReleaseCLEAN" 
!ELSE 
CLEAN :
!ENDIF 
	-@erase "$(INTDIR)\prog1.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase ".\ProgramTest.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MD /W3 /GR /GX /O2 /I "..\..\src" /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /Fp"$(INTDIR)\ProgramTest.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\ProgramTest.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib obdll.lib obconv.lib /nologo /subsystem:console /incremental:no /pdb:"$(OUTDIR)\ProgramTest.pdb" /machine:I386 /out:"ProgramTest.exe" /libpath:"..\OBConv\release" /libpath:"..\OBDLL\release" 
LINK32_OBJS= \
	"$(INTDIR)\prog1.obj" \
	"..\OBConv\Release\OBConv.lib" \
	"..\OBDLL\Release\OBDLL.lib"

".\ProgramTest.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "ProgramTest - Win32 Debug"

OUTDIR=.\Debug
INTDIR=.\Debug
# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

!IF "$(RECURSE)" == "0" 

ALL : "$(OUTDIR)\ProgramTest.exe" "$(OUTDIR)\ProgramTest.bsc"

!ELSE 

ALL : "OBFormats - Win32 Debug" "OBDLL - Win32 Debug" "OBConv - Win32 Debug" "$(OUTDIR)\ProgramTest.exe" "$(OUTDIR)\ProgramTest.bsc"

!ENDIF 

!IF "$(RECURSE)" == "1" 
CLEAN :"OBConv - Win32 DebugCLEAN" "OBDLL - Win32 DebugCLEAN" "OBFormats - Win32 DebugCLEAN" 
!ELSE 
CLEAN :
!ENDIF 
	-@erase "$(INTDIR)\prog1.obj"
	-@erase "$(INTDIR)\prog1.sbr"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vc60.pdb"
	-@erase "$(OUTDIR)\ProgramTest.bsc"
	-@erase "$(OUTDIR)\ProgramTest.exe"
	-@erase "$(OUTDIR)\ProgramTest.ilk"
	-@erase "$(OUTDIR)\ProgramTest.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MDd /W3 /Gm /GR /GX /ZI /Od /I "..\..\src" /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /FR"$(INTDIR)\\" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)\ProgramTest.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\prog1.sbr"

"$(OUTDIR)\ProgramTest.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib obdll.lib obconv.lib /nologo /subsystem:console /incremental:yes /pdb:"$(OUTDIR)\ProgramTest.pdb" /debug /machine:I386 /out:"$(OUTDIR)\ProgramTest.exe" /pdbtype:sept /libpath:"..\OBConv\debug" /libpath:"..\OBDLL\debug" 
LINK32_OBJS= \
	"$(INTDIR)\prog1.obj" \
	"..\OBConv\Debug\OBConv.lib" \
	"..\OBDLL\Debug\OBDLL.lib"

"$(OUTDIR)\ProgramTest.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
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

$(DS_POSTBUILD_DEP) : "OBFormats - Win32 Debug" "OBDLL - Win32 Debug" "OBConv - Win32 Debug" "$(OUTDIR)\ProgramTest.exe" "$(OUTDIR)\ProgramTest.bsc"
   copy  ..\obdll\debug\obdll.dll  .\debug  /Y
	copy  ..\obconv\debug\obconv.dll  .\debug  /Y
	copy  ..\obformats\debug\obformats.obf  .  /Y
	echo Helper for Post-build step > "$(DS_POSTBUILD_DEP)"

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("ProgramTest.dep")
!INCLUDE "ProgramTest.dep"
!ELSE 
!MESSAGE Warning: cannot find "ProgramTest.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "ProgramTest - Win32 Release" || "$(CFG)" == "ProgramTest - Win32 Debug"
SOURCE=.\prog1.cpp

!IF  "$(CFG)" == "ProgramTest - Win32 Release"


"$(INTDIR)\prog1.obj" : $(SOURCE) "$(INTDIR)"


!ELSEIF  "$(CFG)" == "ProgramTest - Win32 Debug"


"$(INTDIR)\prog1.obj"	"$(INTDIR)\prog1.sbr" : $(SOURCE) "$(INTDIR)"


!ENDIF 

!IF  "$(CFG)" == "ProgramTest - Win32 Release"

"OBConv - Win32 Release" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\windows\OBConv"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBConv.mak" CFG="OBConv - Win32 Release" 
   cd "..\ProgramTest"

"OBConv - Win32 ReleaseCLEAN" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\windows\OBConv"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBConv.mak" CFG="OBConv - Win32 Release" RECURSE=1 CLEAN 
   cd "..\ProgramTest"

!ELSEIF  "$(CFG)" == "ProgramTest - Win32 Debug"

"OBConv - Win32 Debug" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\windows\OBConv"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBConv.mak" CFG="OBConv - Win32 Debug" 
   cd "..\ProgramTest"

"OBConv - Win32 DebugCLEAN" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\windows\OBConv"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBConv.mak" CFG="OBConv - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\ProgramTest"

!ENDIF 

!IF  "$(CFG)" == "ProgramTest - Win32 Release"

"OBDLL - Win32 Release" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\windows\OBDLL"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBDLL.mak" CFG="OBDLL - Win32 Release" 
   cd "..\ProgramTest"

"OBDLL - Win32 ReleaseCLEAN" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\windows\OBDLL"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBDLL.mak" CFG="OBDLL - Win32 Release" RECURSE=1 CLEAN 
   cd "..\ProgramTest"

!ELSEIF  "$(CFG)" == "ProgramTest - Win32 Debug"

"OBDLL - Win32 Debug" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\windows\OBDLL"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBDLL.mak" CFG="OBDLL - Win32 Debug" 
   cd "..\ProgramTest"

"OBDLL - Win32 DebugCLEAN" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\windows\OBDLL"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBDLL.mak" CFG="OBDLL - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\ProgramTest"

!ENDIF 

!IF  "$(CFG)" == "ProgramTest - Win32 Release"

"OBFormats - Win32 Release" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\windows\OBFormats"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBFormats.mak" CFG="OBFormats - Win32 Release" 
   cd "..\ProgramTest"

"OBFormats - Win32 ReleaseCLEAN" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\windows\OBFormats"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBFormats.mak" CFG="OBFormats - Win32 Release" RECURSE=1 CLEAN 
   cd "..\ProgramTest"

!ELSEIF  "$(CFG)" == "ProgramTest - Win32 Debug"

"OBFormats - Win32 Debug" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\windows\OBFormats"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBFormats.mak" CFG="OBFormats - Win32 Debug" 
   cd "..\ProgramTest"

"OBFormats - Win32 DebugCLEAN" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\windows\OBFormats"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBFormats.mak" CFG="OBFormats - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\ProgramTest"

!ENDIF 


!ENDIF 

