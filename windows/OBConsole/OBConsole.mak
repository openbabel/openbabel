# Microsoft Developer Studio Generated NMAKE File, Based on OBConsole.dsp
!IF "$(CFG)" == ""
CFG=OBConsole - Win32 Debug
!MESSAGE No configuration specified. Defaulting to OBConsole - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "OBConsole - Win32 Release" && "$(CFG)" != "OBConsole - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "OBConsole.mak" CFG="OBConsole - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "OBConsole - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "OBConsole - Win32 Debug" (based on "Win32 (x86) Console Application")
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

!IF  "$(CFG)" == "OBConsole - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release

!IF "$(RECURSE)" == "0" 

ALL : ".\babel.exe"

!ELSE 

ALL : "OBDLL - Win32 Release" "OBConv - Win32 Release" ".\babel.exe"

!ENDIF 

!IF "$(RECURSE)" == "1" 
CLEAN :"OBConv - Win32 ReleaseCLEAN" "OBDLL - Win32 ReleaseCLEAN" 
!ELSE 
CLEAN :
!ENDIF 
	-@erase "$(INTDIR)\main.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase ".\babel.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP_PROJ=/nologo /MD /W3 /GR /GX /O2 /I "..\OBConversion" /I ".." /D "NDEBUG" /D "WIN32" /D "_CONSOLE" /D "_MBCS" /D "USING_DYNAMIC_LIBS" /D "HAVE_CONFIG_H" /Fp"$(INTDIR)\OBConsole.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\OBConsole.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib obconv.lib /nologo /subsystem:console /incremental:no /pdb:"$(OUTDIR)\babel.pdb" /machine:I386 /out:"babel.exe" /libpath:"..\OBConv\Release" 
LINK32_OBJS= \
	"$(INTDIR)\main.obj" \
	"..\OBConv\Release\OBConv.lib" \
	"..\OBDLL\Release\OBDLL.lib"

".\babel.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

SOURCE="$(InputPath)"
PostBuild_Desc=Copy OBConv.dll, OBDLL.dll and OBFormats2.obf
DS_POSTBUILD_DEP=$(INTDIR)\postbld.dep

ALL : $(DS_POSTBUILD_DEP)

$(DS_POSTBUILD_DEP) : "OBDLL - Win32 Release" "OBConv - Win32 Release" ".\babel.exe"
   copy ..\OBConv\OBConv.dll  . /Y
	copy  ..\OBFormats2\OBFormats2.obf . /Y
	copy ..\OBDLL\OBDLL.dll  .  /Y
	echo Helper for Post-build step > "$(DS_POSTBUILD_DEP)"

!ELSEIF  "$(CFG)" == "OBConsole - Win32 Debug"

OUTDIR=.\Debug
INTDIR=.\Debug
# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

!IF "$(RECURSE)" == "0" 

ALL : "$(OUTDIR)\babel.exe" "$(OUTDIR)\OBConsole.bsc"

!ELSE 

ALL : "OBDLL - Win32 Debug" "OBConv - Win32 Debug" "$(OUTDIR)\babel.exe" "$(OUTDIR)\OBConsole.bsc"

!ENDIF 

!IF "$(RECURSE)" == "1" 
CLEAN :"OBConv - Win32 DebugCLEAN" "OBDLL - Win32 DebugCLEAN" 
!ELSE 
CLEAN :
!ENDIF 
	-@erase "$(INTDIR)\main.obj"
	-@erase "$(INTDIR)\main.sbr"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vc60.pdb"
	-@erase "$(OUTDIR)\babel.exe"
	-@erase "$(OUTDIR)\babel.ilk"
	-@erase "$(OUTDIR)\babel.pdb"
	-@erase "$(OUTDIR)\OBConsole.bsc"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP_PROJ=/nologo /MDd /W3 /Gm /GR /GX /ZI /Od /I "..\OBConversion" /I ".." /D "_DEBUG" /D "WIN32" /D "_CONSOLE" /D "_MBCS" /D "USING_DYNAMIC_LIBS" /D "HAVE_CONFIG_H" /FR"$(INTDIR)\\" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\OBConsole.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\main.sbr"

"$(OUTDIR)\OBConsole.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib obconv.lib /nologo /subsystem:console /incremental:yes /pdb:"$(OUTDIR)\babel.pdb" /debug /machine:I386 /out:"$(OUTDIR)\babel.exe" /pdbtype:sept /libpath:"..\OBConv\debug\\" 
LINK32_OBJS= \
	"$(INTDIR)\main.obj" \
	"..\OBConv\Debug\OBConv.lib" \
	"..\OBDLL\Debug\OBDLL.lib"

"$(OUTDIR)\babel.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

SOURCE="$(InputPath)"
PostBuild_Desc=Copy OBConv.dll, OBDLL.dll and OBFormats2.obf
DS_POSTBUILD_DEP=$(INTDIR)\postbld.dep

ALL : $(DS_POSTBUILD_DEP)

# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

$(DS_POSTBUILD_DEP) : "OBDLL - Win32 Debug" "OBConv - Win32 Debug" "$(OUTDIR)\babel.exe" "$(OUTDIR)\OBConsole.bsc"
   copy ..\OBConv\Debug\OBConv.dll  .\debug /Y
	copy  ..\OBFormats2\debug\OBFormats2D.obf .\debug /Y
	copy ..\OBDLL\debug\OBDLL.dll  .\debug  /Y
	echo Helper for Post-build step > "$(DS_POSTBUILD_DEP)"

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
!IF EXISTS("OBConsole.dep")
!INCLUDE "OBConsole.dep"
!ELSE 
!MESSAGE Warning: cannot find "OBConsole.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "OBConsole - Win32 Release" || "$(CFG)" == "OBConsole - Win32 Debug"
SOURCE=..\..\src\main.cpp

!IF  "$(CFG)" == "OBConsole - Win32 Release"


"$(INTDIR)\main.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBConsole - Win32 Debug"


"$(INTDIR)\main.obj"	"$(INTDIR)\main.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

!IF  "$(CFG)" == "OBConsole - Win32 Release"

"OBConv - Win32 Release" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\CVSforOB2\openbabel\windows\OBConv"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBConv.mak" CFG="OBConv - Win32 Release" 
   cd "..\OBConsole"

"OBConv - Win32 ReleaseCLEAN" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\CVSforOB2\openbabel\windows\OBConv"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBConv.mak" CFG="OBConv - Win32 Release" RECURSE=1 CLEAN 
   cd "..\OBConsole"

!ELSEIF  "$(CFG)" == "OBConsole - Win32 Debug"

"OBConv - Win32 Debug" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\CVSforOB2\openbabel\windows\OBConv"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBConv.mak" CFG="OBConv - Win32 Debug" 
   cd "..\OBConsole"

"OBConv - Win32 DebugCLEAN" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\CVSforOB2\openbabel\windows\OBConv"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBConv.mak" CFG="OBConv - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\OBConsole"

!ENDIF 

!IF  "$(CFG)" == "OBConsole - Win32 Release"

"OBDLL - Win32 Release" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\CVSforOB2\openbabel\windows\OBDLL"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBDLL.mak" CFG="OBDLL - Win32 Release" 
   cd "..\OBConsole"

"OBDLL - Win32 ReleaseCLEAN" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\CVSforOB2\openbabel\windows\OBDLL"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBDLL.mak" CFG="OBDLL - Win32 Release" RECURSE=1 CLEAN 
   cd "..\OBConsole"

!ELSEIF  "$(CFG)" == "OBConsole - Win32 Debug"

"OBDLL - Win32 Debug" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\CVSforOB2\openbabel\windows\OBDLL"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBDLL.mak" CFG="OBDLL - Win32 Debug" 
   cd "..\OBConsole"

"OBDLL - Win32 DebugCLEAN" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\CVSforOB2\openbabel\windows\OBDLL"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBDLL.mak" CFG="OBDLL - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\OBConsole"

!ENDIF 


!ENDIF 

