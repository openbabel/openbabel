# Microsoft Developer Studio Generated NMAKE File, Based on OBGUI.dsp
!IF "$(CFG)" == ""
CFG=OBGUI - Win32 Debug
!MESSAGE No configuration specified. Defaulting to OBGUI - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "OBGUI - Win32 Release" && "$(CFG)" != "OBGUI - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "OBGUI.mak" CFG="OBGUI - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "OBGUI - Win32 Release" (based on "Win32 (x86) Application")
!MESSAGE "OBGUI - Win32 Debug" (based on "Win32 (x86) Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "OBGUI - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release

!IF "$(RECURSE)" == "0" 

ALL : ".\OBGUI.exe"

!ELSE 

ALL : "OBConv - Win32 Release" ".\OBGUI.exe"

!ENDIF 

!IF "$(RECURSE)" == "1" 
CLEAN :"OBConv - Win32 ReleaseCLEAN" 
!ELSE 
CLEAN :
!ENDIF 
	-@erase "$(INTDIR)\DynamicOptions.obj"
	-@erase "$(INTDIR)\OBGUI.obj"
	-@erase "$(INTDIR)\OBGUI.res"
	-@erase "$(INTDIR)\OBGUIDlg.obj"
	-@erase "$(INTDIR)\StdAfx.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase ".\OBGUI.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MD /W3 /GR /GX /I "..\..\src" /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_AFXDLL" /D "GUI" /D "USING_DYNAMIC_LIBS" /Fp"$(INTDIR)\OBGUI.pch" /Yu"stdafx.h" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

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
RSC_PROJ=/l 0x809 /fo"$(INTDIR)\OBGUI.res" /d "NDEBUG" /d "_AFXDLL" 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\OBGUI.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=obconv.lib /nologo /subsystem:windows /incremental:no /pdb:"$(OUTDIR)\OBGUI.pdb" /machine:I386 /out:"OBGUI.exe" /libpath:"..\OBConv\Release\\" 
LINK32_OBJS= \
	"$(INTDIR)\DynamicOptions.obj" \
	"$(INTDIR)\OBGUI.obj" \
	"$(INTDIR)\OBGUIDlg.obj" \
	"$(INTDIR)\StdAfx.obj" \
	"$(INTDIR)\OBGUI.res" \
	"..\OBConv\Release\OBConv.lib"

".\OBGUI.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

SOURCE="$(InputPath)"
PostBuild_Desc=Copy obconv.dll, obdll.dll and obformats.obf
DS_POSTBUILD_DEP=$(INTDIR)\postbld.dep

ALL : $(DS_POSTBUILD_DEP)

$(DS_POSTBUILD_DEP) : "OBConv - Win32 Release" ".\OBGUI.exe"
   Copy  ..\obformats2\obformats2.obf . /Y
	Copy  ..\obconv\obconv.dll . /Y
	Copy  ..\obdll\obdll.dll . /Y
	echo Helper for Post-build step > "$(DS_POSTBUILD_DEP)"

!ELSEIF  "$(CFG)" == "OBGUI - Win32 Debug"

OUTDIR=.\Debug
INTDIR=.\Debug
# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

!IF "$(RECURSE)" == "0" 

ALL : "$(OUTDIR)\OBGUI.exe" "$(OUTDIR)\OBGUI.bsc"

!ELSE 

ALL : "OBConv - Win32 Debug" "$(OUTDIR)\OBGUI.exe" "$(OUTDIR)\OBGUI.bsc"

!ENDIF 

!IF "$(RECURSE)" == "1" 
CLEAN :"OBConv - Win32 DebugCLEAN" 
!ELSE 
CLEAN :
!ENDIF 
	-@erase "$(INTDIR)\DynamicOptions.obj"
	-@erase "$(INTDIR)\DynamicOptions.sbr"
	-@erase "$(INTDIR)\OBGUI.obj"
	-@erase "$(INTDIR)\OBGUI.res"
	-@erase "$(INTDIR)\OBGUI.sbr"
	-@erase "$(INTDIR)\OBGUIDlg.obj"
	-@erase "$(INTDIR)\OBGUIDlg.sbr"
	-@erase "$(INTDIR)\StdAfx.obj"
	-@erase "$(INTDIR)\StdAfx.sbr"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vc60.pdb"
	-@erase "$(OUTDIR)\OBGUI.bsc"
	-@erase "$(OUTDIR)\OBGUI.exe"
	-@erase "$(OUTDIR)\OBGUI.ilk"
	-@erase "$(OUTDIR)\OBGUI.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MDd /W3 /Gm /Gi /GR /GX /ZI /Od /I "..\..\src" /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_AFXDLL" /D "_MBCS" /D "GUI" /D "USING_DYNAMIC_LIBS" /FR"$(INTDIR)\\" /Fp"$(INTDIR)\OBGUI.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

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
RSC_PROJ=/l 0x809 /fo"$(INTDIR)\OBGUI.res" /d "_DEBUG" /d "_AFXDLL" 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\OBGUI.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\DynamicOptions.sbr" \
	"$(INTDIR)\OBGUI.sbr" \
	"$(INTDIR)\OBGUIDlg.sbr" \
	"$(INTDIR)\StdAfx.sbr"

"$(OUTDIR)\OBGUI.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
LINK32_FLAGS=obconv.lib /nologo /subsystem:windows /incremental:yes /pdb:"$(OUTDIR)\OBGUI.pdb" /debug /machine:I386 /out:"$(OUTDIR)\OBGUI.exe" /pdbtype:sept /libpath:"..\OBConv\debug\\" 
LINK32_OBJS= \
	"$(INTDIR)\DynamicOptions.obj" \
	"$(INTDIR)\OBGUI.obj" \
	"$(INTDIR)\OBGUIDlg.obj" \
	"$(INTDIR)\StdAfx.obj" \
	"$(INTDIR)\OBGUI.res" \
	"..\OBConv\Debug\OBConv.lib"

"$(OUTDIR)\OBGUI.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

SOURCE="$(InputPath)"
PostBuild_Desc=Copy obconv.dll, obdll.dll and obformats.obf
DS_POSTBUILD_DEP=$(INTDIR)\postbld.dep

ALL : $(DS_POSTBUILD_DEP)

# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

$(DS_POSTBUILD_DEP) : "OBConv - Win32 Debug" "$(OUTDIR)\OBGUI.exe" "$(OUTDIR)\OBGUI.bsc"
   Copy  ..\obconv\debug\obconv.dll  .\debug  /Y
	Copy  ..\obformats2\debug\obformats2D.obf .\debug /Y
	Copy  ..\obdll\debug\obdll.dll  .\debug  /Y
	echo Helper for Post-build step > "$(DS_POSTBUILD_DEP)"

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("OBGUI.dep")
!INCLUDE "OBGUI.dep"
!ELSE 
!MESSAGE Warning: cannot find "OBGUI.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "OBGUI - Win32 Release" || "$(CFG)" == "OBGUI - Win32 Debug"
SOURCE=.\DynamicOptions.cpp

!IF  "$(CFG)" == "OBGUI - Win32 Release"

CPP_SWITCHES=/nologo /MD /W3 /GR /GX /I "..\..\src" /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_AFXDLL" /D "GUI" /D "USING_DYNAMIC_LIBS" /Fp"$(INTDIR)\OBGUI.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

"$(INTDIR)\DynamicOptions.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) @<<
  $(CPP_SWITCHES) $(SOURCE)
<<


!ELSEIF  "$(CFG)" == "OBGUI - Win32 Debug"

CPP_SWITCHES=/nologo /MDd /W3 /Gm /Gi /GR /GX /ZI /Od /I "..\..\src" /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_AFXDLL" /D "_MBCS" /D "GUI" /D "USING_DYNAMIC_LIBS" /FR"$(INTDIR)\\" /Fp"$(INTDIR)\OBGUI.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

"$(INTDIR)\DynamicOptions.obj"	"$(INTDIR)\DynamicOptions.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) @<<
  $(CPP_SWITCHES) $(SOURCE)
<<


!ENDIF 

SOURCE=.\OBGUI.cpp

!IF  "$(CFG)" == "OBGUI - Win32 Release"

CPP_SWITCHES=/nologo /MD /W3 /GR /GX /I "..\..\src" /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_AFXDLL" /D "GUI" /D "USING_DYNAMIC_LIBS" /Fp"$(INTDIR)\OBGUI.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

"$(INTDIR)\OBGUI.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) @<<
  $(CPP_SWITCHES) $(SOURCE)
<<


!ELSEIF  "$(CFG)" == "OBGUI - Win32 Debug"

CPP_SWITCHES=/nologo /MDd /W3 /Gm /Gi /GR /GX /ZI /Od /I "..\..\src" /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_AFXDLL" /D "_MBCS" /D "GUI" /D "USING_DYNAMIC_LIBS" /FR"$(INTDIR)\\" /Fp"$(INTDIR)\OBGUI.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

"$(INTDIR)\OBGUI.obj"	"$(INTDIR)\OBGUI.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) @<<
  $(CPP_SWITCHES) $(SOURCE)
<<


!ENDIF 

SOURCE=.\OBGUI.rc

"$(INTDIR)\OBGUI.res" : $(SOURCE) "$(INTDIR)"
	$(RSC) $(RSC_PROJ) $(SOURCE)


SOURCE=.\OBGUIDlg.cpp

!IF  "$(CFG)" == "OBGUI - Win32 Release"

CPP_SWITCHES=/nologo /MD /W3 /GR /GX /I "..\..\src" /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_AFXDLL" /D "GUI" /D "USING_DYNAMIC_LIBS" /Fp"$(INTDIR)\OBGUI.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

"$(INTDIR)\OBGUIDlg.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) @<<
  $(CPP_SWITCHES) $(SOURCE)
<<


!ELSEIF  "$(CFG)" == "OBGUI - Win32 Debug"

CPP_SWITCHES=/nologo /MDd /W3 /Gm /Gi /GR /GX /ZI /Od /I "..\..\src" /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_AFXDLL" /D "_MBCS" /D "GUI" /D "USING_DYNAMIC_LIBS" /FR"$(INTDIR)\\" /Fp"$(INTDIR)\OBGUI.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

"$(INTDIR)\OBGUIDlg.obj"	"$(INTDIR)\OBGUIDlg.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) @<<
  $(CPP_SWITCHES) $(SOURCE)
<<


!ENDIF 

SOURCE=.\StdAfx.cpp

!IF  "$(CFG)" == "OBGUI - Win32 Release"

CPP_SWITCHES=/nologo /MD /W3 /GR /GX /I "..\..\src" /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_AFXDLL" /D "GUI" /D "USING_DYNAMIC_LIBS" /Fp"$(INTDIR)\OBGUI.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

"$(INTDIR)\StdAfx.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) @<<
  $(CPP_SWITCHES) $(SOURCE)
<<


!ELSEIF  "$(CFG)" == "OBGUI - Win32 Debug"

CPP_SWITCHES=/nologo /MDd /W3 /Gm /Gi /GR /GX /ZI /Od /I "..\..\src" /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_AFXDLL" /D "_MBCS" /D "GUI" /D "USING_DYNAMIC_LIBS" /FR"$(INTDIR)\\" /Fp"$(INTDIR)\OBGUI.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

"$(INTDIR)\StdAfx.obj"	"$(INTDIR)\StdAfx.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) @<<
  $(CPP_SWITCHES) $(SOURCE)
<<


!ENDIF 

!IF  "$(CFG)" == "OBGUI - Win32 Release"

"OBConv - Win32 Release" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\windows\OBConv"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBConv.mak" CFG="OBConv - Win32 Release" 
   cd "..\OBGUI"

"OBConv - Win32 ReleaseCLEAN" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\windows\OBConv"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBConv.mak" CFG="OBConv - Win32 Release" RECURSE=1 CLEAN 
   cd "..\OBGUI"

!ELSEIF  "$(CFG)" == "OBGUI - Win32 Debug"

"OBConv - Win32 Debug" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\windows\OBConv"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBConv.mak" CFG="OBConv - Win32 Debug" 
   cd "..\OBGUI"

"OBConv - Win32 DebugCLEAN" : 
   cd "\My Documents\MSVC\OpenBabel Ultimate\windows\OBConv"
   $(MAKE) /$(MAKEFLAGS) /F ".\OBConv.mak" CFG="OBConv - Win32 Debug" RECURSE=1 CLEAN 
   cd "..\OBGUI"

!ENDIF 


!ENDIF 

