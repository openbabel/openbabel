# Microsoft Developer Studio Generated NMAKE File, Based on OBExtraFormats.dsp
!IF "$(CFG)" == ""
CFG=OBExtraFormats - Win32 Debug
!MESSAGE No configuration specified. Defaulting to OBExtraFormats - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "OBExtraFormats - Win32 Release" && "$(CFG)" != "OBExtraFormats - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "OBExtraFormats.mak" CFG="OBExtraFormats - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "OBExtraFormats - Win32 Release" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE "OBExtraFormats - Win32 Debug" (based on "Win32 (x86) Dynamic-Link Library")
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

!IF  "$(CFG)" == "OBExtraFormats - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release

ALL : ".\OBExtra.obf"


CLEAN :
	-@erase "$(INTDIR)\fastsearchformat.obj"
	-@erase "$(INTDIR)\finger2.obj"
	-@erase "$(INTDIR)\finger3.obj"
	-@erase "$(INTDIR)\fingerprintformat.obj"
	-@erase "$(INTDIR)\inchiformat.obj"
	-@erase "$(INTDIR)\mpdformat.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\OBExtra.exp"
	-@erase ".\OBExtra.obf"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP_PROJ=/nologo /MD /W3 /GR /GX /O2 /I "..\..\data" /I ".." /I "..\..\src" /D "NDEBUG" /D "USING_DYNAMIC_LIBS" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "INCHI_LINK_AS_DLL" /D "USING_OBDLL" /D "HAVE_CONFIG_H" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 
MTL_PROJ=/nologo /D "NDEBUG" /mktyplib203 /win32 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\OBExtraFormats.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib obconv.lib obdll.lib libinchi.lib /nologo /dll /incremental:no /pdb:"$(OUTDIR)\OBExtra.pdb" /machine:I386 /out:"OBExtra.obf" /implib:"$(OUTDIR)\OBExtra.lib" /libpath:"..\OBDLL\Release" /libpath:"..\OBConv\Release" /libpath:".." 
LINK32_OBJS= \
	"$(INTDIR)\fastsearchformat.obj" \
	"$(INTDIR)\finger2.obj" \
	"$(INTDIR)\finger3.obj" \
	"$(INTDIR)\fingerprintformat.obj" \
	"$(INTDIR)\inchiformat.obj" \
	"$(INTDIR)\mpdformat.obj"

".\OBExtra.obf" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "OBExtraFormats - Win32 Debug"

OUTDIR=.\Debug
INTDIR=.\Debug
# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

ALL : "$(OUTDIR)\OBExtraD.obf" "$(OUTDIR)\OBExtraFormats.bsc"


CLEAN :
	-@erase "$(INTDIR)\fastsearchformat.obj"
	-@erase "$(INTDIR)\fastsearchformat.sbr"
	-@erase "$(INTDIR)\finger2.obj"
	-@erase "$(INTDIR)\finger2.sbr"
	-@erase "$(INTDIR)\finger3.obj"
	-@erase "$(INTDIR)\finger3.sbr"
	-@erase "$(INTDIR)\fingerprintformat.obj"
	-@erase "$(INTDIR)\fingerprintformat.sbr"
	-@erase "$(INTDIR)\inchiformat.obj"
	-@erase "$(INTDIR)\inchiformat.sbr"
	-@erase "$(INTDIR)\mpdformat.obj"
	-@erase "$(INTDIR)\mpdformat.sbr"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vc60.pdb"
	-@erase "$(OUTDIR)\OBExtraD.exp"
	-@erase "$(OUTDIR)\OBExtraD.ilk"
	-@erase "$(OUTDIR)\OBExtraD.obf"
	-@erase "$(OUTDIR)\OBExtraD.pdb"
	-@erase "$(OUTDIR)\OBExtraFormats.bsc"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP_PROJ=/nologo /MDd /W3 /Gm /GR /GX /ZI /Od /I ".." /I "..\..\src" /I "..\..\data" /D "_DEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "INCHI_LINK_AS_DLL" /D "USING_OBDLL" /D "HAVE_CONFIG_H" /D "USING_DYNAMIC_LIBS" /FR"$(INTDIR)\\" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 
MTL_PROJ=/nologo /D "_DEBUG" /mktyplib203 /win32 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\OBExtraFormats.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\fastsearchformat.sbr" \
	"$(INTDIR)\finger2.sbr" \
	"$(INTDIR)\finger3.sbr" \
	"$(INTDIR)\fingerprintformat.sbr" \
	"$(INTDIR)\inchiformat.sbr" \
	"$(INTDIR)\mpdformat.sbr"

"$(OUTDIR)\OBExtraFormats.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib obconv.lib obdll.lib libinchi.lib /nologo /dll /incremental:yes /pdb:"$(OUTDIR)\OBExtraD.pdb" /debug /machine:I386 /out:"$(OUTDIR)\OBExtraD.obf" /implib:"$(OUTDIR)\OBExtraD.lib" /pdbtype:sept /libpath:"..\OBDLL\Debug" /libpath:"..\OBConv\Debug" /libpath:".." 
LINK32_OBJS= \
	"$(INTDIR)\fastsearchformat.obj" \
	"$(INTDIR)\finger2.obj" \
	"$(INTDIR)\finger3.obj" \
	"$(INTDIR)\fingerprintformat.obj" \
	"$(INTDIR)\inchiformat.obj" \
	"$(INTDIR)\mpdformat.obj"

"$(OUTDIR)\OBExtraD.obf" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

SOURCE="$(InputPath)"
PostBuild_Desc=Copy debug versions of obconv.dll, obdll.dll
DS_POSTBUILD_DEP=$(INTDIR)\postbld.dep

ALL : $(DS_POSTBUILD_DEP)

# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

$(DS_POSTBUILD_DEP) : "$(OUTDIR)\OBExtraD.obf" "$(OUTDIR)\OBExtraFormats.bsc"
   Copy  ..\obconv\debug\obconv.dll  .\debug  /Y
	Copy  ..\obdll\debug\obdll.dll  .\debug  /Y
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
!IF EXISTS("OBExtraFormats.dep")
!INCLUDE "OBExtraFormats.dep"
!ELSE 
!MESSAGE Warning: cannot find "OBExtraFormats.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "OBExtraFormats - Win32 Release" || "$(CFG)" == "OBExtraFormats - Win32 Debug"
SOURCE=..\..\src\formats\fastsearchformat.cpp

!IF  "$(CFG)" == "OBExtraFormats - Win32 Release"


"$(INTDIR)\fastsearchformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBExtraFormats - Win32 Debug"


"$(INTDIR)\fastsearchformat.obj"	"$(INTDIR)\fastsearchformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\fingerprints\finger2.cpp

!IF  "$(CFG)" == "OBExtraFormats - Win32 Release"


"$(INTDIR)\finger2.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBExtraFormats - Win32 Debug"


"$(INTDIR)\finger2.obj"	"$(INTDIR)\finger2.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\fingerprints\finger3.cpp

!IF  "$(CFG)" == "OBExtraFormats - Win32 Release"


"$(INTDIR)\finger3.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBExtraFormats - Win32 Debug"


"$(INTDIR)\finger3.obj"	"$(INTDIR)\finger3.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\fingerprintformat.cpp

!IF  "$(CFG)" == "OBExtraFormats - Win32 Release"


"$(INTDIR)\fingerprintformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBExtraFormats - Win32 Debug"


"$(INTDIR)\fingerprintformat.obj"	"$(INTDIR)\fingerprintformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\inchiformat.cpp

!IF  "$(CFG)" == "OBExtraFormats - Win32 Release"


"$(INTDIR)\inchiformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBExtraFormats - Win32 Debug"


"$(INTDIR)\inchiformat.obj"	"$(INTDIR)\inchiformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\mpdformat.cpp

!IF  "$(CFG)" == "OBExtraFormats - Win32 Release"


"$(INTDIR)\mpdformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBExtraFormats - Win32 Debug"


"$(INTDIR)\mpdformat.obj"	"$(INTDIR)\mpdformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 


!ENDIF 

