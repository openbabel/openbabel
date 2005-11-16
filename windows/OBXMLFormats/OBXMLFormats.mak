# Microsoft Developer Studio Generated NMAKE File, Based on OBXMLFormats.dsp
!IF "$(CFG)" == ""
CFG=OBXMLFormats - Win32 Debug
!MESSAGE No configuration specified. Defaulting to OBXMLFormats - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "OBXMLFormats - Win32 Release" && "$(CFG)" != "OBXMLFormats - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "OBXMLFormats.mak" CFG="OBXMLFormats - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "OBXMLFormats - Win32 Release" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE "OBXMLFormats - Win32 Debug" (based on "Win32 (x86) Dynamic-Link Library")
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

!IF  "$(CFG)" == "OBXMLFormats - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release

ALL : ".\OBXML.obf"


CLEAN :
	-@erase "$(INTDIR)\cmlreactlformat.obj"
	-@erase "$(INTDIR)\pubchem.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\xcmlformat.obj"
	-@erase "$(INTDIR)\xml.obj"
	-@erase "$(INTDIR)\xmlformat.obj"
	-@erase "$(OUTDIR)\OBXML.exp"
	-@erase ".\OBXML.obf"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP_PROJ=/nologo /MD /W3 /GR /GX /O2 /I "..\..\data" /I ".." /I "..\..\src" /I "..\..\src\formats" /I "..\..\src\formats\xml" /D "NDEBUG" /D "USING_DYNAMIC_LIBS" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "INCHI_LINK_AS_DLL" /D "USING_OBDLL" /D "HAVE_CONFIG_H" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 
MTL_PROJ=/nologo /D "NDEBUG" /mktyplib203 /win32 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\OBXMLFormats.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib obconv.lib obdll.lib libinchi.lib /nologo /dll /incremental:no /pdb:"$(OUTDIR)\OBXML.pdb" /machine:I386 /out:"OBXML.obf" /implib:"$(OUTDIR)\OBXML.lib" /libpath:"..\OBDLL\Release" /libpath:"..\OBConv\Release" /libpath:".." 
LINK32_OBJS= \
	"$(INTDIR)\cmlreactlformat.obj" \
	"$(INTDIR)\pubchem.obj" \
	"$(INTDIR)\xcmlformat.obj" \
	"$(INTDIR)\xml.obj" \
	"$(INTDIR)\xmlformat.obj" \
	"..\libxml2.lib"

".\OBXML.obf" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "OBXMLFormats - Win32 Debug"

OUTDIR=.\Debug
INTDIR=.\Debug
# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

ALL : "$(OUTDIR)\OBXMLD.obf" "$(OUTDIR)\OBXMLFormats.bsc"


CLEAN :
	-@erase "$(INTDIR)\cmlreactlformat.obj"
	-@erase "$(INTDIR)\cmlreactlformat.sbr"
	-@erase "$(INTDIR)\pubchem.obj"
	-@erase "$(INTDIR)\pubchem.sbr"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vc60.pdb"
	-@erase "$(INTDIR)\xcmlformat.obj"
	-@erase "$(INTDIR)\xcmlformat.sbr"
	-@erase "$(INTDIR)\xml.obj"
	-@erase "$(INTDIR)\xml.sbr"
	-@erase "$(INTDIR)\xmlformat.obj"
	-@erase "$(INTDIR)\xmlformat.sbr"
	-@erase "$(OUTDIR)\OBXMLD.exp"
	-@erase "$(OUTDIR)\OBXMLD.ilk"
	-@erase "$(OUTDIR)\OBXMLD.obf"
	-@erase "$(OUTDIR)\OBXMLD.pdb"
	-@erase "$(OUTDIR)\OBXMLFormats.bsc"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP_PROJ=/nologo /MDd /W3 /Gm /GR /GX /ZI /Od /I ".." /I "..\..\src" /I "..\..\data" /I "..\..\src\formats" /I "..\..\src\formats\xml" /D "_DEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "INCHI_LINK_AS_DLL" /D "USING_OBDLL" /D "HAVE_CONFIG_H" /D "USING_DYNAMIC_LIBS" /FR"$(INTDIR)\\" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 
MTL_PROJ=/nologo /D "_DEBUG" /mktyplib203 /win32 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\OBXMLFormats.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\cmlreactlformat.sbr" \
	"$(INTDIR)\pubchem.sbr" \
	"$(INTDIR)\xcmlformat.sbr" \
	"$(INTDIR)\xml.sbr" \
	"$(INTDIR)\xmlformat.sbr"

"$(OUTDIR)\OBXMLFormats.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib obconv.lib obdll.lib libinchi.lib /nologo /dll /incremental:yes /pdb:"$(OUTDIR)\OBXMLD.pdb" /debug /machine:I386 /out:"$(OUTDIR)\OBXMLD.obf" /implib:"$(OUTDIR)\OBXMLD.lib" /pdbtype:sept /libpath:"..\OBDLL\Debug" /libpath:"..\OBConv\Debug" /libpath:".." 
LINK32_OBJS= \
	"$(INTDIR)\cmlreactlformat.obj" \
	"$(INTDIR)\pubchem.obj" \
	"$(INTDIR)\xcmlformat.obj" \
	"$(INTDIR)\xml.obj" \
	"$(INTDIR)\xmlformat.obj" \
	"..\libxml2.lib"

"$(OUTDIR)\OBXMLD.obf" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
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

$(DS_POSTBUILD_DEP) : "$(OUTDIR)\OBXMLD.obf" "$(OUTDIR)\OBXMLFormats.bsc"
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
!IF EXISTS("OBXMLFormats.dep")
!INCLUDE "OBXMLFormats.dep"
!ELSE 
!MESSAGE Warning: cannot find "OBXMLFormats.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "OBXMLFormats - Win32 Release" || "$(CFG)" == "OBXMLFormats - Win32 Debug"
SOURCE=..\..\src\formats\xml\cmlreactlformat.cpp

!IF  "$(CFG)" == "OBXMLFormats - Win32 Release"


"$(INTDIR)\cmlreactlformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBXMLFormats - Win32 Debug"


"$(INTDIR)\cmlreactlformat.obj"	"$(INTDIR)\cmlreactlformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\xml\pubchem.cpp

!IF  "$(CFG)" == "OBXMLFormats - Win32 Release"


"$(INTDIR)\pubchem.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBXMLFormats - Win32 Debug"


"$(INTDIR)\pubchem.obj"	"$(INTDIR)\pubchem.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\xml\xcmlformat.cpp

!IF  "$(CFG)" == "OBXMLFormats - Win32 Release"


"$(INTDIR)\xcmlformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBXMLFormats - Win32 Debug"


"$(INTDIR)\xcmlformat.obj"	"$(INTDIR)\xcmlformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\xml\xml.cpp

!IF  "$(CFG)" == "OBXMLFormats - Win32 Release"


"$(INTDIR)\xml.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBXMLFormats - Win32 Debug"


"$(INTDIR)\xml.obj"	"$(INTDIR)\xml.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\xml\xmlformat.cpp

!IF  "$(CFG)" == "OBXMLFormats - Win32 Release"


"$(INTDIR)\xmlformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBXMLFormats - Win32 Debug"


"$(INTDIR)\xmlformat.obj"	"$(INTDIR)\xmlformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 


!ENDIF 

