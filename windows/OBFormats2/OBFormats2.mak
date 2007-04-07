# Microsoft Developer Studio Generated NMAKE File, Based on OBFormats2.dsp
!IF "$(CFG)" == ""
CFG=OBFormats2 - Win32 Debug
!MESSAGE No configuration specified. Defaulting to OBFormats2 - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "OBFormats2 - Win32 Release" && "$(CFG)" != "OBFormats2 - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "OBFormats2.mak" CFG="OBFormats2 - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "OBFormats2 - Win32 Release" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE "OBFormats2 - Win32 Debug" (based on "Win32 (x86) Dynamic-Link Library")
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

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release

ALL : ".\OBFormats2.obf"


CLEAN :
	-@erase "$(INTDIR)\alchemyformat.obj"
	-@erase "$(INTDIR)\amberformat.obj"
	-@erase "$(INTDIR)\APIInterface.obj"
	-@erase "$(INTDIR)\balstformat.obj"
	-@erase "$(INTDIR)\bgfformat.obj"
	-@erase "$(INTDIR)\boxformat.obj"
	-@erase "$(INTDIR)\cacaoformat.obj"
	-@erase "$(INTDIR)\cacheformat.obj"
	-@erase "$(INTDIR)\carformat.obj"
	-@erase "$(INTDIR)\cccformat.obj"
	-@erase "$(INTDIR)\chem3dformat.obj"
	-@erase "$(INTDIR)\chemdrawformat.obj"
	-@erase "$(INTDIR)\chemtoolformat.obj"
	-@erase "$(INTDIR)\CRKformat.obj"
	-@erase "$(INTDIR)\CSRformat.obj"
	-@erase "$(INTDIR)\cssrformat.obj"
	-@erase "$(INTDIR)\dmolformat.obj"
	-@erase "$(INTDIR)\featformat.obj"
	-@erase "$(INTDIR)\fhformat.obj"
	-@erase "$(INTDIR)\freefracformat.obj"
	-@erase "$(INTDIR)\gamessformat.obj"
	-@erase "$(INTDIR)\gaussformat.obj"
	-@erase "$(INTDIR)\ghemicalformat.obj"
	-@erase "$(INTDIR)\gromos96format.obj"
	-@erase "$(INTDIR)\hinformat.obj"
	-@erase "$(INTDIR)\jaguarformat.obj"
	-@erase "$(INTDIR)\mdlformat.obj"
	-@erase "$(INTDIR)\mmodformat.obj"
	-@erase "$(INTDIR)\mol2format.obj"
	-@erase "$(INTDIR)\mopacformat.obj"
	-@erase "$(INTDIR)\mpqcformat.obj"
	-@erase "$(INTDIR)\nwchemformat.obj"
	-@erase "$(INTDIR)\pcmodelformat.obj"
	-@erase "$(INTDIR)\pdbformat.obj"
	-@erase "$(INTDIR)\povrayformat.obj"
	-@erase "$(INTDIR)\PQSformat.obj"
	-@erase "$(INTDIR)\qchemformat.obj"
	-@erase "$(INTDIR)\reportformat.obj"
	-@erase "$(INTDIR)\rxnformat.obj"
	-@erase "$(INTDIR)\shelxformat.obj"
	-@erase "$(INTDIR)\smilesformat.obj"
	-@erase "$(INTDIR)\tinkerformat.obj"
	-@erase "$(INTDIR)\turbomoleformat.obj"
	-@erase "$(INTDIR)\unichemformat.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\viewmolformat.obj"
	-@erase "$(INTDIR)\xedformat.obj"
	-@erase "$(INTDIR)\xyzformat.obj"
	-@erase "$(INTDIR)\yasaraformat.obj"
	-@erase "$(INTDIR)\zindoformat.obj"
	-@erase "$(OUTDIR)\OBFormats2.exp"
	-@erase ".\OBFormats2.obf"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP_PROJ=/nologo /MD /W3 /GR /GX /O2 /I "..\..\data" /I ".." /I "..\..\src" /D "NDEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "USING_DYNAMIC_LIBS" /D "USING_OBDLL" /D "HAVE_CONFIG_H" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 
MTL_PROJ=/nologo /D "NDEBUG" /mktyplib203 /win32 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\OBFormats2.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib obdll.lib obconv.lib /nologo /dll /incremental:no /pdb:"$(OUTDIR)\OBFormats2.pdb" /machine:I386 /out:"OBFormats2.obf" /implib:"$(OUTDIR)\OBFormats2.lib" /libpath:"..\OBDLL\Release" /libpath:"..\OBConv\Release" 
LINK32_OBJS= \
	"$(INTDIR)\alchemyformat.obj" \
	"$(INTDIR)\amberformat.obj" \
	"$(INTDIR)\APIInterface.obj" \
	"$(INTDIR)\balstformat.obj" \
	"$(INTDIR)\bgfformat.obj" \
	"$(INTDIR)\boxformat.obj" \
	"$(INTDIR)\cacaoformat.obj" \
	"$(INTDIR)\cacheformat.obj" \
	"$(INTDIR)\carformat.obj" \
	"$(INTDIR)\cccformat.obj" \
	"$(INTDIR)\chem3dformat.obj" \
	"$(INTDIR)\chemdrawformat.obj" \
	"$(INTDIR)\chemtoolformat.obj" \
	"$(INTDIR)\CRKformat.obj" \
	"$(INTDIR)\CSRformat.obj" \
	"$(INTDIR)\cssrformat.obj" \
	"$(INTDIR)\dmolformat.obj" \
	"$(INTDIR)\featformat.obj" \
	"$(INTDIR)\fhformat.obj" \
	"$(INTDIR)\freefracformat.obj" \
	"$(INTDIR)\gamessformat.obj" \
	"$(INTDIR)\gaussformat.obj" \
	"$(INTDIR)\ghemicalformat.obj" \
	"$(INTDIR)\gromos96format.obj" \
	"$(INTDIR)\hinformat.obj" \
	"$(INTDIR)\jaguarformat.obj" \
	"$(INTDIR)\mdlformat.obj" \
	"$(INTDIR)\mmodformat.obj" \
	"$(INTDIR)\mol2format.obj" \
	"$(INTDIR)\mopacformat.obj" \
	"$(INTDIR)\mpqcformat.obj" \
	"$(INTDIR)\nwchemformat.obj" \
	"$(INTDIR)\pdbformat.obj" \
	"$(INTDIR)\povrayformat.obj" \
	"$(INTDIR)\PQSformat.obj" \
	"$(INTDIR)\qchemformat.obj" \
	"$(INTDIR)\reportformat.obj" \
	"$(INTDIR)\rxnformat.obj" \
	"$(INTDIR)\shelxformat.obj" \
	"$(INTDIR)\smilesformat.obj" \
	"$(INTDIR)\tinkerformat.obj" \
	"$(INTDIR)\turbomoleformat.obj" \
	"$(INTDIR)\unichemformat.obj" \
	"$(INTDIR)\viewmolformat.obj" \
	"$(INTDIR)\xedformat.obj" \
	"$(INTDIR)\xyzformat.obj" \
	"$(INTDIR)\yasaraformat.obj" \
	"$(INTDIR)\zindoformat.obj" \
	"$(INTDIR)\pcmodelformat.obj"

".\OBFormats2.obf" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"

OUTDIR=.\Debug
INTDIR=.\Debug
# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

ALL : "$(OUTDIR)\OBFormats2D.obf" "$(OUTDIR)\OBFormats2.bsc"


CLEAN :
	-@erase "$(INTDIR)\alchemyformat.obj"
	-@erase "$(INTDIR)\alchemyformat.sbr"
	-@erase "$(INTDIR)\amberformat.obj"
	-@erase "$(INTDIR)\amberformat.sbr"
	-@erase "$(INTDIR)\APIInterface.obj"
	-@erase "$(INTDIR)\APIInterface.sbr"
	-@erase "$(INTDIR)\balstformat.obj"
	-@erase "$(INTDIR)\balstformat.sbr"
	-@erase "$(INTDIR)\bgfformat.obj"
	-@erase "$(INTDIR)\bgfformat.sbr"
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
	-@erase "$(INTDIR)\chem3dformat.obj"
	-@erase "$(INTDIR)\chem3dformat.sbr"
	-@erase "$(INTDIR)\chemdrawformat.obj"
	-@erase "$(INTDIR)\chemdrawformat.sbr"
	-@erase "$(INTDIR)\chemtoolformat.obj"
	-@erase "$(INTDIR)\chemtoolformat.sbr"
	-@erase "$(INTDIR)\CRKformat.obj"
	-@erase "$(INTDIR)\CRKformat.sbr"
	-@erase "$(INTDIR)\CSRformat.obj"
	-@erase "$(INTDIR)\CSRformat.sbr"
	-@erase "$(INTDIR)\cssrformat.obj"
	-@erase "$(INTDIR)\cssrformat.sbr"
	-@erase "$(INTDIR)\dmolformat.obj"
	-@erase "$(INTDIR)\dmolformat.sbr"
	-@erase "$(INTDIR)\featformat.obj"
	-@erase "$(INTDIR)\featformat.sbr"
	-@erase "$(INTDIR)\fhformat.obj"
	-@erase "$(INTDIR)\fhformat.sbr"
	-@erase "$(INTDIR)\freefracformat.obj"
	-@erase "$(INTDIR)\freefracformat.sbr"
	-@erase "$(INTDIR)\gamessformat.obj"
	-@erase "$(INTDIR)\gamessformat.sbr"
	-@erase "$(INTDIR)\gaussformat.obj"
	-@erase "$(INTDIR)\gaussformat.sbr"
	-@erase "$(INTDIR)\ghemicalformat.obj"
	-@erase "$(INTDIR)\ghemicalformat.sbr"
	-@erase "$(INTDIR)\gromos96format.obj"
	-@erase "$(INTDIR)\gromos96format.sbr"
	-@erase "$(INTDIR)\hinformat.obj"
	-@erase "$(INTDIR)\hinformat.sbr"
	-@erase "$(INTDIR)\jaguarformat.obj"
	-@erase "$(INTDIR)\jaguarformat.sbr"
	-@erase "$(INTDIR)\mdlformat.obj"
	-@erase "$(INTDIR)\mdlformat.sbr"
	-@erase "$(INTDIR)\mmodformat.obj"
	-@erase "$(INTDIR)\mmodformat.sbr"
	-@erase "$(INTDIR)\mol2format.obj"
	-@erase "$(INTDIR)\mol2format.sbr"
	-@erase "$(INTDIR)\mopacformat.obj"
	-@erase "$(INTDIR)\mopacformat.sbr"
	-@erase "$(INTDIR)\mpqcformat.obj"
	-@erase "$(INTDIR)\mpqcformat.sbr"
	-@erase "$(INTDIR)\nwchemformat.obj"
	-@erase "$(INTDIR)\nwchemformat.sbr"
	-@erase "$(INTDIR)\pcmodelformat.obj"
	-@erase "$(INTDIR)\pcmodelformat.sbr"
	-@erase "$(INTDIR)\pdbformat.obj"
	-@erase "$(INTDIR)\pdbformat.sbr"
	-@erase "$(INTDIR)\povrayformat.obj"
	-@erase "$(INTDIR)\povrayformat.sbr"
	-@erase "$(INTDIR)\PQSformat.obj"
	-@erase "$(INTDIR)\PQSformat.sbr"
	-@erase "$(INTDIR)\qchemformat.obj"
	-@erase "$(INTDIR)\qchemformat.sbr"
	-@erase "$(INTDIR)\reportformat.obj"
	-@erase "$(INTDIR)\reportformat.sbr"
	-@erase "$(INTDIR)\rxnformat.obj"
	-@erase "$(INTDIR)\rxnformat.sbr"
	-@erase "$(INTDIR)\shelxformat.obj"
	-@erase "$(INTDIR)\shelxformat.sbr"
	-@erase "$(INTDIR)\smilesformat.obj"
	-@erase "$(INTDIR)\smilesformat.sbr"
	-@erase "$(INTDIR)\tinkerformat.obj"
	-@erase "$(INTDIR)\tinkerformat.sbr"
	-@erase "$(INTDIR)\turbomoleformat.obj"
	-@erase "$(INTDIR)\turbomoleformat.sbr"
	-@erase "$(INTDIR)\unichemformat.obj"
	-@erase "$(INTDIR)\unichemformat.sbr"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vc60.pdb"
	-@erase "$(INTDIR)\viewmolformat.obj"
	-@erase "$(INTDIR)\viewmolformat.sbr"
	-@erase "$(INTDIR)\xedformat.obj"
	-@erase "$(INTDIR)\xedformat.sbr"
	-@erase "$(INTDIR)\xyzformat.obj"
	-@erase "$(INTDIR)\xyzformat.sbr"
	-@erase "$(INTDIR)\yasaraformat.obj"
	-@erase "$(INTDIR)\yasaraformat.sbr"
	-@erase "$(INTDIR)\zindoformat.obj"
	-@erase "$(INTDIR)\zindoformat.sbr"
	-@erase "$(OUTDIR)\OBFormats2.bsc"
	-@erase "$(OUTDIR)\OBFormats2D.exp"
	-@erase "$(OUTDIR)\OBFormats2D.ilk"
	-@erase "$(OUTDIR)\OBFormats2D.obf"
	-@erase "$(OUTDIR)\OBFormats2D.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP_PROJ=/nologo /MDd /W3 /Gm /GR /GX /ZI /Od /I ".." /I "..\..\src" /I "..\..\data" /D "_DEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "USING_DYNAMIC_LIBS" /D "USING_OBDLL" /D "HAVE_CONFIG_H" /FR"$(INTDIR)\\" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 
MTL_PROJ=/nologo /D "_DEBUG" /mktyplib203 /win32 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\OBFormats2.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\alchemyformat.sbr" \
	"$(INTDIR)\amberformat.sbr" \
	"$(INTDIR)\APIInterface.sbr" \
	"$(INTDIR)\balstformat.sbr" \
	"$(INTDIR)\bgfformat.sbr" \
	"$(INTDIR)\boxformat.sbr" \
	"$(INTDIR)\cacaoformat.sbr" \
	"$(INTDIR)\cacheformat.sbr" \
	"$(INTDIR)\carformat.sbr" \
	"$(INTDIR)\cccformat.sbr" \
	"$(INTDIR)\chem3dformat.sbr" \
	"$(INTDIR)\chemdrawformat.sbr" \
	"$(INTDIR)\chemtoolformat.sbr" \
	"$(INTDIR)\CRKformat.sbr" \
	"$(INTDIR)\CSRformat.sbr" \
	"$(INTDIR)\cssrformat.sbr" \
	"$(INTDIR)\dmolformat.sbr" \
	"$(INTDIR)\featformat.sbr" \
	"$(INTDIR)\fhformat.sbr" \
	"$(INTDIR)\freefracformat.sbr" \
	"$(INTDIR)\gamessformat.sbr" \
	"$(INTDIR)\gaussformat.sbr" \
	"$(INTDIR)\ghemicalformat.sbr" \
	"$(INTDIR)\gromos96format.sbr" \
	"$(INTDIR)\hinformat.sbr" \
	"$(INTDIR)\jaguarformat.sbr" \
	"$(INTDIR)\mdlformat.sbr" \
	"$(INTDIR)\mmodformat.sbr" \
	"$(INTDIR)\mol2format.sbr" \
	"$(INTDIR)\mopacformat.sbr" \
	"$(INTDIR)\mpqcformat.sbr" \
	"$(INTDIR)\nwchemformat.sbr" \
	"$(INTDIR)\pdbformat.sbr" \
	"$(INTDIR)\povrayformat.sbr" \
	"$(INTDIR)\PQSformat.sbr" \
	"$(INTDIR)\qchemformat.sbr" \
	"$(INTDIR)\reportformat.sbr" \
	"$(INTDIR)\rxnformat.sbr" \
	"$(INTDIR)\shelxformat.sbr" \
	"$(INTDIR)\smilesformat.sbr" \
	"$(INTDIR)\tinkerformat.sbr" \
	"$(INTDIR)\turbomoleformat.sbr" \
	"$(INTDIR)\unichemformat.sbr" \
	"$(INTDIR)\viewmolformat.sbr" \
	"$(INTDIR)\xedformat.sbr" \
	"$(INTDIR)\xyzformat.sbr" \
	"$(INTDIR)\yasaraformat.sbr" \
	"$(INTDIR)\zindoformat.sbr" \
	"$(INTDIR)\pcmodelformat.sbr"

"$(OUTDIR)\OBFormats2.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib obconv.lib obdll.lib /nologo /dll /incremental:yes /pdb:"$(OUTDIR)\OBFormats2D.pdb" /debug /machine:I386 /out:"$(OUTDIR)\OBFormats2D.obf" /implib:"$(OUTDIR)\OBFormats2D.lib" /pdbtype:sept /libpath:"..\OBDLL\Debug" /libpath:"..\OBConv\Debug" 
LINK32_OBJS= \
	"$(INTDIR)\alchemyformat.obj" \
	"$(INTDIR)\amberformat.obj" \
	"$(INTDIR)\APIInterface.obj" \
	"$(INTDIR)\balstformat.obj" \
	"$(INTDIR)\bgfformat.obj" \
	"$(INTDIR)\boxformat.obj" \
	"$(INTDIR)\cacaoformat.obj" \
	"$(INTDIR)\cacheformat.obj" \
	"$(INTDIR)\carformat.obj" \
	"$(INTDIR)\cccformat.obj" \
	"$(INTDIR)\chem3dformat.obj" \
	"$(INTDIR)\chemdrawformat.obj" \
	"$(INTDIR)\chemtoolformat.obj" \
	"$(INTDIR)\CRKformat.obj" \
	"$(INTDIR)\CSRformat.obj" \
	"$(INTDIR)\cssrformat.obj" \
	"$(INTDIR)\dmolformat.obj" \
	"$(INTDIR)\featformat.obj" \
	"$(INTDIR)\fhformat.obj" \
	"$(INTDIR)\freefracformat.obj" \
	"$(INTDIR)\gamessformat.obj" \
	"$(INTDIR)\gaussformat.obj" \
	"$(INTDIR)\ghemicalformat.obj" \
	"$(INTDIR)\gromos96format.obj" \
	"$(INTDIR)\hinformat.obj" \
	"$(INTDIR)\jaguarformat.obj" \
	"$(INTDIR)\mdlformat.obj" \
	"$(INTDIR)\mmodformat.obj" \
	"$(INTDIR)\mol2format.obj" \
	"$(INTDIR)\mopacformat.obj" \
	"$(INTDIR)\mpqcformat.obj" \
	"$(INTDIR)\nwchemformat.obj" \
	"$(INTDIR)\pdbformat.obj" \
	"$(INTDIR)\povrayformat.obj" \
	"$(INTDIR)\PQSformat.obj" \
	"$(INTDIR)\qchemformat.obj" \
	"$(INTDIR)\reportformat.obj" \
	"$(INTDIR)\rxnformat.obj" \
	"$(INTDIR)\shelxformat.obj" \
	"$(INTDIR)\smilesformat.obj" \
	"$(INTDIR)\tinkerformat.obj" \
	"$(INTDIR)\turbomoleformat.obj" \
	"$(INTDIR)\unichemformat.obj" \
	"$(INTDIR)\viewmolformat.obj" \
	"$(INTDIR)\xedformat.obj" \
	"$(INTDIR)\xyzformat.obj" \
	"$(INTDIR)\yasaraformat.obj" \
	"$(INTDIR)\zindoformat.obj" \
	"$(INTDIR)\pcmodelformat.obj"

"$(OUTDIR)\OBFormats2D.obf" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
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

$(DS_POSTBUILD_DEP) : "$(OUTDIR)\OBFormats2D.obf" "$(OUTDIR)\OBFormats2.bsc"
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
!IF EXISTS("OBFormats2.dep")
!INCLUDE "OBFormats2.dep"
!ELSE 
!MESSAGE Warning: cannot find "OBFormats2.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "OBFormats2 - Win32 Release" || "$(CFG)" == "OBFormats2 - Win32 Debug"
SOURCE=..\..\src\formats\alchemyformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\alchemyformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\alchemyformat.obj"	"$(INTDIR)\alchemyformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\amberformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\amberformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\amberformat.obj"	"$(INTDIR)\amberformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\APIInterface.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\APIInterface.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\APIInterface.obj"	"$(INTDIR)\APIInterface.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\balstformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\balstformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\balstformat.obj"	"$(INTDIR)\balstformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\bgfformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\bgfformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\bgfformat.obj"	"$(INTDIR)\bgfformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\boxformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\boxformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\boxformat.obj"	"$(INTDIR)\boxformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\cacaoformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\cacaoformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\cacaoformat.obj"	"$(INTDIR)\cacaoformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\cacheformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\cacheformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\cacheformat.obj"	"$(INTDIR)\cacheformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\carformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\carformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\carformat.obj"	"$(INTDIR)\carformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\cccformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\cccformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\cccformat.obj"	"$(INTDIR)\cccformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\chem3dformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\chem3dformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\chem3dformat.obj"	"$(INTDIR)\chem3dformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\chemdrawformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\chemdrawformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\chemdrawformat.obj"	"$(INTDIR)\chemdrawformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\chemtoolformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\chemtoolformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\chemtoolformat.obj"	"$(INTDIR)\chemtoolformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\CRKformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\CRKformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\CRKformat.obj"	"$(INTDIR)\CRKformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\CSRformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\CSRformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\CSRformat.obj"	"$(INTDIR)\CSRformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\cssrformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\cssrformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\cssrformat.obj"	"$(INTDIR)\cssrformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\dmolformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\dmolformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\dmolformat.obj"	"$(INTDIR)\dmolformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\featformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\featformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\featformat.obj"	"$(INTDIR)\featformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\fhformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\fhformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\fhformat.obj"	"$(INTDIR)\fhformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\freefracformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\freefracformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\freefracformat.obj"	"$(INTDIR)\freefracformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\gamessformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\gamessformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\gamessformat.obj"	"$(INTDIR)\gamessformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\gaussformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\gaussformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\gaussformat.obj"	"$(INTDIR)\gaussformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\ghemicalformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\ghemicalformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\ghemicalformat.obj"	"$(INTDIR)\ghemicalformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\gromos96format.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\gromos96format.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\gromos96format.obj"	"$(INTDIR)\gromos96format.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\hinformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\hinformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\hinformat.obj"	"$(INTDIR)\hinformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\jaguarformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\jaguarformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\jaguarformat.obj"	"$(INTDIR)\jaguarformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\mdlformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\mdlformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\mdlformat.obj"	"$(INTDIR)\mdlformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\mmodformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\mmodformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\mmodformat.obj"	"$(INTDIR)\mmodformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\mol2format.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\mol2format.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\mol2format.obj"	"$(INTDIR)\mol2format.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\mopacformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\mopacformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\mopacformat.obj"	"$(INTDIR)\mopacformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\mpqcformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\mpqcformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\mpqcformat.obj"	"$(INTDIR)\mpqcformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\nwchemformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\nwchemformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\nwchemformat.obj"	"$(INTDIR)\nwchemformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\pcmodelformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\pcmodelformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\pcmodelformat.obj"	"$(INTDIR)\pcmodelformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\pdbformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\pdbformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\pdbformat.obj"	"$(INTDIR)\pdbformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\povrayformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\povrayformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\povrayformat.obj"	"$(INTDIR)\povrayformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\PQSformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\PQSformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\PQSformat.obj"	"$(INTDIR)\PQSformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\qchemformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\qchemformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\qchemformat.obj"	"$(INTDIR)\qchemformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\reportformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\reportformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\reportformat.obj"	"$(INTDIR)\reportformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\rxnformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\rxnformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\rxnformat.obj"	"$(INTDIR)\rxnformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\shelxformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\shelxformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\shelxformat.obj"	"$(INTDIR)\shelxformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\smilesformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\smilesformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\smilesformat.obj"	"$(INTDIR)\smilesformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\tinkerformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\tinkerformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\tinkerformat.obj"	"$(INTDIR)\tinkerformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\turbomoleformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\turbomoleformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\turbomoleformat.obj"	"$(INTDIR)\turbomoleformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\unichemformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\unichemformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\unichemformat.obj"	"$(INTDIR)\unichemformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\viewmolformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\viewmolformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\viewmolformat.obj"	"$(INTDIR)\viewmolformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\xedformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\xedformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\xedformat.obj"	"$(INTDIR)\xedformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\xyzformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\xyzformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\xyzformat.obj"	"$(INTDIR)\xyzformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\yasaraformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\yasaraformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\yasaraformat.obj"	"$(INTDIR)\yasaraformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\src\formats\zindoformat.cpp

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"


"$(INTDIR)\zindoformat.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"


"$(INTDIR)\zindoformat.obj"	"$(INTDIR)\zindoformat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 


!ENDIF 

