# Microsoft Developer Studio Generated NMAKE File, Based on OBConv.dsp
!IF "$(CFG)" == ""
CFG=OBConv - Win32 Debug
!MESSAGE No configuration specified. Defaulting to OBConv - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "OBConv - Win32 Release" && "$(CFG)" != "OBConv - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "OBConv.mak" CFG="OBConv - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "OBConv - Win32 Release" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE "OBConv - Win32 Debug" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

!IF  "$(CFG)" == "OBConv - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release

ALL : ".\OBConv.dll"


CLEAN :
	-@erase "$(INTDIR)\dlhandler_win32..obj"
	-@erase "$(INTDIR)\OBConvDLLinit.obj"
	-@erase "$(INTDIR)\obconversion.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\OBConv.exp"
	-@erase "$(OUTDIR)\OBConv.lib"
	-@erase ".\OBConv.dll"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MD /W3 /GR /GX /O2 /I ".." /D "NDEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "OBCONV_EXPORTS" /D "USING_DYNAMIC_LIBS" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

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
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\OBConv.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib /nologo /dll /incremental:no /pdb:"$(OUTDIR)\OBConv.pdb" /machine:I386 /out:"OBConv.dll" /implib:"$(OUTDIR)\OBConv.lib" 
LINK32_OBJS= \
	"$(INTDIR)\dlhandler_win32..obj" \
	"$(INTDIR)\OBConvDLLinit.obj" \
	"$(INTDIR)\obconversion.obj"

".\OBConv.dll" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "OBConv - Win32 Debug"

OUTDIR=.\Debug
INTDIR=.\Debug
# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

ALL : "$(OUTDIR)\OBConv.dll" "$(OUTDIR)\OBConv.bsc"


CLEAN :
	-@erase "$(INTDIR)\dlhandler_win32..obj"
	-@erase "$(INTDIR)\dlhandler_win32..sbr"
	-@erase "$(INTDIR)\obconversion.obj"
	-@erase "$(INTDIR)\obconversion.sbr"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vc60.pdb"
	-@erase "$(OUTDIR)\OBConv.bsc"
	-@erase "$(OUTDIR)\OBConv.dll"
	-@erase "$(OUTDIR)\OBConv.exp"
	-@erase "$(OUTDIR)\OBConv.ilk"
	-@erase "$(OUTDIR)\OBConv.lib"
	-@erase "$(OUTDIR)\OBConv.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MDd /W3 /Gm /GR /GX /ZI /Od /I ".." /D "_DEBUG" /D "USING_DYNAMIC_LIBS" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "OBCONV_EXPORTS" /FR"$(INTDIR)\\" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

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
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\OBConv.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\dlhandler_win32..sbr" \
	"$(INTDIR)\obconversion.sbr"

"$(OUTDIR)\OBConv.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /incremental:yes /pdb:"$(OUTDIR)\OBConv.pdb" /debug /machine:I386 /out:"$(OUTDIR)\OBConv.dll" /implib:"$(OUTDIR)\OBConv.lib" /pdbtype:sept 
LINK32_OBJS= \
	"$(INTDIR)\dlhandler_win32..obj" \
	"$(INTDIR)\obconversion.obj"

"$(OUTDIR)\OBConv.dll" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("OBConv.dep")
!INCLUDE "OBConv.dep"
!ELSE 
!MESSAGE Warning: cannot find "OBConv.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "OBConv - Win32 Release" || "$(CFG)" == "OBConv - Win32 Debug"
SOURCE=..\..\src\dlhandler_win32..cpp

!IF  "$(CFG)" == "OBConv - Win32 Release"


"$(INTDIR)\dlhandler_win32..obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBConv - Win32 Debug"


"$(INTDIR)\dlhandler_win32..obj"	"$(INTDIR)\dlhandler_win32..sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=.\OBConvDLLinit.cpp

!IF  "$(CFG)" == "OBConv - Win32 Release"


"$(INTDIR)\OBConvDLLinit.obj" : $(SOURCE) "$(INTDIR)"


!ELSEIF  "$(CFG)" == "OBConv - Win32 Debug"

!ENDIF 

SOURCE=..\..\src\obconversion.cpp

!IF  "$(CFG)" == "OBConv - Win32 Release"


"$(INTDIR)\obconversion.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "OBConv - Win32 Debug"


"$(INTDIR)\obconversion.obj"	"$(INTDIR)\obconversion.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 


!ENDIF 

