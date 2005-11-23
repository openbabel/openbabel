# Microsoft Developer Studio Generated NMAKE File, Based on obgrep.dsp
!IF "$(CFG)" == ""
CFG=obgrep - Win32 Debug
!MESSAGE No configuration specified. Defaulting to obgrep - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "obgrep - Win32 Release" && "$(CFG)" != "obgrep - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "obgrep.mak" CFG="obgrep - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "obgrep - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "obgrep - Win32 Debug" (based on "Win32 (x86) Console Application")
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

!IF  "$(CFG)" == "obgrep - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release

ALL : ".\obgrep.exe"


CLEAN :
	-@erase "$(INTDIR)\obgrep.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\XGetopt.obj"
	-@erase ".\obgrep.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP_PROJ=/nologo /MD /W3 /GR /GX /O2 /I ".." /I "..\..\src" /D "NDEBUG" /D "WIN32" /D "_CONSOLE" /D "_MBCS" /D "HAVE_CONFIG_H" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\obgrep.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib OBDLL.lib OBConv.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /incremental:no /pdb:"$(OUTDIR)\obgrep.pdb" /machine:I386 /out:"obgrep.exe" /libpath:"..\OBConv\Release" /libpath:"..\OBDLL\Release" 
LINK32_OBJS= \
	"$(INTDIR)\obgrep.obj" \
	"$(INTDIR)\XGetopt.obj"

".\obgrep.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "obgrep - Win32 Debug"

OUTDIR=.\Debug
INTDIR=.\Debug
# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

ALL : "$(OUTDIR)\obgrep.exe"


CLEAN :
	-@erase "$(INTDIR)\obgrep.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vc60.pdb"
	-@erase "$(INTDIR)\XGetopt.obj"
	-@erase "$(OUTDIR)\obgrep.exe"
	-@erase "$(OUTDIR)\obgrep.ilk"
	-@erase "$(OUTDIR)\obgrep.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP_PROJ=/nologo /MDd /W3 /Gm /GR /GX /ZI /Od /I ".." /I "..\..\src" /D "_DEBUG" /D "WIN32" /D "_CONSOLE" /D "_MBCS" /D "HAVE_CONFIG_H" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\obgrep.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib OBDLL.lib OBConv.lib /nologo /subsystem:console /incremental:yes /pdb:"$(OUTDIR)\obgrep.pdb" /debug /machine:I386 /out:"$(OUTDIR)\obgrep.exe" /pdbtype:sept /libpath:"..\OBConv\debug" /libpath:"..\OBDLL\debug" 
LINK32_OBJS= \
	"$(INTDIR)\obgrep.obj" \
	"$(INTDIR)\XGetopt.obj"

"$(OUTDIR)\obgrep.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
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
!IF EXISTS("obgrep.dep")
!INCLUDE "obgrep.dep"
!ELSE 
!MESSAGE Warning: cannot find "obgrep.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "obgrep - Win32 Release" || "$(CFG)" == "obgrep - Win32 Debug"
SOURCE=..\..\tools\obgrep.cpp

"$(INTDIR)\obgrep.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\XGetopt.cpp

"$(INTDIR)\XGetopt.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)



!ENDIF 

