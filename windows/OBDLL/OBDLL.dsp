# Microsoft Developer Studio Project File - Name="OBDLL" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Dynamic-Link Library" 0x0102

CFG=OBDLL - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "OBDLL.mak".
!MESSAGE 
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

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "OBDLL - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MT /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "OBDLL_EXPORTS" /Yu"stdafx.h" /FD /c
# ADD CPP /nologo /MD /W3 /GR /GX /I "..\..\src" /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "__KCC" /D "USING_DYNAMIC_LIBS" /D "OBDLL_EXPORTS" /FD /c
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x809 /d "NDEBUG"
# ADD RSC /l 0x809 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /machine:I386
# ADD LINK32 kernel32.lib user32.lib /nologo /dll /machine:I386 /out:"OBDLL.dll"

!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MTd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "OBDLL_EXPORTS" /Yu"stdafx.h" /FD /GZ /c
# ADD CPP /nologo /MDd /Gm /GR /GX /ZI /Od /I "..\math ..\src" /I "..\..\src" /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "__KCC" /D "USING_DYNAMIC_LIBS" /D "OBDLL_EXPORTS" /FR /FD /GZ /c
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x809 /d "_DEBUG"
# ADD RSC /l 0x809 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /map:"OBDLL.map" /debug /machine:I386 /pdbtype:sept
# SUBTRACT LINK32 /pdb:none
# Begin Special Build Tool
SOURCE="$(InputPath)"
PostBuild_Desc=Construct new obdll.def from obdll.map
PostBuild_Cmds=map2def2.exe obdll.map
# End Special Build Tool

!ENDIF 

# Begin Target

# Name "OBDLL - Win32 Release"
# Name "OBDLL - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\src\atom.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\base.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\binary.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\bitvec.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\bond.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\cache.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\chains.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\chiral.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\data.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\generic.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\grid.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\matrix.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\math\matrix3x3.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\mol.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\molchrg.cpp
# End Source File
# Begin Source File

SOURCE=.\obdll.def

!IF  "$(CFG)" == "OBDLL - Win32 Release"

!ELSEIF  "$(CFG)" == "OBDLL - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\src\oberror.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\obutil.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\parsmart.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\patty.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\phmodel.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\rand.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\report.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\residue.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\ring.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\rotor.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\tokenst.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\transform.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\typer.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\math\vector3.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter ""
# End Group
# Begin Source File

SOURCE=..\..\src\aromatic.h
# End Source File
# Begin Source File

SOURCE=..\..\src\atomtyp.h
# End Source File
# Begin Source File

SOURCE=..\..\src\babelconfig.h
# End Source File
# Begin Source File

SOURCE=..\..\src\base.h
# End Source File
# Begin Source File

SOURCE=..\..\src\binary.h
# End Source File
# Begin Source File

SOURCE=..\..\src\bitvec.h
# End Source File
# Begin Source File

SOURCE=..\..\src\bondtyp.h
# End Source File
# Begin Source File

SOURCE=..\..\src\bondtyper.h
# End Source File
# Begin Source File

SOURCE=..\..\src\chains.h
# End Source File
# Begin Source File

SOURCE=..\..\src\chiral.h
# End Source File
# Begin Source File

SOURCE=..\..\src\data.h
# End Source File
# Begin Source File

SOURCE=..\..\src\element.h
# End Source File
# Begin Source File

SOURCE=..\..\src\extable.h
# End Source File
# Begin Source File

SOURCE=..\..\src\generic.h
# End Source File
# Begin Source File

SOURCE=..\..\src\grid.h
# End Source File
# Begin Source File

SOURCE=..\..\src\isotope.h
# End Source File
# Begin Source File

SOURCE=..\..\src\matrix.h
# End Source File
# Begin Source File

SOURCE=..\..\src\math\matrix3x3.h
# End Source File
# Begin Source File

SOURCE=..\..\src\mol.h
# End Source File
# Begin Source File

SOURCE=..\..\src\molchrg.h
# End Source File
# Begin Source File

SOURCE=..\..\src\oberror.h
# End Source File
# Begin Source File

SOURCE=..\..\src\obifstream.h
# End Source File
# Begin Source File

SOURCE=..\..\src\obutil.h
# End Source File
# Begin Source File

SOURCE=..\..\src\parsmart.h
# End Source File
# Begin Source File

SOURCE=..\..\src\patty.h
# End Source File
# Begin Source File

SOURCE=..\..\src\phmodel.h
# End Source File
# Begin Source File

SOURCE=..\..\src\phmodeldata.h
# End Source File
# Begin Source File

SOURCE=..\..\src\resdata.h
# End Source File
# Begin Source File

SOURCE=..\..\src\ring.h
# End Source File
# Begin Source File

SOURCE=..\..\src\rotor.h
# End Source File
# Begin Source File

SOURCE=..\..\src\snprintf.h
# End Source File
# Begin Source File

SOURCE=..\..\src\typer.h
# End Source File
# Begin Source File

SOURCE=..\..\src\types.h
# End Source File
# Begin Source File

SOURCE=..\..\src\math\vector3.h
# End Source File
# End Target
# End Project
