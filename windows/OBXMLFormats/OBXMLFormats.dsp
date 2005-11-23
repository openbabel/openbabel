# Microsoft Developer Studio Project File - Name="OBXMLFormats" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Dynamic-Link Library" 0x0102

CFG=OBXMLFormats - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "OBXMLFormats.mak".
!MESSAGE 
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

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "OBXMLFormats - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 1
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MT /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "OBXMLFormats_EXPORTS" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GR /GX /O2 /I "..\..\data" /I ".." /I "..\..\src" /I "..\..\src\formats" /I "..\..\src\formats\xml" /D "NDEBUG" /D "USING_DYNAMIC_LIBS" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "INCHI_LINK_AS_DLL" /D "USING_OBDLL" /D "HAVE_CONFIG_H" /FD /c
# SUBTRACT CPP /YX
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x809 /d "NDEBUG"
# ADD RSC /l 0x809 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /machine:I386
# ADD LINK32 gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib obconv.lib obdll.lib libinchi.lib /nologo /dll /machine:I386 /out:"OBXML.obf" /libpath:"..\OBDLL\Release" /libpath:"..\OBConv\Release" /libpath:".."

!ELSEIF  "$(CFG)" == "OBXMLFormats - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 1
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MTd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "OBXMLFormats_EXPORTS" /YX /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /Gm /GR /GX /ZI /Od /I ".." /I "..\..\src" /I "..\..\data" /I "..\..\src\formats" /I "..\..\src\formats\xml" /D "_DEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "INCHI_LINK_AS_DLL" /D "USING_OBDLL" /D "HAVE_CONFIG_H" /D "USING_DYNAMIC_LIBS" /FR /FD /GZ /c
# SUBTRACT CPP /YX
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x809 /d "_DEBUG"
# ADD RSC /l 0x809 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib obconv.lib obdll.lib libinchi.lib /nologo /dll /debug /machine:I386 /out:"Debug/OBXMLD.obf" /pdbtype:sept /libpath:"..\OBDLL\Debug" /libpath:"..\OBConv\Debug" /libpath:".."
# Begin Special Build Tool
SOURCE="$(InputPath)"
PostBuild_Desc=Copy debug versions of obconv.dll, obdll.dll
PostBuild_Cmds=Copy  ..\obconv\debug\obconv.dll  .\debug  /Y	Copy  ..\obdll\debug\obdll.dll  .\debug  /Y
# End Special Build Tool

!ENDIF 

# Begin Target

# Name "OBXMLFormats - Win32 Release"
# Name "OBXMLFormats - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\src\formats\xml\cmlreactlformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\xml\pubchem.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\xml\xcmlformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\xml\xml.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\xml\xmlformat.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\src\formats\xml\iconv.h
# End Source File
# Begin Source File

SOURCE=..\..\src\mol.h
# End Source File
# Begin Source File

SOURCE=..\..\src\obconversion.h
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\obmolecformat.h
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\xml\xml.h
# End Source File
# End Group
# Begin Source File

SOURCE=..\libxml2.lib
# End Source File
# End Target
# End Project
