# Microsoft Developer Studio Project File - Name="OBFormats2" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Dynamic-Link Library" 0x0102

CFG=OBFormats2 - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "OBFormats2.mak".
!MESSAGE 
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

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "OBFormats2 - Win32 Release"

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
# ADD BASE CPP /nologo /MT /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "OBFormats2_EXPORTS" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GR /GX /O2 /I "..\..\data" /I ".." /I "..\..\src" /D "NDEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "USING_DYNAMIC_LIBS" /D "USING_OBDLL" /D "HAVE_CONFIG_H" /FD /c
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
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib obdll.lib obconv.lib /nologo /dll /machine:I386 /out:"OBFormats2.obf" /libpath:"..\OBDLL\Release" /libpath:"..\OBConv\Release"

!ELSEIF  "$(CFG)" == "OBFormats2 - Win32 Debug"

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
# ADD BASE CPP /nologo /MTd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "OBFormats2_EXPORTS" /YX /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /Gm /GR /GX /ZI /Od /I ".." /I "..\..\src" /I "..\..\data" /D "_DEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "USING_DYNAMIC_LIBS" /D "USING_OBDLL" /D "HAVE_CONFIG_H" /FR /FD /GZ /c
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
# ADD LINK32 kernel32.lib user32.lib obconv.lib obdll.lib /nologo /dll /debug /machine:I386 /out:"Debug/OBFormats2D.obf" /pdbtype:sept /libpath:"..\OBDLL\Debug" /libpath:"..\OBConv\Debug"
# Begin Special Build Tool
SOURCE="$(InputPath)"
PostBuild_Desc=Copy debug versions of obconv.dll, obdll.dll
PostBuild_Cmds=Copy  ..\obconv\debug\obconv.dll  .\debug  /Y	Copy  ..\obdll\debug\obdll.dll  .\debug  /Y
# End Special Build Tool

!ENDIF 

# Begin Target

# Name "OBFormats2 - Win32 Release"
# Name "OBFormats2 - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\src\formats\alchemyformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\amberformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\APIInterface.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\balstformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\bgfformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\boxformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\cacaoformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\cacheformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\carformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\cccformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\chem3dformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\chemdrawformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\chemtoolformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\CRKformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\CSRformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\cssrformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\dmolformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\featformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\fhformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\freefracformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\gamessformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\gaussformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\ghemicalformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\gromos96format.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\hinformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\jaguarformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\mdlformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\mmodformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\mol2format.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\mopacformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\mpqcformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\nwchemformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\pcmodelformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\pdbformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\povrayformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\PQSformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\qchemformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\reportformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\rxnformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\shelxformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\smilesformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\tinkerformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\turbomoleformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\unichemformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\viewmolformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\xedformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\xyzformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\yasaraformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\zindoformat.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\src\mol.h
# End Source File
# Begin Source File

SOURCE=..\..\src\obconversion.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# Begin Source File

SOURCE=.\Stereochemistry.txt
# End Source File
# End Target
# End Project
