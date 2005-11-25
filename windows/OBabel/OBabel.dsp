# Microsoft Developer Studio Project File - Name="OBabel" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=OBabel - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "OBabel.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "OBabel.mak" CFG="OBabel - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "OBabel - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "OBabel - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "OBabel - Win32 Release"

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
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /MT /W3 /GR /GX /O1 /I "..\..\src" /I ".." /I "../../data" /I "..\..\src\formats" /I "..\..\src\formats\xml" /D "NDEBUG" /D "WIN32" /D "_CONSOLE" /D "_MBCS" /D "INCHI_LINK_AS_DLL" /D "HAVE_CONFIG_H" /YX /FD /c
# SUBTRACT CPP /Fr
# ADD BASE RSC /l 0x809 /d "NDEBUG"
# ADD RSC /l 0x809 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 libinchi.lib /nologo /subsystem:console /machine:I386 /out:"babel.exe" /libpath:".."
# SUBTRACT LINK32 /debug

!ELSEIF  "$(CFG)" == "OBabel - Win32 Debug"

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
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /MTd /W3 /Gm /GR /GX /ZI /Od /I "..\..\src" /I ".." /I "../../data" /I "..\..\src\formats" /I "..\..\src\formats\xml" /D "_DEBUG" /D "WIN32" /D "_CONSOLE" /D "_MBCS" /D "INCHI_LINK_AS_DLL" /D "HAVE_CONFIG_H" /FR /FD /GZ /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x809 /d "_DEBUG"
# ADD RSC /l 0x809 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib libinchi.lib /nologo /subsystem:console /debug /machine:I386 /out:"Debug/babel.exe" /pdbtype:sept /libpath:".."
# SUBTRACT LINK32 /profile

!ENDIF 

# Begin Target

# Name "OBabel - Win32 Release"
# Name "OBabel - Win32 Debug"
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

SOURCE=..\..\src\atom.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\balstformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\base.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\bgfformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\bitvec.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\bond.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\bondtyper.cpp
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

SOURCE=..\..\src\chains.cpp
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

SOURCE=..\..\src\chiral.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\xml\cmlreactlformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\copyformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\crkformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\CSRformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\cssrformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\data.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\dlhandler_win32.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\dmolformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\fastsearchformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\featformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\fhformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\fingerprints\finger2.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\fingerprints\finger3.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\fingerprint.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\fingerprintformat.cpp
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

SOURCE=..\..\src\generic.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\ghemicalformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\grid.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\gromos96format.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\hinformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\inchiformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\jaguarformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\kekulize.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\main.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\matrix.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\math\matrix3x3.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\mdlformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\mmodformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\mol.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\mol2format.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\molchrg.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\mopacformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\mpdformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\mpqcformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\nwchemformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\obconversion.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\oberror.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\obiter.cpp
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

SOURCE=..\..\src\formats\pcmodelformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\pdbformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\phmodel.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\povrayformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\PQSformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\xml\pubchem.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\qchemformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\rand.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\reportformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\residue.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\ring.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\rotamer.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\rotor.cpp
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

SOURCE=..\..\src\tokenst.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\transform.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\turbomoleformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\typer.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\unichemformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\math\vector3.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\viewmolformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\xml\xcmlformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\xedformat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\xml\xml.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\xml\xmlformat.cpp
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

SOURCE=..\babelconfig.h
# End Source File
# Begin Source File

SOURCE=..\..\src\base.h
# End Source File
# Begin Source File

SOURCE="C:\Program Files\Microsoft Visual Studio\VC98\Include\BASETSD.H"
# End Source File
# Begin Source File

SOURCE=..\..\src\bitvec.h
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

SOURCE=..\..\src\crk.h
# End Source File
# Begin Source File

SOURCE=..\..\src\data.h
# End Source File
# Begin Source File

SOURCE=..\..\src\dlhandler.h
# End Source File
# Begin Source File

SOURCE=..\..\src\extable.h
# End Source File
# Begin Source File

SOURCE=..\..\src\fingerprint.h
# End Source File
# Begin Source File

SOURCE=..\..\src\generic.h
# End Source File
# Begin Source File

SOURCE=..\..\src\grid.h
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\inchi_api.h
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

SOURCE=..\..\src\obconversion.h
# End Source File
# Begin Source File

SOURCE=..\..\src\oberror.h
# End Source File
# Begin Source File

SOURCE=..\..\src\obifstream.h
# End Source File
# Begin Source File

SOURCE=..\..\src\obiter.h
# End Source File
# Begin Source File

SOURCE=..\..\src\obmolecformat.h
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

SOURCE=..\..\src\reaction.h
# End Source File
# Begin Source File

SOURCE=..\..\src\ring.h
# End Source File
# Begin Source File

SOURCE=..\..\src\rotamer.h
# End Source File
# Begin Source File

SOURCE=..\..\src\rotor.h
# End Source File
# Begin Source File

SOURCE=..\..\src\formats\smi.h
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
# Begin Source File

SOURCE=..\..\src\formats\xml\xml.h
# End Source File
# End Group
# Begin Source File

SOURCE=..\libxml2.lib
# End Source File
# Begin Source File

SOURCE=..\zdll.lib
# End Source File
# End Target
# End Project
