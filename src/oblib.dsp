# Microsoft Developer Studio Project File - Name="oblib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Dynamic-Link Library" 0x0102

CFG=oblib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "oblib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "oblib.mak" CFG="oblib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "oblib - Win32 Release" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE "oblib - Win32 Debug" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "oblib - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MT /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "oblib_EXPORTS" /YX /FD /c
# ADD CPP /nologo /MT /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "oblib_EXPORTS" /YX /FD /c
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /machine:I386

!ELSEIF  "$(CFG)" == "oblib - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MTd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "oblib_EXPORTS" /YX /FD /GZ /c
# ADD CPP /nologo /MTd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "oblib_EXPORTS" /YX /FD /GZ /c
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "oblib - Win32 Release"
# Name "oblib - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\alchemy.cpp
# End Source File
# Begin Source File

SOURCE=.\amber.cpp
# End Source File
# Begin Source File

SOURCE=.\atom.cpp
# End Source File
# Begin Source File

SOURCE=.\balst.cpp
# End Source File
# Begin Source File

SOURCE=.\base.cpp
# End Source File
# Begin Source File

SOURCE=.\bgf.cpp
# End Source File
# Begin Source File

SOURCE=.\binary.cpp
# End Source File
# Begin Source File

SOURCE=.\binary_io.cpp
# End Source File
# Begin Source File

SOURCE=.\bitgrid.cpp
# End Source File
# Begin Source File

SOURCE=.\bitvec.cpp
# End Source File
# Begin Source File

SOURCE=.\bond.cpp
# End Source File
# Begin Source File

SOURCE=.\box.cpp
# End Source File
# Begin Source File

SOURCE=.\c3d.cpp
# End Source File
# Begin Source File

SOURCE=.\cacao.cpp
# End Source File
# Begin Source File

SOURCE=.\cache.cpp
# End Source File
# Begin Source File

SOURCE=.\car.cpp
# End Source File
# Begin Source File

SOURCE=.\ccc.cpp
# End Source File
# Begin Source File

SOURCE=.\chains.cpp
# End Source File
# Begin Source File

SOURCE=.\chdrw.cpp
# End Source File
# Begin Source File

SOURCE=.\chiral.cpp
# End Source File
# Begin Source File

SOURCE=.\commandline.cpp
# End Source File
# Begin Source File

SOURCE=.\csr.cpp
# End Source File
# Begin Source File

SOURCE=.\cssr.cpp
# End Source File
# Begin Source File

SOURCE=.\ctransform.cpp
# End Source File
# Begin Source File

SOURCE=.\cwrap.cpp
# End Source File
# Begin Source File

SOURCE=.\data.cpp
# End Source File
# Begin Source File

SOURCE=.\dmol.cpp
# End Source File
# Begin Source File

SOURCE=.\feat.cpp
# End Source File
# Begin Source File

SOURCE=.\fh.cpp
# End Source File
# Begin Source File

SOURCE=.\fileformat.cpp
# End Source File
# Begin Source File

SOURCE=.\gamess.cpp
# End Source File
# Begin Source File

SOURCE=.\gaussian.cpp
# End Source File
# Begin Source File

SOURCE=.\generic.cpp
# End Source File
# Begin Source File

SOURCE=.\ghemical.cpp
# End Source File
# Begin Source File

SOURCE=.\grid.cpp
# End Source File
# Begin Source File

SOURCE=.\gromos96.cpp
# End Source File
# Begin Source File

SOURCE=.\hin.cpp
# End Source File
# Begin Source File

SOURCE=.\jaguar.cpp
# End Source File
# Begin Source File

SOURCE=.\matrix.cpp
# End Source File
# Begin Source File

SOURCE=.\mdl.cpp
# End Source File
# Begin Source File

SOURCE=.\mmod.cpp
# End Source File
# Begin Source File

SOURCE=.\mol.cpp
# End Source File
# Begin Source File

SOURCE=.\mol2.cpp
# End Source File
# Begin Source File

SOURCE=.\molchrg.cpp
# End Source File
# Begin Source File

SOURCE=.\molvector.cpp
# End Source File
# Begin Source File

SOURCE=.\mopac.cpp
# End Source File
# Begin Source File

SOURCE=.\mpqc.cpp
# End Source File
# Begin Source File

SOURCE=.\newmain.cpp
# End Source File
# Begin Source File

SOURCE=.\nwchem.cpp
# End Source File
# Begin Source File

SOURCE=.\obutil.cpp
# End Source File
# Begin Source File

SOURCE=.\parsmart.cpp
# End Source File
# Begin Source File

SOURCE=.\parsmi.cpp
# End Source File
# Begin Source File

SOURCE=.\patty.cpp
# End Source File
# Begin Source File

SOURCE=.\pdb.cpp
# End Source File
# Begin Source File

SOURCE=.\phmodel.cpp
# End Source File
# Begin Source File

SOURCE=.\qchem.cpp
# End Source File
# Begin Source File

SOURCE=.\quat.c
# End Source File
# Begin Source File

SOURCE=.\rand.cpp
# End Source File
# Begin Source File

SOURCE=.\report.cpp
# End Source File
# Begin Source File

SOURCE=.\ring.cpp
# End Source File
# Begin Source File

SOURCE=.\rotor.cpp
# End Source File
# Begin Source File

SOURCE=.\smarts.cpp
# End Source File
# Begin Source File

SOURCE=.\smi.cpp
# End Source File
# Begin Source File

SOURCE=.\tinker.cpp
# End Source File
# Begin Source File

SOURCE=.\tokenst.cpp
# End Source File
# Begin Source File

SOURCE=.\typer.cpp
# End Source File
# Begin Source File

SOURCE=.\unichem.cpp
# End Source File
# Begin Source File

SOURCE=.\Vector.cpp
# End Source File
# Begin Source File

SOURCE=.\xed.cpp
# End Source File
# Begin Source File

SOURCE=.\xyz.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\aromatic.h
# End Source File
# Begin Source File

SOURCE=.\atomtyp.h
# End Source File
# Begin Source File

SOURCE=.\base.h
# End Source File
# Begin Source File

SOURCE=.\binary.h
# End Source File
# Begin Source File

SOURCE=.\binary_io.h
# End Source File
# Begin Source File

SOURCE=.\bitgrid.h
# End Source File
# Begin Source File

SOURCE=.\bitvec.h
# End Source File
# Begin Source File

SOURCE=.\chains.h
# End Source File
# Begin Source File

SOURCE=.\chiral.h
# End Source File
# Begin Source File

SOURCE=.\commandline.h
# End Source File
# Begin Source File

SOURCE=.\ctransform.h
# End Source File
# Begin Source File

SOURCE=.\cwrap.h
# End Source File
# Begin Source File

SOURCE=.\data.h
# End Source File
# Begin Source File

SOURCE=.\element.h
# End Source File
# Begin Source File

SOURCE=.\extable.h
# End Source File
# Begin Source File

SOURCE=.\fileformat.h
# End Source File
# Begin Source File

SOURCE=.\generic.h
# End Source File
# Begin Source File

SOURCE=.\grid.h
# End Source File
# Begin Source File

SOURCE=.\matrix.h
# End Source File
# Begin Source File

SOURCE=.\mol.h
# End Source File
# Begin Source File

SOURCE=.\molchrg.h
# End Source File
# Begin Source File

SOURCE=.\molvector.h
# End Source File
# Begin Source File

SOURCE=.\obifstream.h
# End Source File
# Begin Source File

SOURCE=.\obutil.h
# End Source File
# Begin Source File

SOURCE=.\parsmart.h
# End Source File
# Begin Source File

SOURCE=.\patty.h
# End Source File
# Begin Source File

SOURCE=.\phmodel.h
# End Source File
# Begin Source File

SOURCE=.\phmodeldata.h
# End Source File
# Begin Source File

SOURCE=.\resdata.h
# End Source File
# Begin Source File

SOURCE=.\Resource.h
# End Source File
# Begin Source File

SOURCE=.\ring.h
# End Source File
# Begin Source File

SOURCE=.\rotor.h
# End Source File
# Begin Source File

SOURCE=.\smarts.h
# End Source File
# Begin Source File

SOURCE=.\smi.h
# End Source File
# Begin Source File

SOURCE=.\typer.h
# End Source File
# Begin Source File

SOURCE=.\types.h
# End Source File
# Begin Source File

SOURCE=.\Vector.h
# End Source File
# Begin Source File

SOURCE=.\winversion.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
