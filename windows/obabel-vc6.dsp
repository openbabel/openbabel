# Microsoft Developer Studio Project File - Name="obabel" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=obabel - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "obabel.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "obabel.mak" CFG="obabel - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "obabel - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "obabel - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "obabel - Win32 Release"

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
# ADD CPP /nologo /W3 /GX /O2 /I ".. ..\src" /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# SUBTRACT CPP /X
# ADD BASE RSC /l 0x809 /d "NDEBUG"
# ADD RSC /l 0x809 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386 /out:"../src/obabel.exe"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

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
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /I "..\src" /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /FR /FD /GZ /c
# SUBTRACT CPP /X /YX
# ADD BASE RSC /l 0x809 /d "_DEBUG"
# ADD RSC /l 0x809 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "obabel - Win32 Release"
# Name "obabel - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\src\alchemy.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\amber.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\atom.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\balst.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\base.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\bgf.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\binary.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\bitvec.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\bond.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\box.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\c3d.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\cacao.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\cache.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\car.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\ccc.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\chains.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\chdrw.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\chemtool.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\chiral.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\cml.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\crk.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\csr.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\cssr.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\data.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\dmol.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\feat.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\fh.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\fileformat.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\gamess.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\gaussian.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\generic.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\ghemical.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\grid.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\gromos96.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\hin.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\jaguar.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\main.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\matrix.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\math\matrix3x3.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\mdl.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\mm3.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\mmod.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\mol.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\mol2.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\molchrg.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\molvector.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\mopac.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\mpqc.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\nwchem.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\oberror.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\obutil.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\parsmart.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\parsmi.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\patty.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\pdb.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\phmodel.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\povray.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\pqs.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\qchem.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\rand.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\report.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\residue.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\ring.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\rotor.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\shelx.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\smi.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\snprintf.c

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\strncasecmp.c

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\tinker.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\tokenst.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\typer.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\unichem.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\math\vector3.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\viewmol.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\xed.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\xyz.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\zindo.cpp

!IF  "$(CFG)" == "obabel - Win32 Release"

# ADD CPP /I "..\src"
# SUBTRACT CPP /I ".. ..\src"

!ELSEIF  "$(CFG)" == "obabel - Win32 Debug"

!ENDIF 

# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\src\aromatic.h
# End Source File
# Begin Source File

SOURCE=..\src\atomtyp.h
# End Source File
# Begin Source File

SOURCE=.\babelconfig.h
# End Source File
# Begin Source File

SOURCE=..\src\base.h
# End Source File
# Begin Source File

SOURCE=..\src\binary.h
# End Source File
# Begin Source File

SOURCE=..\src\bitvec.h
# End Source File
# Begin Source File

SOURCE=..\src\bondtyp.h
# End Source File
# Begin Source File

SOURCE=..\src\bondtyper.h
# End Source File
# Begin Source File

SOURCE=..\src\chains.h
# End Source File
# Begin Source File

SOURCE=..\src\chiral.h
# End Source File
# Begin Source File

SOURCE=..\src\crk.h
# End Source File
# Begin Source File

SOURCE=..\src\data.h
# End Source File
# Begin Source File

SOURCE=..\src\element.h
# End Source File
# Begin Source File

SOURCE=..\src\extable.h
# End Source File
# Begin Source File

SOURCE=..\src\fileformat.h
# End Source File
# Begin Source File

SOURCE=..\src\generic.h
# End Source File
# Begin Source File

SOURCE=..\src\grid.h
# End Source File
# Begin Source File

SOURCE=..\src\isotope.h
# End Source File
# Begin Source File

SOURCE=..\src\matrix.h
# End Source File
# Begin Source File

SOURCE=..\src\math\matrix3x3.h
# End Source File
# Begin Source File

SOURCE=..\src\mol.h
# End Source File
# Begin Source File

SOURCE=..\src\molchrg.h
# End Source File
# Begin Source File

SOURCE=..\src\molvector.h
# End Source File
# Begin Source File

SOURCE=..\src\oberror.h
# End Source File
# Begin Source File

SOURCE=..\src\obifstream.h
# End Source File
# Begin Source File

SOURCE=..\src\obutil.h
# End Source File
# Begin Source File

SOURCE=..\src\parsmart.h
# End Source File
# Begin Source File

SOURCE=..\src\patty.h
# End Source File
# Begin Source File

SOURCE=..\src\phmodel.h
# End Source File
# Begin Source File

SOURCE=..\src\phmodeldata.h
# End Source File
# Begin Source File

SOURCE=..\src\resdata.h
# End Source File
# Begin Source File

SOURCE=..\src\ring.h
# End Source File
# Begin Source File

SOURCE=..\src\rotor.h
# End Source File
# Begin Source File

SOURCE=..\src\smi.h
# End Source File
# Begin Source File

SOURCE=..\src\snprintf.h
# End Source File
# Begin Source File

SOURCE=..\src\typer.h
# End Source File
# Begin Source File

SOURCE=..\src\types.h
# End Source File
# Begin Source File

SOURCE=..\src\math\vector3.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# Begin Source File

SOURCE=".\How to Compile with VC++.txt"
# End Source File
# End Target
# End Project
