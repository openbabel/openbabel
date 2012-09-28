;NSIS Modern User Interface
;Installer for Open Babel

;--------------------------------
;Include Modern UI

  !include "MUI.nsh"

; ======================================
; The main part of this script starts at
; around line 455
; ======================================

;--------------------------------
;AddToPath functions taken from "Path Manipulation" on NSIS wiki

!ifndef _AddToPath_nsh
!define _AddToPath_nsh

!verbose 3
!include "WinMessages.NSH"
!verbose 4
 
!ifndef WriteEnvStr_RegKey
  !ifdef ALL_USERS
    !define WriteEnvStr_RegKey \
       'HKLM "SYSTEM\CurrentControlSet\Control\Session Manager\Environment"'
  !else
    !define WriteEnvStr_RegKey 'HKCU "Environment"'
  !endif
!endif
 
; AddToPath - Adds the given dir to the search path.
;        Input - head of the stack
;        Note - Win9x systems requires reboot
 
Function AddToPath
  Exch $0
  Push $1
  Push $2
  Push $3
 
  # don't add if the path doesn't exist
  IfFileExists "$0\*.*" "" AddToPath_done
 
  ReadEnvStr $1 PATH
  Push "$1;"
  Push "$0;"
  Call StrStr
  Pop $2
  StrCmp $2 "" "" AddToPath_done
  Push "$1;"
  Push "$0\;"
  Call StrStr
  Pop $2
  StrCmp $2 "" "" AddToPath_done
  GetFullPathName /SHORT $3 $0
  Push "$1;"
  Push "$3;"
  Call StrStr
  Pop $2
  StrCmp $2 "" "" AddToPath_done
  Push "$1;"
  Push "$3\;"
  Call StrStr
  Pop $2
  StrCmp $2 "" "" AddToPath_done
 
  Call IsNT
  Pop $1
  StrCmp $1 1 AddToPath_NT
    ; Not on NT
    StrCpy $1 $WINDIR 2
    FileOpen $1 "$1\autoexec.bat" a
    FileSeek $1 -1 END
    FileReadByte $1 $2
    IntCmp $2 26 0 +2 +2 # DOS EOF
      FileSeek $1 -1 END # write over EOF
    FileWrite $1 "$\r$\nSET PATH=%PATH%;$3$\r$\n"
    FileClose $1
    SetRebootFlag true
    Goto AddToPath_done
 
  AddToPath_NT:
    ReadRegStr $1 ${WriteEnvStr_RegKey} "PATH"
    StrCmp $1 "" AddToPath_NTdoIt
      Push $1
      Call Trim
      Pop $1
      ;StrCpy $0 "$1;$0" Edited by Noel - add to front of PATH
      StrCpy $0 "$0;$1"
    AddToPath_NTdoIt:
      WriteRegExpandStr ${WriteEnvStr_RegKey} "PATH" $0
      SendMessage ${HWND_BROADCAST} ${WM_WININICHANGE} 0 "STR:Environment" /TIMEOUT=5000
 
  AddToPath_done:
    Pop $3
    Pop $2
    Pop $1
    Pop $0
FunctionEnd
 
; RemoveFromPath - Remove a given dir from the path
;     Input: head of the stack
 
Function un.RemoveFromPath
  Exch $0
  Push $1
  Push $2
  Push $3
  Push $4
  Push $5
  Push $6
 
  IntFmt $6 "%c" 26 # DOS EOF
 
  Call un.IsNT
  Pop $1
  StrCmp $1 1 unRemoveFromPath_NT
    ; Not on NT
    StrCpy $1 $WINDIR 2
    FileOpen $1 "$1\autoexec.bat" r
    GetTempFileName $4
    FileOpen $2 $4 w
    GetFullPathName /SHORT $0 $0
    StrCpy $0 "SET PATH=%PATH%;$0"
    Goto unRemoveFromPath_dosLoop
 
    unRemoveFromPath_dosLoop:
      FileRead $1 $3
      StrCpy $5 $3 1 -1 # read last char
      StrCmp $5 $6 0 +2 # if DOS EOF
        StrCpy $3 $3 -1 # remove DOS EOF so we can compare
      StrCmp $3 "$0$\r$\n" unRemoveFromPath_dosLoopRemoveLine
      StrCmp $3 "$0$\n" unRemoveFromPath_dosLoopRemoveLine
      StrCmp $3 "$0" unRemoveFromPath_dosLoopRemoveLine
      StrCmp $3 "" unRemoveFromPath_dosLoopEnd
      FileWrite $2 $3
      Goto unRemoveFromPath_dosLoop
      unRemoveFromPath_dosLoopRemoveLine:
        SetRebootFlag true
        Goto unRemoveFromPath_dosLoop
 
    unRemoveFromPath_dosLoopEnd:
      FileClose $2
      FileClose $1
      StrCpy $1 $WINDIR 2
      Delete "$1\autoexec.bat"
      CopyFiles /SILENT $4 "$1\autoexec.bat"
      Delete $4
      Goto unRemoveFromPath_done
 
  unRemoveFromPath_NT:
    ReadRegStr $1 ${WriteEnvStr_RegKey} "PATH"
    StrCpy $5 $1 1 -1 # copy last char
    StrCmp $5 ";" +2 # if last char != ;
      StrCpy $1 "$1;" # append ;
    Push $1
    Push "$0;"
    Call un.StrStr ; Find `$0;` in $1
    Pop $2 ; pos of our dir
    StrCmp $2 "" unRemoveFromPath_done
      ; else, it is in path
      # $0 - path to add
      # $1 - path var
      StrLen $3 "$0;"
      StrLen $4 $2
      StrCpy $5 $1 -$4 # $5 is now the part before the path to remove
      StrCpy $6 $2 "" $3 # $6 is now the part after the path to remove
      StrCpy $3 $5$6
 
      StrCpy $5 $3 1 -1 # copy last char
      StrCmp $5 ";" 0 +2 # if last char == ;
        StrCpy $3 $3 -1 # remove last char
 
      WriteRegExpandStr ${WriteEnvStr_RegKey} "PATH" $3
      SendMessage ${HWND_BROADCAST} ${WM_WININICHANGE} 0 "STR:Environment" /TIMEOUT=5000
 
  unRemoveFromPath_done:
    Pop $6
    Pop $5
    Pop $4
    Pop $3
    Pop $2
    Pop $1
    Pop $0
FunctionEnd
 
 
 
; AddToEnvVar - Adds the given value to the given environment var
;        Input - head of the stack $0 environement variable $1=value to add
;        Note - Win9x systems requires reboot
 
Function AddToEnvVar
 
  Exch $1 ; $1 has environment variable value
  Exch
  Exch $0 ; $0 has environment variable name
 
  DetailPrint "Adding $1 to $0"
  Push $2
  Push $3
  Push $4
 
 
  ReadEnvStr $2 $0
  Push "$2;"
  Push "$1;"
  Call StrStr
  Pop $3
  StrCmp $3 "" "" AddToEnvVar_done
 
  Push "$2;"
  Push "$1\;"
  Call StrStr
  Pop $3
  StrCmp $3 "" "" AddToEnvVar_done
  
 
  Call IsNT
  Pop $2
  StrCmp $2 1 AddToEnvVar_NT
    ; Not on NT
    StrCpy $2 $WINDIR 2
    FileOpen $2 "$2\autoexec.bat" a
    FileSeek $2 -1 END
    FileReadByte $2 $3
    IntCmp $3 26 0 +2 +2 # DOS EOF
      FileSeek $2 -1 END # write over EOF
    FileWrite $2 "$\r$\nSET $0=%$0%;$4$\r$\n"
    FileClose $2
    SetRebootFlag true
    Goto AddToEnvVar_done
 
  AddToEnvVar_NT:
    ReadRegStr $2 ${WriteEnvStr_RegKey} $0
    StrCpy $3 $2 1 -1 # copy last char
    StrCmp $3 ";" 0 +2 # if last char == ;
      StrCpy $2 $2 -1 # remove last char
    StrCmp $2 "" AddToEnvVar_NTdoIt
      ; StrCpy $1 "$2;$1" Noel - Don't add - set the value
      StrCpy $1 "$2"
    AddToEnvVar_NTdoIt:
      WriteRegExpandStr ${WriteEnvStr_RegKey} $0 $1
      SendMessage ${HWND_BROADCAST} ${WM_WININICHANGE} 0 "STR:Environment" /TIMEOUT=5000
 
  AddToEnvVar_done:
    Pop $4
    Pop $3
    Pop $2
    Pop $0
    Pop $1
 
FunctionEnd
 
; RemoveFromEnvVar - Remove a given value from a environment var
;     Input: head of the stack
 
Function un.RemoveFromEnvVar
 
  Exch $1 ; $1 has environment variable value
  Exch
  Exch $0 ; $0 has environment variable name
 
  DetailPrint "Removing $1 from $0"
  Push $2
  Push $3
  Push $4
  Push $5
  Push $6
  Push $7
 
  IntFmt $7 "%c" 26 # DOS EOF
 
  Call un.IsNT
  Pop $2
  StrCmp $2 1 unRemoveFromEnvVar_NT
    ; Not on NT
    StrCpy $2 $WINDIR 2
    FileOpen $2 "$2\autoexec.bat" r
    GetTempFileName $5
    FileOpen $3 $5 w
    GetFullPathName /SHORT $1 $1
    StrCpy $1 "SET $0=%$0%;$1"
    Goto unRemoveFromEnvVar_dosLoop
 
    unRemoveFromEnvVar_dosLoop:
      FileRead $2 $4
      StrCpy $6 $4 1 -1 # read last char
      StrCmp $6 $7 0 +2 # if DOS EOF
        StrCpy $4 $4 -1 # remove DOS EOF so we can compare
      StrCmp $4 "$1$\r$\n" unRemoveFromEnvVar_dosLoopRemoveLine
      StrCmp $4 "$1$\n" unRemoveFromEnvVar_dosLoopRemoveLine
      StrCmp $4 "$1" unRemoveFromEnvVar_dosLoopRemoveLine
      StrCmp $4 "" unRemoveFromEnvVar_dosLoopEnd
      FileWrite $3 $4
      Goto unRemoveFromEnvVar_dosLoop
      unRemoveFromEnvVar_dosLoopRemoveLine:
        SetRebootFlag true
        Goto unRemoveFromEnvVar_dosLoop
 
    unRemoveFromEnvVar_dosLoopEnd:
      FileClose $3
      FileClose $2
      StrCpy $2 $WINDIR 2
      Delete "$2\autoexec.bat"
      CopyFiles /SILENT $5 "$2\autoexec.bat"
      Delete $5
      Goto unRemoveFromEnvVar_done
 
  unRemoveFromEnvVar_NT:
    ReadRegStr $2 ${WriteEnvStr_RegKey} $0
    StrCpy $6 $2 1 -1 # copy last char
    StrCmp $6 ";" +2 # if last char != ;
      StrCpy $2 "$2;" # append ;
    Push $2
    Push "$1;"
    Call un.StrStr ; Find `$1;` in $2
    Pop $3 ; pos of our dir
    StrCmp $3 "" unRemoveFromEnvVar_done
      ; else, it is in path
      # $1 - path to add
      # $2 - path var
      StrLen $4 "$1;"
      StrLen $5 $3
      StrCpy $6 $2 -$5 # $6 is now the part before the path to remove
      StrCpy $7 $3 "" $4 # $7 is now the part after the path to remove
      StrCpy $4 $6$7
 
      StrCpy $6 $4 1 -1 # copy last char
      StrCmp $6 ";" 0 +2 # if last char == ;
      StrCpy $4 $4 -1 # remove last char
 
      WriteRegExpandStr ${WriteEnvStr_RegKey} $0 $4
      SendMessage ${HWND_BROADCAST} ${WM_WININICHANGE} 0 "STR:Environment" /TIMEOUT=5000
 
  unRemoveFromEnvVar_done:
    Pop $7
    Pop $6
    Pop $5
    Pop $4
    Pop $3
    Pop $2
    Pop $1
    Pop $0
FunctionEnd
 
 
 
 
!ifndef IsNT_KiCHiK
!define IsNT_KiCHiK
 
###########################################
#            Utility Functions            #
###########################################
 
; IsNT
; no input
; output, top of the stack = 1 if NT or 0 if not
;
; Usage:
;   Call IsNT
;   Pop $R0
;  ($R0 at this point is 1 or 0)
 
!macro IsNT un
Function ${un}IsNT
  Push $0
  ReadRegStr $0 HKLM "SOFTWARE\Microsoft\Windows NT\CurrentVersion" CurrentVersion
  StrCmp $0 "" 0 IsNT_yes
  ; we are not NT.
  Pop $0
  Push 0
  Return
 
  IsNT_yes:
    ; NT!!!
    Pop $0
    Push 1
FunctionEnd
!macroend
!insertmacro IsNT ""
!insertmacro IsNT "un."
 
!endif ; IsNT_KiCHiK
 
; StrStr
; input, top of stack = string to search for
;        top of stack-1 = string to search in
; output, top of stack (replaces with the portion of the string remaining)
; modifies no other variables.
;
; Usage:
;   Push "this is a long ass string"
;   Push "ass"
;   Call StrStr
;   Pop $R0
;  ($R0 at this point is "ass string")
 
!macro StrStr un
Function ${un}StrStr
Exch $R1 ; st=haystack,old$R1, $R1=needle
  Exch    ; st=old$R1,haystack
  Exch $R2 ; st=old$R1,old$R2, $R2=haystack
  Push $R3
  Push $R4
  Push $R5
  StrLen $R3 $R1
  StrCpy $R4 0
  ; $R1=needle
  ; $R2=haystack
  ; $R3=len(needle)
  ; $R4=cnt
  ; $R5=tmp
  loop:
    StrCpy $R5 $R2 $R3 $R4
    StrCmp $R5 $R1 done
    StrCmp $R5 "" done
    IntOp $R4 $R4 + 1
    Goto loop
done:
  StrCpy $R1 $R2 "" $R4
  Pop $R5
  Pop $R4
  Pop $R3
  Pop $R2
  Exch $R1
FunctionEnd
!macroend
!insertmacro StrStr ""
!insertmacro StrStr "un."
 
!endif ; _AddToPath_nsh
 
Function Trim ; Added by Pelaca
	Exch $R1
	Push $R2
Loop:
	StrCpy $R2 "$R1" 1 -1
	StrCmp "$R2" " " RTrim
	StrCmp "$R2" "$\n" RTrim
	StrCmp "$R2" "$\r" RTrim
	StrCmp "$R2" ";" RTrim
	GoTo Done
RTrim:	
	StrCpy $R1 "$R1" -1
	Goto Loop
Done:
	Pop $R2
	Exch $R1
FunctionEnd

;--------------------------------
;General

  ;OpenBabel version
  !define OBVersion 2.3.2

  ;Name and file
  Name "OpenBabel ${OBVERSION}"
  OutFile "OpenBabel${OBVERSION}_Windows_Installer.exe"

  ;Default installation folder 
  InstallDir "$PROGRAMFILES\OpenBabel-${OBVERSION}"
  
  ;Get installation folder from registry if available
  InstallDirRegKey HKCU "Software\OpenBabel ${OBVERSION}" ""

;--------------------------------
;Variables

  Var MUI_TEMP
  Var STARTMENU_FOLDER

;--------------------------------
;Interface Settings

  ;!define MUI_ICON "babel.ico"
  ;!define MUI_UNICON "babel.ico"
  !define MUI_ABORTWARNING
  !define MUI_HEADERIMAGE
  !define MUI_HEADERIMAGE_BITMAP "logo_small.bmp"
  !define MUI_WELCOMEFINISHPAGE_BITMAP "logo_big.bmp"
  !define MUI_FINISHPAGE_RUN "$INSTDIR/OBGUI.exe"
;  !define MUI_FINISHPAGE_SHOWREADME "http://openbabel.org/wiki/OpenBabelGUI"

;--------------------------------
;Pages

  !insertmacro MUI_PAGE_WELCOME
  !insertmacro MUI_PAGE_LICENSE "../../COPYING"
  !insertmacro MUI_PAGE_DIRECTORY

  ;Start Menu Folder Page Configuration
  !define MUI_STARTMENUPAGE_REGISTRY_ROOT "HKCU" 
  !define MUI_STARTMENUPAGE_REGISTRY_KEY "Software\OpenBabel ${OBVERSION}" 
  !define MUI_STARTMENUPAGE_REGISTRY_VALUENAME "Start Menu Folder"
  !insertmacro MUI_PAGE_STARTMENU Application $STARTMENU_FOLDER

  !insertmacro MUI_PAGE_INSTFILES
  !insertmacro MUI_PAGE_FINISH
  
  !insertmacro MUI_UNPAGE_CONFIRM
  !insertmacro MUI_UNPAGE_INSTFILES
  
;--------------------------------
;Languages
 
  !insertmacro MUI_LANGUAGE "English"

;--------------------------------
;Installer Sections

Section "Dummy Section" SecDummy

  SetShellVarContext all

  ;The data will be in (writable) subfolder of %APPDATA%, if it exists
  Var /GLOBAL DataBase
  ClearErrors
  ReadEnvStr $0 APPDATA
  IfErrors 0 +3    
  StrCpy $DataBase "$INSTDIR"
  Goto +2
  StrCpy $DataBase "$%APPDATA%\OpenBabel-${OBVERSION}"
   
  SetOutPath "$DataBase\data"
  File /r /x .svn /x *.h ..\..\data\*.*

  SetOutPath "$DataBase\doc"
  File ..\..\doc\OpenBabelGUI.html
  File ToolsPrograms.txt

  SetOutPath "$INSTDIR"
  File /oname=License.txt ..\..\COPYING
  File ..\sdf.bat
  File ..\obdepict.bat
  
  File /r /x test_*.* ..\build\bin\Release\*.exe
  File /r ..\build\bin\Release\*.obf
  File ..\build\bin\Release\openbabel-2.dll
  
  File vcredist_x86.exe

  ;Java and CSharp bindings
  File ..\..\scripts\java\openbabel.jar
  File ..\build\bin\Release\openbabel_java.dll
  File ..\build\bin\Release\openbabel_csharp.dll
  File ..\..\scripts\csharp\OBDotNet.dll

  File ..\libs\i386\*.dll

  SetOutPath "$DataBase\examples"
  File /r /x .svn ExampleFiles\*.*
  
  ;Files for user to build own C++ programs using precompiled OpenBabel
  ;Not working correctly
  ;SetOutPath "$DataBase\obbuild"
  ;File ..\build\src\Release\openbabel-2.lib   
  ;File /r /x .svn obbuild\*.*
  ;SetOutPath "$DataBase\obbuild\openbabel"
  ;File /r /x .svn ..\..\include\openbabel\*.*
  ;File ..\build\include\openbabel\babelconfig.h  
  
  ;Store installation folder
  WriteRegStr HKCU "Software\OpenBabel ${OBVERSION}" "" $INSTDIR
  
  ;Install VC++ 2010 redistributable
  ExecWait '"$INSTDIR/vcredist_x86.exe" /q:a'
  Delete "$INSTDIR\vcredist_x86.exe"

  ;Create uninstaller
  WriteUninstaller "$INSTDIR\Uninstall.exe"

  ;Create shortcuts
  !insertmacro MUI_STARTMENU_WRITE_BEGIN Application
    CreateDirectory "$SMPROGRAMS\$STARTMENU_FOLDER"
    CreateShortCut "$SMPROGRAMS\$STARTMENU_FOLDER\Open Babel GUI.lnk" "$INSTDIR\OBGUI.exe"
    CreateShortCut "$SMPROGRAMS\$STARTMENU_FOLDER\Guide to using Open Babel GUI.lnk" "$DataBase\doc\OpenBabelGUI.html"
    CreateShortCut "$SMPROGRAMS\$STARTMENU_FOLDER\Guide to using babel (web).lnk" "http://openbabel.org/docs/2.3.1"
    CreateShortCut "$SMPROGRAMS\$STARTMENU_FOLDER\Uninstall.lnk" "$INSTDIR\Uninstall.exe"
  !insertmacro MUI_STARTMENU_WRITE_END

  ;Add to PATH
  Push $INSTDIR
  Call AddToPath

  ; Set env var BABEL_DATADIR
  ; First remove any existing value
  DeleteRegValue        HKCU "Environment" "BABEL_DATADIR"

  ; Old way: WriteRegStr HKCU "Environment" "BABEL_DATADIR" "$INSTDIR" 
  ; New way: Works immediately
  Push "BABEL_DATADIR"
  Push "$DataBase\data"
  Call AddToEnvVar
  
  ; Entry for Add/Remove Programs
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\OpenBabel-${OBVERSION}" "DisplayName" "OpenBabel-${OBVERSION}"
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\OpenBabel-${OBVERSION}" "UninstallString" "$INSTDIR\Uninstall.exe"
 
  ;Add Firefox to PATH so GUI can use it for displaying structures
  EnumRegKey $6 HKLM "SOFTWARE\Mozilla\Mozilla Firefox" 0
  ReadRegStr $7 HKLM "SOFTWARE\Mozilla\Mozilla Firefox\$6\Main" "Install Directory"
  Push $7
  Call AddToPath

  ;Put OBDepict.bat in context menu
  WriteRegStr HKCR "*\shell\OBDepict\command" "" "$INSTDIR\obdepict.bat %1"
  
  ; Convenience shortcut to OBDepict
  ;CreateShortCut "$INSTDIR\OBDepict.lnk" "$INSTDIR\obdepict.bat" "" "" 0 SW_SHOWMINIMIZED 
   
SectionEnd

;--------------------------------
;Descriptions

  ;Language strings
  LangString DESC_SecDummy ${LANG_ENGLISH} "A test section."

  ;Assign language strings to sections
  !insertmacro MUI_FUNCTION_DESCRIPTION_BEGIN
    !insertmacro MUI_DESCRIPTION_TEXT ${SecDummy} $(DESC_SecDummy)
  !insertmacro MUI_FUNCTION_DESCRIPTION_END

;--------------------------------
;Uninstaller Section

Section "Uninstall"
  Delete "$INSTDIR\ob*.exe"
  Delete "$INSTDIR\babel.exe"
  Delete "$INSTDIR\obdepict.bat"
  Delete "$INSTDIR\OBDepict.lnk"
  Delete "$INSTDIR\sdf.bat"
  Delete "$INSTDIR\*.obf"
  
  ;May not be installed, but it doesn't matter.
  Delete "$INSTDIR\Sieve.exe"
  
  Delete "$INSTDIR\openbabel-2.dll"
  Delete "$INSTDIR\iconv.dll"
  Delete "$INSTDIR\libinchi.dll"
  Delete "$INSTDIR\libxml2.dll"
  Delete "$INSTDIR\zlib1.dll"
  Delete "$INSTDIR\xdr-0.dll"
  Delete "$INSTDIR\License.txt"
  Delete "$INSTDIR\OBGUI.exe"
  Delete "$INSTDIR\OpenBabelGUI.html"
  Delete "$INSTDIR\openbabel.jar"
  Delete "$INSTDIR\openbabel_java.dll"
  Delete "$INSTDIR\OBDotNet.dll"
  Delete "$INSTDIR\openbabel_csharp.dll"
  Delete "$INSTDIR\libpng14-14.dll"
  Delete "$INSTDIR\freetype6.dll"
  Delete "$INSTDIR\libcairo-2.dll"
  Delete "$INSTDIR\libexpat-1.dll"
  Delete "$INSTDIR\libfontconfig-1.dll"

  Delete "$INSTDIR\Uninstall.exe"

  RMDir "$INSTDIR"
  
  SetShellVarContext all
  !insertmacro MUI_STARTMENU_GETFOLDER Application $MUI_TEMP
    
  ;Delete "$SMPROGRAMS\$MUI_TEMP\*.lnk"
  RMDir /r "$SMPROGRAMS\$MUI_TEMP"

  ;Remove from PATH
  push $INSTDIR
  Call un.RemoveFromPath
  
  ;Delete empty start menu parent diretories
  ;StrCpy $MUI_TEMP "$SMPROGRAMS\$MUI_TEMP"
 
;  startMenuDeleteLoop:
;	ClearErrors
;    RMDir $MUI_TEMP
;    GetFullPathName $MUI_TEMP "$MUI_TEMP\.."
;    
;    IfErrors startMenuDeleteLoopDone
;  
;    StrCmp $MUI_TEMP $SMPROGRAMS startMenuDeleteLoopDone startMenuDeleteLoop
;  startMenuDeleteLoopDone:

  DeleteRegKey /ifempty HKCU "Software\OpenBabel ${OBVERSION}"
  DeleteRegKey          HKCU "Software\OpenBabelGUI"
  DeleteRegKey          HKCR "*\shell\OBDepict"
  
  ; Remove entry in Add/Remove Programs
  DeleteRegKey HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\OpenBabel-${OBVERSION}" 

  ; Remove entry in Add/Remove Programs
  DeleteRegKey HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\OpenBabel-${OBVERSION}" 

  ; Remove env var
  push "BABEL_DATADIR"
  push $INSTDIR
  Call un.RemoveFromEnvVar
  DeleteRegValue        HKCU "Environment" "BABEL_DATADIR"

  ClearErrors
  ReadEnvStr $0 APPDATA
  IfErrors 0 +3    
  StrCpy $DataBase "$INSTDIR"
  Goto +2
  StrCpy $DataBase "$%APPDATA%\OpenBabel-${OBVERSION}"
  RMDir /r "$DataBase\data"
  RMDir /r "$DataBase\examples"
  RMDir /r "$DataBase\doc"
  ;RMDir /r "$DataBase" is not safe 
  RMDir "$DataBase"
  
SectionEnd
