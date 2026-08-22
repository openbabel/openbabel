; Open Babel Windows installer
; Derived from the legacy maintenance installer while keeping the
; current native x64 payload and dependency layout.

!include "MUI.nsh"

!ifndef SourceDir
!define SourceDir "."
!endif
!ifndef BuildDir
!define BuildDir "."
!endif
!ifndef DepsDir
!define DepsDir "."
!endif
!ifndef VCRedist
!define VCRedist "vc_redist.x64.exe"
!endif
!ifndef myOutFile
!define myOutFile "OpenBabel-Installer.exe"
!endif
!ifndef OBVersion
!define OBVersion "0.0.0"
!endif

Name "OpenBabel ${OBVersion}"
Caption "OpenBabel ${OBVersion} Setup"
OutFile "${myOutFile}"
InstallDir "$PROGRAMFILES64\OpenBabel-${OBVersion}"
InstallDirRegKey HKCU "Software\OpenBabel ${OBVERSION}" ""
RequestExecutionLevel admin

Var STARTMENU_FOLDER

!ifndef WriteEnvStr_RegKey
!ifdef ALL_USERS
!define WriteEnvStr_RegKey 'HKLM "SYSTEM\CurrentControlSet\Control\Session Manager\Environment"'
!else
!define WriteEnvStr_RegKey 'HKCU "Environment"'
!endif
!endif

!macro IsNT un
Function ${un}IsNT
  Push $0
  ReadRegStr $0 HKLM "SOFTWARE\Microsoft\Windows NT\CurrentVersion" CurrentVersion
  StrCmp $0 "" 0 IsNT_yes
  Pop $0
  Push 0
  Return
IsNT_yes:
  Pop $0
  Push 1
FunctionEnd
!macroend
!insertmacro IsNT ""
!insertmacro IsNT "un."

!macro StrStr un
Function ${un}StrStr
  Exch $R1
  Exch
  Exch $R2
  Push $R3
  Push $R4
  Push $R5
  StrLen $R3 $R1
  StrCpy $R4 0
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

Function Trim
  Exch $R1
  Push $R2
Loop:
  StrCpy $R2 "$R1" 1 -1
  StrCmp "$R2" " " RTrim
  StrCmp "$R2" "$\n" RTrim
  StrCmp "$R2" "$\r" RTrim
  StrCmp "$R2" ";" RTrim
  Goto Done
RTrim:
  StrCpy $R1 "$R1" -1
  Goto Loop
Done:
  Pop $R2
  Exch $R1
FunctionEnd

Function AddToPath
  Exch $0
  Push $1
  Push $2
  Push $3
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
    StrCpy $1 $WINDIR 2
    FileOpen $1 "$1\autoexec.bat" a
    FileSeek $1 -1 END
    FileReadByte $1 $2
    IntCmp $2 26 0 +2 +2
      FileSeek $1 -1 END
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

Function un.RemoveFromPath
  Exch $0
  Push $1
  Push $2
  Push $3
  Push $4
  Push $5
  Push $6
  Call un.IsNT
  Pop $1
  StrCmp $1 1 unRemoveFromPath_NT
    StrCpy $1 $WINDIR 2
    FileOpen $1 "$1\autoexec.bat" r
    GetTempFileName $4
    FileOpen $2 $4 w
    GetFullPathName /SHORT $0 $0
    StrCpy $0 "SET PATH=%PATH%;$0"
    Goto unRemoveFromPath_dosLoop
unRemoveFromPath_dosLoop:
  FileRead $1 $3
  StrCpy $5 $3 1 -1
  StrCmp $5 $6 0 +2
    StrCpy $3 $3 -1
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
  StrCpy $5 $1 1 -1
  StrCmp $5 ";" +2
    StrCpy $1 "$1;"
  Push $1
  Push "$0;"
  Call un.StrStr
  Pop $2
  StrCmp $2 "" unRemoveFromPath_done
  StrLen $3 "$0;"
  StrLen $4 $2
  StrCpy $5 $1 -$4
  StrCpy $6 $2 "" $3
  StrCpy $3 $5$6
  StrCpy $5 $3 1 -1
  StrCmp $5 ";" 0 +2
    StrCpy $3 $3 -1
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

!define MUI_ABORTWARNING
!define MUI_FINISHPAGE_RUN "$INSTDIR\obgui.exe"
!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_LICENSE "${SourceDir}\COPYING"
!insertmacro MUI_PAGE_DIRECTORY
!define MUI_STARTMENUPAGE_REGISTRY_ROOT "HKCU"
!define MUI_STARTMENUPAGE_REGISTRY_KEY "Software\OpenBabel ${OBVERSION}"
!define MUI_STARTMENUPAGE_REGISTRY_VALUENAME "Start Menu Folder"
!insertmacro MUI_PAGE_STARTMENU Application $STARTMENU_FOLDER
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
!insertmacro MUI_LANGUAGE "English"

Section "Open Babel" SecOpenBabel
  SectionIn RO
  SetOutPath "$INSTDIR"
  File /r "${BuildDir}\bin\Release\*.*"
  File /nonfatal "${DepsDir}\libs-common\x64\*.dll"
  File /nonfatal "${DepsDir}\libs-vs12\x64\*.dll"
  File "${VCRedist}"
  WriteRegStr HKCU "Software\OpenBabel ${OBVERSION}" "" "$INSTDIR"
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\OpenBabel-${OBVERSION}" "DisplayName" "OpenBabel-${OBVERSION}"
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\OpenBabel-${OBVERSION}" "UninstallString" '"$INSTDIR\Uninstall.exe"'
  Push $INSTDIR
  Call AddToPath
  WriteUninstaller "$INSTDIR\Uninstall.exe"
  !insertmacro MUI_STARTMENU_WRITE_BEGIN Application
    CreateDirectory "$SMPROGRAMS\$STARTMENU_FOLDER"
    CreateShortCut "$SMPROGRAMS\$STARTMENU_FOLDER\Open Babel GUI.lnk" "$INSTDIR\obgui.exe"
    CreateShortCut "$SMPROGRAMS\$STARTMENU_FOLDER\Uninstall.lnk" "$INSTDIR\Uninstall.exe"
  !insertmacro MUI_STARTMENU_WRITE_END
  ExecWait '"$INSTDIR\vc_redist.x64.exe" /quiet'
  Delete "$INSTDIR\vc_redist.x64.exe"
SectionEnd

Section "Uninstall"
  Push $INSTDIR
  Call un.RemoveFromPath
  Delete "$SMPROGRAMS\$STARTMENU_FOLDER\Open Babel GUI.lnk"
  Delete "$SMPROGRAMS\$STARTMENU_FOLDER\Uninstall.lnk"
  RMDir "$SMPROGRAMS\$STARTMENU_FOLDER"
  DeleteRegKey HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\OpenBabel-${OBVERSION}"
  DeleteRegKey HKCU "Software\OpenBabel ${OBVERSION}"
  RMDir /r "$INSTDIR"
SectionEnd
