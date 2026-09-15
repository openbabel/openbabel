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

  ; Update only the current user's PATH. EnVar handles an unset/empty PATH
  ; without copying the process's merged system PATH into HKCU.
  EnVar::SetHKCU
  Pop $0
  EnVar::AddValue "PATH" "$INSTDIR"
  Pop $0
  StrCmp $0 "0" path_added
    MessageBox MB_OK|MB_ICONSTOP "Failed to add Open Babel to the user PATH. The installation will be aborted. (EnVar error: $0)"
    Abort
path_added:

  WriteUninstaller "$INSTDIR\Uninstall.exe"
  !insertmacro MUI_STARTMENU_WRITE_BEGIN Application
    CreateDirectory "$SMPROGRAMS\$STARTMENU_FOLDER"
    CreateShortCut "$SMPROGRAMS\$STARTMENU_FOLDER\Open Babel GUI.lnk" "$INSTDIR\obgui.exe"
    CreateShortCut "$SMPROGRAMS\$STARTMENU_FOLDER\Uninstall.lnk" "$INSTDIR\Uninstall.exe"
  !insertmacro MUI_STARTMENU_WRITE_END
  ExecWait '"$INSTDIR\vc_redist.x64.exe" /quiet' $0
  Delete "$INSTDIR\vc_redist.x64.exe"
  
  ${If} $0 != 0
      Abort "Visual C++ Redistributable installation failed with exit code $0."
  ${EndIf}
SectionEnd

Section "Uninstall"
  ; Remove only this installation directory from the current user's PATH.
  EnVar::SetHKCU
  Pop $0
  EnVar::DeleteValue "PATH" "$INSTDIR"
  Pop $0

  !insertmacro MUI_STARTMENU_GETFOLDER Application $STARTMENU_FOLDER
  Delete "$SMPROGRAMS\$STARTMENU_FOLDER\Open Babel GUI.lnk"
  Delete "$SMPROGRAMS\$STARTMENU_FOLDER\Uninstall.lnk"
  RMDir "$SMPROGRAMS\$STARTMENU_FOLDER"
  DeleteRegKey HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\OpenBabel-${OBVERSION}"
  DeleteRegKey HKCU "Software\OpenBabel ${OBVERSION}"
  RMDir /r "$INSTDIR"
SectionEnd
