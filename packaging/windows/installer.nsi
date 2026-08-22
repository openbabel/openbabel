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

; Keep the legacy page sequence: welcome, license, directory, start menu,
; installation, finish, then the standard uninstall pages.
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

  ; Native x64 Open Babel build output.
  File /r "${BuildDir}\bin\Release\*.*"

  ; Runtime dependencies. zlib1.dll comes from the current dependency
  ; snapshot and is intentionally not replaced by the legacy zlib DLL.
  File /nonfatal "${DepsDir}\libs-common\x64\*.dll"
  File /nonfatal "${DepsDir}\libs-vs12\x64\*.dll"

  ; Microsoft Visual C++ runtime.
  File "${VCRedist}"

  WriteRegStr HKCU "Software\OpenBabel ${OBVERSION}" "" "$INSTDIR"
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\OpenBabel-${OBVERSION}" "DisplayName" "OpenBabel-${OBVERSION}"
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\OpenBabel-${OBVERSION}" "UninstallString" '"$INSTDIR\Uninstall.exe"'

  ; Add the Open Babel installation directory to the machine PATH.
  ; WM_SETTINGCHANGE makes already-running applications aware that the
  ; environment has changed; newly started processes inherit it normally.
  EnVar::AddValue "PATH" "$INSTDIR"

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
  ; Remove only the value we added; preserve PATH entries belonging to other
  ; applications.
  EnVar::DeleteValue "PATH" "$INSTDIR"

  Delete "$SMPROGRAMS\$STARTMENU_FOLDER\Open Babel GUI.lnk"
  Delete "$SMPROGRAMS\$STARTMENU_FOLDER\Uninstall.lnk"
  RMDir "$SMPROGRAMS\$STARTMENU_FOLDER"
  DeleteRegKey HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\OpenBabel-${OBVERSION}"
  DeleteRegKey HKCU "Software\OpenBabel ${OBVERSION}"
  RMDir /r "$INSTDIR"
SectionEnd
