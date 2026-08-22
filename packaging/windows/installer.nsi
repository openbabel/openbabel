!include "MUI2.nsh"

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

Name "Open Babel"
Caption "Open Babel ${OBVersion} Setup"
OutFile "${myOutFile}"
InstallDir "$PROGRAMFILES64\OpenBabel"
InstallDirRegKey HKLM "Software\OpenBabel" "InstallDir"
RequestExecutionLevel admin
Unicode True

!define MUI_ABORTWARNING
!define MUI_ICON "${NSISDIR}\Contrib\Graphics\Icons\modern-install.ico"
!define MUI_UNICON "${NSISDIR}\Contrib\Graphics\Icons\modern-uninstall.ico"

!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH
!insertmacro MUI_LANGUAGE "English"

Section "Open Babel" SecOpenBabel
  SectionIn RO
  SetOutPath "$INSTDIR"

  ; Install the complete native Open Babel release tree. This deliberately
  ; replaces the old hand-maintained list of executables and plugins.
  File /r "${BuildDir}\bin\Release\*.*"

  ; Runtime DLLs supplied by the legacy MSVC dependency bundle and wxWidgets.
  File /nonfatal "${DepsDir}\libs-common\x64\*.dll"
  File /nonfatal "${DepsDir}\libs-vs12\x64\*.dll"

  ; Microsoft Visual C++ runtime.
  File "${VCRedist}"

  WriteRegStr HKLM "Software\OpenBabel" "InstallDir" "$INSTDIR"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Open Babel" "DisplayName" "Open Babel"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Open Babel" "DisplayVersion" "${OBVersion}"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Open Babel" "InstallLocation" "$INSTDIR"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Open Babel" "UninstallString" '"$INSTDIR\Uninstall.exe"'
  WriteUninstaller "$INSTDIR\Uninstall.exe"

  CreateDirectory "$SMPROGRAMS\Open Babel"
  CreateShortCut "$SMPROGRAMS\Open Babel\Open Babel GUI.lnk" "$INSTDIR\obgui.exe"
  CreateShortCut "$DESKTOP\Open Babel GUI.lnk" "$INSTDIR\obgui.exe"

  ExecWait '"$INSTDIR\vc_redist.x64.exe" /install /quiet /norestart'
  Delete "$INSTDIR\vc_redist.x64.exe"
SectionEnd

Section "Uninstall"
  Delete "$DESKTOP\Open Babel GUI.lnk"
  Delete "$SMPROGRAMS\Open Babel\Open Babel GUI.lnk"
  RMDir "$SMPROGRAMS\Open Babel"
  DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Open Babel"
  DeleteRegKey HKLM "Software\OpenBabel"
  RMDir /r "$INSTDIR"
SectionEnd
