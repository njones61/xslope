; XSLOPE Studio - Windows installer (NSIS 3.x)
;
; Per-user install, no administrator rights: everything lands under
; %LOCALAPPDATA%\Programs\XSLOPE Studio and every registry key is written to
; HKCU. That keeps SmartScreen/UAC out of the way for classroom machines and
; lets the in-app updater re-run this installer silently.
;
; Built by CI (see .github/workflows/release-installers.yml) from the PyInstaller
; onedir output:
;
;   makensis /DVERSION=0.1.58 ^
;            /DSRCDIR="..\..\packaging\dist\XSLOPE Studio" ^
;            /DOUTFILE="..\..\packaging\dist\XSLOPE-Studio-0.1.58-windows-x64-setup.exe" ^
;            installer.nsi
;
; Silent upgrade in place (what the updater runs):
;
;   "XSLOPE-Studio-0.1.58-windows-x64-setup.exe" /S
;
; /S suppresses all UI, reuses the recorded install directory, closes a running
; copy of the app, replaces the payload and rewrites the uninstall entry.

Unicode true

!include "MUI2.nsh"
!include "LogicLib.nsh"
!include "FileFunc.nsh"   ; ${GetSize}, for the Add/Remove Programs size

;--------------------------------------------------------------------------
; Parameters (all overridable from the makensis command line)
;--------------------------------------------------------------------------
!ifndef VERSION
  !define VERSION "0.0.0"
!endif
!ifndef SRCDIR
  !define SRCDIR "${__FILEDIR__}\..\dist\XSLOPE Studio"
!endif
!ifndef OUTFILE
  !define OUTFILE "${__FILEDIR__}\..\dist\XSLOPE-Studio-${VERSION}-windows-x64-setup.exe"
!endif

!define APPNAME    "XSLOPE Studio"
!define SHORTNAME  "XSLOPEStudio"
!define EXENAME    "XSLOPE Studio.exe"
!define PUBLISHER  "Norman L. Jones"
!define HOMEPAGE   "https://github.com/njones61/xslope"
!define APPKEY     "Software\${SHORTNAME}"
!define UNINSTKEY  "Software\Microsoft\Windows\CurrentVersion\Uninstall\${SHORTNAME}"

Name "${APPNAME} ${VERSION}"
OutFile "${OUTFILE}"
InstallDir "$LOCALAPPDATA\Programs\${APPNAME}"
InstallDirRegKey HKCU "${APPKEY}" "InstallDir"
RequestExecutionLevel user      ; per-user: never prompts for admin
SetCompressor /SOLID lzma
ShowInstDetails show
ShowUninstDetails show

;--------------------------------------------------------------------------
; Version resource on the installer .exe itself, so Explorer, the release page
; and the updater all see the same number that Studio reports.
;--------------------------------------------------------------------------
; VIProductVersion insists on x.x.x.x, so split "0.1.58" into its parts. The
; !ifndef fallbacks keep a non-numeric version (e.g. "0.1.58rc1") from failing
; the compile - the human-readable strings below still carry the exact version.
!searchparse /noerrors "${VERSION}" "" VER_MAJOR "." VER_MINOR "." VER_PATCH
!ifndef VER_MAJOR
  !define VER_MAJOR 0
!endif
!ifndef VER_MINOR
  !define VER_MINOR 0
!endif
!ifndef VER_PATCH
  !define VER_PATCH 0
!endif
VIProductVersion "${VER_MAJOR}.${VER_MINOR}.${VER_PATCH}.0"
VIAddVersionKey "ProductName"     "${APPNAME}"
VIAddVersionKey "ProductVersion"  "${VERSION}"
VIAddVersionKey "FileVersion"     "${VERSION}"
VIAddVersionKey "FileDescription" "${APPNAME} installer"
VIAddVersionKey "CompanyName"     "${PUBLISHER}"
VIAddVersionKey "LegalCopyright"  "Apache-2.0"

;--------------------------------------------------------------------------
; UI (skipped entirely under /S)
;--------------------------------------------------------------------------
!define MUI_ABORTWARNING
; build_app.py renders packaging\build\xslope.ico from the app PNG when Pillow is
; installed. Without it the installer simply uses the default NSIS icon.
!define ICONFILE "${__FILEDIR__}\..\build\xslope.ico"
!if /FileExists "${ICONFILE}"
  !define MUI_ICON "${ICONFILE}"
  !define MUI_UNICON "${ICONFILE}"
!endif
!define MUI_FINISHPAGE_RUN "$INSTDIR\${EXENAME}"
!define MUI_FINISHPAGE_RUN_TEXT "Launch ${APPNAME}"

!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
!insertmacro MUI_LANGUAGE "English"

;--------------------------------------------------------------------------
; Helpers
;--------------------------------------------------------------------------
; Close a running copy so an upgrade can overwrite the payload. Only this
; user's processes are touched, which needs no elevation.
!macro CloseRunningApp
  DetailPrint "Closing any running ${APPNAME}..."
  nsExec::Exec 'taskkill /IM "${EXENAME}" /FI "USERNAME eq $%USERNAME%"'
  Pop $0
  Sleep 500
  nsExec::Exec 'taskkill /F /IM "${EXENAME}" /FI "USERNAME eq $%USERNAME%"'
  Pop $0
  Sleep 500
!macroend

;--------------------------------------------------------------------------
; Install
;--------------------------------------------------------------------------
Section "${APPNAME}" SecMain
  SectionIn RO

  !insertmacro CloseRunningApp

  ; Upgrade in place: wipe the old payload directory first, otherwise stale
  ; modules from a previous version stay behind in _internal and can shadow the
  ; new ones. User data never lives here, so this is safe.
  RMDir /r "$INSTDIR\_internal"
  Delete "$INSTDIR\${EXENAME}"

  SetOutPath "$INSTDIR"
  File /r "${SRCDIR}\*"

  WriteUninstaller "$INSTDIR\Uninstall.exe"

  ; Start-menu shortcut (per-user).
  CreateDirectory "$SMPROGRAMS\${APPNAME}"
  CreateShortCut "$SMPROGRAMS\${APPNAME}\${APPNAME}.lnk" "$INSTDIR\${EXENAME}" \
    "" "$INSTDIR\${EXENAME}" 0
  CreateShortCut "$SMPROGRAMS\${APPNAME}\Uninstall ${APPNAME}.lnk" \
    "$INSTDIR\Uninstall.exe"

  ; What the app itself reads to know where and what it is.
  WriteRegStr HKCU "${APPKEY}" "InstallDir" "$INSTDIR"
  WriteRegStr HKCU "${APPKEY}" "Version" "${VERSION}"

  ; Add/Remove Programs entry (per-user hive: no admin needed).
  WriteRegStr HKCU "${UNINSTKEY}" "DisplayName"     "${APPNAME}"
  WriteRegStr HKCU "${UNINSTKEY}" "DisplayVersion"  "${VERSION}"
  WriteRegStr HKCU "${UNINSTKEY}" "Publisher"       "${PUBLISHER}"
  WriteRegStr HKCU "${UNINSTKEY}" "URLInfoAbout"    "${HOMEPAGE}"
  WriteRegStr HKCU "${UNINSTKEY}" "DisplayIcon"     "$INSTDIR\${EXENAME}"
  WriteRegStr HKCU "${UNINSTKEY}" "InstallLocation" "$INSTDIR"
  WriteRegStr HKCU "${UNINSTKEY}" "UninstallString" '"$INSTDIR\Uninstall.exe"'
  WriteRegStr HKCU "${UNINSTKEY}" "QuietUninstallString" '"$INSTDIR\Uninstall.exe" /S'
  WriteRegDWORD HKCU "${UNINSTKEY}" "NoModify" 1
  WriteRegDWORD HKCU "${UNINSTKEY}" "NoRepair" 1

  ; Size shown in Add/Remove Programs (KB).
  ${GetSize} "$INSTDIR" "/S=0K" $0 $1 $2
  IntFmt $0 "0x%08X" $0
  WriteRegDWORD HKCU "${UNINSTKEY}" "EstimatedSize" "$0"
SectionEnd

;--------------------------------------------------------------------------
; Uninstall
;--------------------------------------------------------------------------
Section "Uninstall"
  !insertmacro CloseRunningApp

  Delete "$SMPROGRAMS\${APPNAME}\${APPNAME}.lnk"
  Delete "$SMPROGRAMS\${APPNAME}\Uninstall ${APPNAME}.lnk"
  RMDir "$SMPROGRAMS\${APPNAME}"

  RMDir /r "$INSTDIR\_internal"
  Delete "$INSTDIR\${EXENAME}"
  Delete "$INSTDIR\Uninstall.exe"
  RMDir /r "$INSTDIR"

  DeleteRegKey HKCU "${UNINSTKEY}"
  DeleteRegKey HKCU "${APPKEY}"
SectionEnd
