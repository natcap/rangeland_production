; Variables needed at the command line:
; VERSION         - the version of InVEST we're building (example: 3.4.5)
;                   This string must not contain characters that are
;                   problematic in Windows paths (no : , etc)
; BINDIR          - The local folder of binaries to include.
; ARCHITECTURE    - The architecture we're building for.  Generally this is x86.
; FORKNAME        - The username of the InVEST fork we're building off of.

!include nsProcess.nsh
!include LogicLib.nsh
; HM NIS Edit Wizard helper defines
!define PRODUCT_NAME "RPM"
!define PRODUCT_VERSION "${VERSION} ${ARCHITECTURE}"
!define PRODUCT_PUBLISHER "The Natural Capital Project"
!define PRODUCT_WEB_SITE "https://www.naturalcapitalproject.org"
!define MUI_COMPONENTSPAGE_NODESC
!define PACKAGE_NAME "${PRODUCT_NAME} ${PRODUCT_VERSION}"

SetCompressor zlib
!define MUI_WELCOMEFINISHPAGE_BITMAP "rpm-vertical.bmp"
!define MUI_UNWELCOMEFINISHPAGE_BITMAP "rpm-vertical.bmp"
!define MUI_HEADERIMAGE
!define MUI_HEADERIMAGE_BITMAP "rpm-horizontal.bmp"
!define MUI_UNHEADERIMAGE_BITMAP "rpm-horizontal.bmp"
!define MUI_UNICON "${NSISDIR}\Contrib\Graphics\Icons\orange-uninstall.ico"

; MUI 1.67 compatible ------
!include "MUI2.nsh"
!include "LogicLib.nsh"
!include "x64.nsh"
!include "FileFunc.nsh"
!include "nsDialogs.nsh"
!include "WinVer.nsh"

; MUI Settings
!define MUI_ABORTWARNING
!define MUI_ICON "rpm.ico"

; Add an advanced options control for the welcome page.
!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_LICENSE "..\..\LICENSE.txt"
!insertmacro MUI_PAGE_COMPONENTS
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH

; MUI Uninstaller settings---------------
!insertmacro MUI_UNPAGE_WELCOME
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
!insertmacro MUI_UNPAGE_FINISH

; Language files
!insertmacro MUI_LANGUAGE "English"

; MUI end ------

Name "${PRODUCT_NAME} ${PRODUCT_VERSION}"
OutFile ..\..\dist\Rangeland_Production_${FORKNAME}${VERSION}_${ARCHITECTURE}_Setup.exe
InstallDir "C:\RPM_${VERSION}_${ARCHITECTURE}"
ShowInstDetails show
RequestExecutionLevel admin

; This function allows us to test to see if a process is currently running.
; If the process name passed in is actually found, a message box is presented
; and the uninstaller should quit.
!macro CheckProgramRunning process_name
    ${nsProcess::FindProcess} "${process_name}.exe" $R0
    Pop $R0

    StrCmp $R0 603 +3
        MessageBox MB_OK|MB_ICONEXCLAMATION "RPM is still running.  Please close RPM and try again."
        Abort
!macroend

; NSIS 3.x defines these variables, NSIS 2.x does not.  This supports both versions.
!ifndef LVM_GETITEMCOUNT
    !define LVM_GETITEMCOUNT 0x1004
!endif
!ifndef LVM_GETITEMTEXT
    !define LVM_GETITEMTEXT 0x102D
!endif

Function DumpLog
    Exch $5
    Push $0
    Push $1
    Push $2
    Push $3
    Push $4
    Push $6

    FindWindow $0 "#32770" "" $HWNDPARENT
    GetDlgItem $0 $0 1016
    StrCmp $0 0 exit
    FileOpen $5 $5 "w"
    StrCmp $5 "" exit
        SendMessage $0 ${LVM_GETITEMCOUNT} 0 0 $6
        System::Alloc ${NSIS_MAX_STRLEN}
        Pop $3
        StrCpy $2 0
        System::Call "*(i, i, i, i, i, i, i, i, i) i \
            (0, 0, 0, 0, 0, r3, ${NSIS_MAX_STRLEN}) .r1"
        loop: StrCmp $2 $6 done
            System::Call "User32::SendMessageA(i, i, i, i) i \
            ($0, ${LVM_GETITEMTEXT}, $2, r1)"
            System::Call "*$3(&t${NSIS_MAX_STRLEN} .r4)"
            FileWrite $5 "$4$\r$\n"
            IntOp $2 $2 + 1
            Goto loop
        done:
            FileClose $5
            System::Free $1
            System::Free $3
    exit:
        Pop $6
        Pop $4
        Pop $3
        Pop $2
        Pop $1
        Pop $0
        Exch $5
FunctionEnd

Function Un.onInit
    !insertmacro CheckProgramRunning "rangeland_production"
FunctionEnd

; Copied into the invest folder later in the NSIS script
!define RPM_BINARIES "$INSTDIR\RPM-x86"
!define RPM_ICON "${RPM_BINARIES}\rpm.ico"

Section "RPM" Section_RPM_Tool
    AddSize 230793  ; This size is based on Windows build of InVEST 3.4.0
    SetShellVarContext all
    SectionIn RO ;require this section

    ; Write the uninstaller to disk
    SetOutPath "$INSTDIR"
    !define UNINSTALL_PATH "$INSTDIR\Uninstall_${VERSION}.exe"
    writeUninstaller "${UNINSTALL_PATH}"

    ; Create start  menu shortcuts.
    ; These shortcut paths are set in the appropriate places based on the SetShellVarConext flag.
    ; This flag is automatically set based on the MULTIUSER installation mode selected by the user.
    !define SMPATH "$SMPROGRAMS\${PACKAGE_NAME}"
    CreateDirectory "${SMPATH}"
    CreateShortCut "${SMPATH}\Rangeland Production ${VERSION}.lnk" "${RPM_BINARIES}\rangeland_production.exe" "" "${RPM_ICON}"

    ; Write registry keys for convenient uninstallation via add/remove programs.
    ; Inspired by the example at
    ; nsis.sourceforge.net/A_simple_installer_with_start_menu_shortcut_and_uninstaller
    !define REGISTRY_PATH "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_PUBLISHER} ${PRODUCT_NAME} ${PRODUCT_VERSION}"
    WriteRegStr HKLM "${REGISTRY_PATH}" "DisplayName"          "${PRODUCT_NAME} ${PRODUCT_VERSION}"
    WriteRegStr HKLM "${REGISTRY_PATH}" "UninstallString"      "${UNINSTALL_PATH}"
    WriteRegStr HKLM "${REGISTRY_PATH}" "QuietUninstallString" "${UNINSTALL_PATH} /S"
    WriteRegStr HKLM "${REGISTRY_PATH}" "InstallLocation"      "$INSTDIR"
    WriteRegStr HKLM "${REGISTRY_PATH}" "DisplayIcon"          "${RPM_ICON}"
    WriteRegStr HKLM "${REGISTRY_PATH}" "Publisher"            "${PRODUCT_PUBLISHER}"
    WriteRegStr HKLM "${REGISTRY_PATH}" "URLInfoAbout"         "${PRODUCT_WEB_SITE}"
    WriteRegStr HKLM "${REGISTRY_PATH}" "DisplayVersion"       "${PRODUCT_VERSION}"
    WriteRegDWORD HKLM "${REGISTRY_PATH}" "NoModify" 1
    WriteRegDWORD HKLM "${REGISTRY_PATH}" "NoRepair" 1


    ; Actually install the information we want to disk.
    SetOutPath "$INSTDIR"
    File ..\..\LICENSE.txt
    file ..\..\HISTORY.rst

    SetOutPath "${RPM_BINARIES}"
    File /r /x *.hg* /x *.svn* ..\..\${BINDIR}\*
    File rpm.ico

    ; Write the install log to a text file on disk.
    StrCpy $0 "$INSTDIR\install_log.txt"
    Push $0
    Call DumpLog

SectionEnd

; Only add this section if we're running the installer on Windows 7 or below.
; See InVEST Issue #3515 (https://bitbucket.org/natcap/invest/issues/3515)
; This section is disabled in .onInit if we're running Windows 8 or later.
Section "MSVCRT 2008 Runtime (Recommended)" Sec_VCRedist2008
    SetOutPath "$INSTDIR"
    File ..\..\build\vcredist_x86.exe
    ExecWait "vcredist_x86.exe /q"
SectionEnd

Section "uninstall"
  ; Need to enforce execution level as admin.  See
  ; nsis.sourceforge.net/Shortcuts_removal_fails_on_Windows_Vista
  SetShellVarContext all
  rmdir /r "$SMPROGRAMS\${PACKAGE_NAME}"

  ; Delete the installation directory on disk
  rmdir /r "$INSTDIR"

  ; Delete the entire registry key for this version of RIOS.
  DeleteRegKey HKLM "${REGISTRY_PATH}"
SectionEnd

Function .onInit
 ${GetOptions} $CMDLINE "/?" $0
 IfErrors skiphelp showhelp
 showhelp:
     MessageBox MB_OK "Rangeland Production Model$\r$\n\
     $\r$\n\
     For more information about RPM or the Natural Capital Project, visit our \
     website: http://naturalcapitalproject.org/invest$\r$\n\
     $\r$\n\
     Command-Line Options:$\r$\n\
         /?$\t$\t=$\tDisplay this help and exit$\r$\n\
         /S$\t$\t=$\tSilently install RPM.$\r$\n\
         "
     abort
 skiphelp:

 System::Call 'kernel32::CreateMutexA(i 0, i 0, t "RPM ${VERSION}") i .r1 ?e'
 Pop $R0

 StrCmp $R0 0 +3
   MessageBox MB_OK|MB_ICONEXCLAMATION "An RPM ${VERSION} installer is already running."
   Abort

  ${ifNot} ${AtMostWin7}
    ; disable the section if we're not running on Windows 7 or earlier.
    ; This section should not execute for Windows 8 or later.
    SectionGetFlags ${Sec_VCRedist2008} $0
    IntOp $0 $0 & ${SECTION_OFF}
    SectionSetFlags ${Sec_VCRedist2008} $0
    SectionSetText ${Sec_VCRedist2008} ""
  ${endIf}
FunctionEnd
