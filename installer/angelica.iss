#define MyAppName      "Angelica"
#define MyAppVersion   "0.1.0"
#define MyAppPublisher "Simon Rodriguez"
#define MyAppURL       "https://scimonCFD.github.io/Angelica_website/"
#define MyAppExeName   "Angelica.exe"

[Setup]
AppId={{5A8B7C3D-E2F1-4096-B8A3-9D7C6E5F4B2A}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
AppPublisherURL={#MyAppURL}
AppSupportURL={#MyAppURL}
AppUpdatesURL={#MyAppURL}
DefaultDirName={autopf}\{#MyAppName}
DefaultGroupName={#MyAppName}
AllowNoIcons=yes
OutputDir=Output
OutputBaseFilename=AngelicaSetup-{#MyAppVersion}
Compression=lzma
SolidCompression=yes
WizardStyle=modern
SetupIconFile=angelica.ico
WizardImageFile=wizard_left.png
WizardSmallImageFile=wizard_small.png
; Install per-user by default so no admin rights are required
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
ArchitecturesAllowed=x64
ArchitecturesInstallIn64BitMode=x64

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[Files]
Source: "..\dist\Angelica\*";  DestDir: "{app}";           Flags: ignoreversion recursesubdirs createallsubdirs
Source: "..\tutorials\*";      DestDir: "{app}\tutorials"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\{#MyAppName}";                                          Filename: "{app}\{#MyAppExeName}"; IconFilename: "{app}\{#MyAppExeName}"
Name: "{group}\{cm:UninstallProgram,{#MyAppName}}";                    Filename: "{uninstallexe}"
Name: "{autodesktop}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; IconFilename: "{app}\{#MyAppExeName}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#MyAppExeName}"; Description: "{cm:LaunchProgram,{#StringChange(MyAppName,'&','&&')}}"; Flags: nowait postinstall skipifsilent
