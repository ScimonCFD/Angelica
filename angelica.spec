# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_data_files

# sv_ttk ships TCL theme files and PNG sprites that must travel with the exe
sv_ttk_datas = collect_data_files("sv_ttk")

a = Analysis(
    ["launcher.py"],
    pathex=["src"],
    binaries=[],
    datas=sv_ttk_datas,
    hiddenimports=[
        "sv_ttk",
        "openpyxl",
        "openpyxl.cell._writer",
        "openpyxl.styles.builtins",
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    # matplotlib is only used by the optional CLI reporting function, not the GUI
    excludes=["matplotlib", "scipy", "IPython", "pandas"],
    noarchive=False,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="Angelica",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,          # no black console window behind the GUI
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon="installer\\angelica.ico",
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name="Angelica",
)
