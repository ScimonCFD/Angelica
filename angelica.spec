# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_data_files, collect_all

# sv_ttk ships TCL theme files and PNG sprites that must travel with the exe
sv_ttk_datas = collect_data_files("sv_ttk")

# thermo + its dependencies (chemicals, fluids) use lazy imports and data files
# that PyInstaller cannot discover through static analysis alone.
thermo_datas,    thermo_bins,    thermo_imports    = collect_all("thermo")
chemicals_datas, chemicals_bins, chemicals_imports = collect_all("chemicals")
fluids_datas,    fluids_bins,    fluids_imports    = collect_all("fluids")

a = Analysis(
    ["launcher.py"],
    pathex=["src"],
    binaries=[*thermo_bins, *chemicals_bins, *fluids_bins],
    datas=[
        *sv_ttk_datas,
        *thermo_datas,
        *chemicals_datas,
        *fluids_datas,
        ("installer/angelica_32.png", "."),
    ],
    hiddenimports=[
        "sv_ttk",
        "openpyxl",
        "openpyxl.cell._writer",
        "openpyxl.styles.builtins",
        "scipy.sparse",
        "scipy.sparse.linalg",
        "scipy.sparse._csr",
        "scipy.sparse.linalg._dsolve",
        "scipy.sparse.linalg._dsolve.linsolve",
        *thermo_imports,
        *chemicals_imports,
        *fluids_imports,
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    # matplotlib is only used by the optional CLI reporting function, not the GUI
    excludes=["matplotlib", "IPython"],
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
