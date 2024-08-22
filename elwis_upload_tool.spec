# -*- mode: python ; coding: utf-8 -*-

block_cipher = None

a = Analysis(
    ['app/main.py'],
    pathex=[],
    binaries=[],  # Initialize binaries here
    datas=[
        ('burana.ico', '.'),  # Icon file
        ('app/styles.kv', '.'),   # Kivy style file
        ('config/*.json', 'config'),  # Configuration files
        ('README.md', '.'),   # Readme file
        ('app/libs/libopenbabel.so', '.'),  # Include the copied Open Babel library
    ],
    hiddenimports=[
        'rdkit',
        'openbabel',
        'requests',
        'dotenv',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='elwis_upload_tool',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    # The icon will be ignored on Linux, so this can be removed if desired
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='elwis_upload_tool',
)
