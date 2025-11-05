# -*- mode: python ; coding: utf-8 -*-
"""
Advanced PyInstaller Spec File for Mutalyze
Includes OpenMM, Streamlit, and all dependencies
"""

import sys
from pathlib import Path

# Paths
block_cipher = None
SCRIPT_DIR = Path('.').absolute()

# Collect all data files
added_files = [
    ('models', 'models'),
    ('data', 'data'),
    ('forcefields', 'forcefields'),
    ('assets', 'assets'),
    ('core', 'core'),
    ('app.py', '.'),
]

# Hidden imports (critical for Streamlit + OpenMM)
hidden_imports = [
    # Streamlit core
    'streamlit',
    'streamlit.web.cli',
    'streamlit.web.server',
    'streamlit.runtime',
    'streamlit.runtime.scriptrunner',
    'streamlit.runtime.scriptrunner.magic_funcs',
    'streamlit.components.v1',
    'streamlit.elements',
    
    # OpenMM (CRITICAL)
    'openmm',
    'openmm.app',
    'openmm.unit',
    'openmm.openmm',
    
    # BioPython
    'Bio',
    'Bio.PDB',
    'Bio.PDB.PDBIO',
    'Bio.PDB.PDBParser',
    'Bio.PDB.Select',
    'Bio.Align',
    'Bio.Blast',
    'Bio.SeqIO',
    
    # ML
    'sklearn',
    'sklearn.ensemble',
    'sklearn.ensemble._forest',
    'sklearn.tree',
    'sklearn.tree._tree',
    
    # Visualization
    'py3Dmol',
    'matplotlib',
    'matplotlib.pyplot',
    'seaborn',
    
    # Scientific
    'numpy',
    'pandas',
    'scipy',
    'scipy.stats',
    
    # Others
    'tqdm',
    'requests',
    'pillow',
    'altair',
    'packaging',
]

# Analysis
a = Analysis(
    ['mutalyze_launcher.py'],
    pathex=[str(SCRIPT_DIR)],
    binaries=[],
    datas=added_files,
    hiddenimports=hidden_imports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['pytest', 'unittest'],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

# Remove duplicate entries
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

# Create executable
exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='Mutalyze',
    debug=False,
    bootloader_ignore_signals=False,
    strip=True,
    upx=False,  # Better compatibility
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,  # Show console for debugging (set False for release)
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon='assets/icon.ico' if Path('assets/icon.ico').exists() else None,
)
