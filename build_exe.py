"""
Mutalyze Executable Builder
Creates standalone .exe with OpenMM support for Windows
"""

import PyInstaller.__main__
import sys
import os
from pathlib import Path

# Get absolute paths
SCRIPT_DIR = Path(__file__).parent.absolute()
CORE_DIR = SCRIPT_DIR / 'core'
MODELS_DIR = SCRIPT_DIR / 'models'
DATA_DIR = SCRIPT_DIR / 'data'
FORCEFIELDS_DIR = SCRIPT_DIR / 'forcefields'
ASSETS_DIR = SCRIPT_DIR / 'assets'

# Build the executable
PyInstaller.__main__.run([
    # Main script
    str(SCRIPT_DIR / 'app.py'),
    
    # Output settings
    '--name=Mutalyze',
    '--onefile',  # Single executable file
    '--windowed',  # No console window (comment out for debugging)
    '--icon=assets/icon.ico' if (ASSETS_DIR / 'icon.ico').exists() else '--noconfirm',
    
    # Include data files
    f'--add-data={MODELS_DIR};models',
    f'--add-data={DATA_DIR};data',
    f'--add-data={FORCEFIELDS_DIR};forcefields',
    f'--add-data={ASSETS_DIR};assets',
    
    # Include core modules
    f'--add-data={CORE_DIR};core',
    
    # Hidden imports (modules not auto-detected)
    '--hidden-import=streamlit',
    '--hidden-import=streamlit.web.cli',
    '--hidden-import=streamlit.runtime.scriptrunner.magic_funcs',
    '--hidden-import=streamlit_option_menu',
    '--hidden-import=openmm',
    '--hidden-import=openmm.app',
    '--hidden-import=openmm.unit',
    '--hidden-import=sklearn',
    '--hidden-import=sklearn.ensemble',
    '--hidden-import=Bio',
    '--hidden-import=Bio.PDB',
    '--hidden-import=Bio.Align',
    '--hidden-import=Bio.Blast',
    '--hidden-import=py3Dmol',
    '--hidden-import=stmol',
    '--hidden-import=scipy',
    '--hidden-import=scipy.stats',
    '--hidden-import=matplotlib',
    '--hidden-import=seaborn',
    
    # Collect all submodules
    '--collect-all=streamlit',
    '--collect-all=openmm',
    '--collect-all=biopython',
    '--collect-all=sklearn',
    
    # Copy metadata
    '--copy-metadata=streamlit',
    '--copy-metadata=openmm',
    
    # Clean previous builds
    '--clean',
    
    # Confirm overwrites
    '--noconfirm',
    
    # Optimization
    '--strip',  # Strip symbols (smaller size)
    '--noupx',  # Don't use UPX (better compatibility)
])

print("\n" + "="*70)
print("✅ BUILD COMPLETE!")
print("="*70)
print(f"\n📦 Executable location: {SCRIPT_DIR / 'dist' / 'Mutalyze.exe'}")
print("\n📋 What's included:")
print("  ✅ Full OpenMM support (AMBER ff19SB + GBSA)")
print("  ✅ Random Forest v5 calibration model")
print("  ✅ All protein analysis features")
print("  ✅ Structure visualization")
print("  ✅ Batch processing")
print("  ✅ Conservation analysis")
print("  ✅ Interface analysis")
print("\n🚀 To run: Double-click Mutalyze.exe")
print("="*70 + "\n")
