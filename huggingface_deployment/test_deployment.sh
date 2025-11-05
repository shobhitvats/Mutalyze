#!/bin/bash
# Quick test script for HuggingFace deployment
# Run this locally before pushing to HuggingFace

echo "=================================================="
echo " Testing HuggingFace Deployment Locally"
echo "=================================================="
echo ""

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ ERROR: Conda not found!"
    echo "Please install Miniconda: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo "✅ Conda found"
echo ""

# Check if environment.yml exists
if [ ! -f "environment.yml" ]; then
    echo "❌ ERROR: environment.yml not found!"
    echo "Make sure you're in the huggingface_deployment folder"
    exit 1
fi

echo "✅ environment.yml found"
echo ""

# Create/update conda environment
echo "Creating conda environment from environment.yml..."
echo "(This may take 5-10 minutes on first run)"
echo ""

conda env create -f environment.yml -n mutalyze_hf_test 2>&1 | grep -v "^Collecting package metadata" || \
conda env update -f environment.yml -n mutalyze_hf_test

# Activate environment
echo ""
echo "Activating environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mutalyze_hf_test

# Test imports
echo ""
echo "Testing critical imports..."
echo ""

python -c "
import sys
success = True

tests = {
    'streamlit': 'Streamlit',
    'openmm': 'OpenMM (CRITICAL)',
    'openmm.app': 'OpenMM app',
    'Bio.PDB': 'BioPython',
    'sklearn': 'scikit-learn',
    'numpy': 'NumPy',
    'pandas': 'Pandas',
}

for module, name in tests.items():
    try:
        __import__(module)
        print(f'  ✅ {name}')
    except ImportError as e:
        print(f'  ❌ {name} - {e}')
        success = False

if not success:
    print()
    print('❌ Some imports failed!')
    sys.exit(1)
else:
    print()
    print('✅ All critical imports successful!')
" || exit 1

# Test OpenMM specifically
echo ""
echo "Testing OpenMM functionality..."
python -c "
import openmm
from openmm import app, unit

print(f'  ✅ OpenMM version: {openmm.version.version}')

# Test force field loading
try:
    ff = app.ForceField('amber14-all.xml', 'implicit/gbn2.xml')
    print('  ✅ AMBER force field loaded')
except Exception as e:
    print(f'  ❌ Force field failed: {e}')
    exit(1)
" || exit 1

# Test core modules
echo ""
echo "Testing core modules..."
python -c "
import sys
sys.path.insert(0, 'core')

modules = [
    'pdb_utils',
    'mutation_builder', 
    'energy_calc',
    'simple_energy',
]

for mod in modules:
    try:
        __import__(mod)
        print(f'  ✅ {mod}')
    except Exception as e:
        print(f'  ❌ {mod}: {e}')
        exit(1)

# Test OPENMM_AVAILABLE flag
from energy_calc import OPENMM_AVAILABLE
if OPENMM_AVAILABLE:
    print('  ✅ OPENMM_AVAILABLE = True')
else:
    print('  ❌ OPENMM_AVAILABLE = False (should be True!)')
    exit(1)
" || exit 1

# Test model file
echo ""
echo "Testing ML model..."
python -c "
import pickle
from pathlib import Path

model_path = Path('models/calibration_v5.pkl')
if not model_path.exists():
    print(f'  ❌ Model not found: {model_path}')
    exit(1)

try:
    with open(model_path, 'rb') as f:
        model = pickle.load(f)
    print(f'  ✅ Model loaded: {type(model).__name__}')
except Exception as e:
    print(f'  ❌ Model load failed: {e}')
    exit(1)
" || exit 1

echo ""
echo "=================================================="
echo " ✅ ALL TESTS PASSED!"
echo "=================================================="
echo ""
echo "Your HuggingFace deployment is ready to upload!"
echo ""
echo "Next steps:"
echo "  1. Create a Space on HuggingFace"
echo "  2. Upload all files from this folder"
echo "  3. HuggingFace will use environment.yml automatically"
echo "  4. Wait 10-15 minutes for first build"
echo "  5. Your app will be live!"
echo ""
echo "Optional: Run Streamlit locally to test:"
echo "  streamlit run app.py"
echo ""
