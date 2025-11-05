#!/bin/bash
# Mutalyze Executable Builder for Linux/Mac
# Builds a standalone executable with full OpenMM support

set -e  # Exit on error

echo "========================================"
echo " MUTALYZE EXECUTABLE BUILDER"
echo "========================================"
echo ""

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ ERROR: Conda not found!"
    echo "Please install Miniconda from: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Check if environment exists
if ! conda env list | grep -q "mutalyze_build"; then
    echo "Creating conda environment..."
    conda create -n mutalyze_build python=3.10 -y
fi

echo ""
echo "Activating environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mutalyze_build

echo ""
echo "Installing OpenMM..."
conda install -c conda-forge openmm -y

echo ""
echo "Installing Python dependencies..."
pip install -r requirements_exe.txt

echo ""
echo "Verifying installations..."
python -c "import openmm; print(f'✅ OpenMM {openmm.version.version} OK')"
python -c "import PyInstaller; print(f'✅ PyInstaller {PyInstaller.__version__} OK')"

echo ""
echo "========================================"
echo " BUILDING EXECUTABLE"
echo "========================================"
echo ""
echo "This may take 5-10 minutes..."
echo ""

# Build using spec file
pyinstaller build_exe_advanced.spec

echo ""
if [ -f "dist/Mutalyze" ]; then
    echo "========================================"
    echo " ✅ BUILD SUCCESSFUL!"
    echo "========================================"
    echo ""
    echo "Executable location: dist/Mutalyze"
    echo ""
    echo "File size:"
    ls -lh dist/Mutalyze | awk '{print $5, $9}'
    echo ""
    echo "To test, run: ./dist/Mutalyze"
    echo ""
else
    echo "========================================"
    echo " ❌ BUILD FAILED!"
    echo "========================================"
    echo ""
    echo "Check the output above for errors."
    echo "Common issues:"
    echo "  - OpenMM not installed"
    echo "  - Missing dependencies"
    echo "  - File path issues"
    echo ""
    exit 1
fi
