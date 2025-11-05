#!/bin/bash

# Mutalyze Setup Script
# Automates the installation and setup of Mutalyze

echo "🧬 Mutalyze Setup Script"
echo "========================================"
echo ""

# Check Python version
echo "Checking Python version..."
python_version=$(python3 --version 2>&1 | awk '{print $2}')
echo "Found Python $python_version"
echo ""

# Create virtual environment (optional)
read -p "Create a virtual environment? (recommended) [y/N]: " create_venv
if [[ $create_venv == "y" || $create_venv == "Y" ]]; then
    echo "Creating virtual environment..."
    python3 -m venv venv
    source venv/bin/activate
    echo "✓ Virtual environment created and activated"
    echo ""
fi

# Install Python dependencies
echo "Installing Python dependencies..."
pip install --upgrade pip
pip install -r requirements.txt
echo "✓ Python dependencies installed"
echo ""

# Check for conda
if command -v conda &> /dev/null; then
    echo "Conda detected!"
    read -p "Install OpenMM via conda? (recommended for energy calculations) [y/N]: " install_openmm
    
    if [[ $install_openmm == "y" || $install_openmm == "Y" ]]; then
        echo "Installing OpenMM..."
        conda install -c conda-forge openmm -y
        echo "✓ OpenMM installed"
    fi
    
    read -p "Install BLAST+ and MAFFT for conservation analysis? [y/N]: " install_tools
    
    if [[ $install_tools == "y" || $install_tools == "Y" ]]; then
        echo "Installing BLAST+ and MAFFT..."
        conda install -c bioconda blast mafft -y
        echo "✓ BLAST+ and MAFFT installed"
    fi
else
    echo "⚠️  Conda not detected. OpenMM requires conda for installation."
    echo "   Please install Miniconda or Anaconda: https://docs.conda.io/en/latest/miniconda.html"
fi

echo ""
echo "========================================"
echo "✓ Setup complete!"
echo ""
echo "To run Mutalyze:"
echo "  1. Activate virtual environment (if created): source venv/bin/activate"
echo "  2. Run the app: streamlit run app.py"
echo ""
echo "For documentation, see README.md and Instructions.txt"
echo "========================================"
