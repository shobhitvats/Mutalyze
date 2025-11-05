# 🏗️ Mutalyze Executable Build Instructions

Complete guide to building a standalone `.exe` file with **full OpenMM support**.

---

## 📋 Prerequisites

### 1. Install Conda (Required for OpenMM)

OpenMM **requires** Conda - it cannot be installed via pip.

**Download Miniconda:**
- Windows: https://docs.conda.io/en/latest/miniconda.html
- Install for "Just Me" or "All Users"
- ✅ Check "Add to PATH" during installation

### 2. Create Conda Environment

```bash
# Create new environment with Python 3.10
conda create -n mutalyze_build python=3.10 -y

# Activate environment
conda activate mutalyze_build
```

### 3. Install OpenMM (via Conda)

```bash
# Install OpenMM from conda-forge
conda install -c conda-forge openmm -y

# Verify installation
python -c "import openmm; print(f'OpenMM {openmm.version.version} installed ✅')"
```

### 4. Install Python Dependencies

```bash
# Install PyInstaller and other requirements
pip install -r requirements_exe.txt

# Verify PyInstaller
pyinstaller --version
```

---

## 🔨 Building the Executable

### Method 1: Quick Build (Recommended)

```bash
# Navigate to Mutalyze directory
cd /path/to/Mutalyze

# Run build script
python build_exe.py
```

**Output:** `dist/Mutalyze.exe` (~500-800 MB)

---

### Method 2: Advanced Build (Custom Spec)

```bash
# Use the advanced spec file for more control
pyinstaller build_exe_advanced.spec
```

**Output:** `dist/Mutalyze.exe`

---

### Method 3: Manual Build

```bash
pyinstaller \
  --name=Mutalyze \
  --onefile \
  --windowed \
  --add-data="models;models" \
  --add-data="data;data" \
  --add-data="forcefields;forcefields" \
  --add-data="assets;assets" \
  --add-data="core;core" \
  --add-data="app.py;." \
  --hidden-import=openmm \
  --hidden-import=openmm.app \
  --hidden-import=openmm.unit \
  --hidden-import=streamlit \
  --hidden-import=streamlit.web.cli \
  --collect-all=streamlit \
  --collect-all=openmm \
  --clean \
  --noconfirm \
  mutalyze_launcher.py
```

---

## 📦 What Gets Included

The executable includes:

- ✅ **Full OpenMM** (AMBER ff19SB force field + GBSA solvation)
- ✅ **Streamlit** web interface
- ✅ **BioPython** for PDB handling
- ✅ **scikit-learn** with Random Forest v5 calibration model
- ✅ **py3Dmol** for 3D structure visualization
- ✅ **All core modules** (mutation builder, energy calc, analysis, etc.)
- ✅ **Force fields** (AMBER parameters)
- ✅ **Calibration model** (models/calibration_v5.pkl)
- ✅ **Sample data** and templates

**Total Size:** 500-800 MB (due to OpenMM and dependencies)

---

## 🚀 Running the Executable

### Windows

1. **Double-click** `Mutalyze.exe`
2. Browser opens automatically at `http://localhost:<port>`
3. Full functionality with OpenMM-powered calculations

### Command Line (Optional)

```bash
# Run executable
./Mutalyze.exe

# Console will show:
# ✅ Starting server on port 8501...
# 🔬 OpenMM support: Enabled (full accuracy)
# 🌐 Browser will open at: http://localhost:8501
```

---

## ⚙️ Build Configuration

### Optimization Options

Edit `build_exe_advanced.spec`:

```python
# For smaller file (slower startup)
exe = EXE(
    ...
    upx=True,  # Enable UPX compression
    ...
)

# For debugging
exe = EXE(
    ...
    console=True,   # Show console window
    debug=True,     # Enable debug mode
    ...
)

# For production release
exe = EXE(
    ...
    console=False,  # Hide console
    debug=False,    # Disable debug
    upx=False,      # Better compatibility
    ...
)
```

---

## 🧪 Testing the Executable

### 1. Test OpenMM Integration

Run the exe and in the app:
1. Upload a PDB file (e.g., `1CRN.pdb`)
2. Go to **Single Point Mutation** tab
3. Select a mutation (e.g., A5G)
4. Click **Calculate**
5. Verify ΔΔG shows a value (should NOT be "disabled")

### 2. Check Console Output

The executable should show:
```
🔬 OpenMM support: Enabled (full accuracy)
✅ Using AMBER ff19SB force field
✅ Random Forest v5 calibration loaded
```

### 3. Verify Accuracy

Compare results with web deployment:
- **EXE (OpenMM):** High accuracy (r=0.837)
- **Web (Simplified):** Approximate predictions

---

## 🐛 Troubleshooting

### Issue: "OpenMM not found"

**Solution:**
```bash
# Verify OpenMM in build environment
conda activate mutalyze_build
python -c "import openmm; print('OK')"

# If fails, reinstall
conda install -c conda-forge openmm --force-reinstall
```

---

### Issue: "Streamlit not starting"

**Solution:**
```bash
# Add --console to see errors
pyinstaller build_exe_advanced.spec

# Check dist/Mutalyze.exe output
```

---

### Issue: "Executable too large (>1 GB)"

**Solutions:**
1. Enable UPX compression:
   ```python
   upx=True  # in .spec file
   ```

2. Exclude unnecessary packages:
   ```python
   excludes=['pytest', 'unittest', 'tkinter']
   ```

3. Use `--onedir` instead of `--onefile`:
   ```bash
   # Creates folder with exe + dependencies
   pyinstaller --onedir mutalyze_launcher.py
   ```

---

### Issue: "Missing DLL errors on other computers"

**Solution:**
Include Visual C++ Redistributable:
1. Download: https://aka.ms/vs/17/release/vc_redist.x64.exe
2. Install on target machine
3. Or bundle with installer

---

## 📊 Performance Comparison

| Metric | Web (Streamlit Cloud) | Desktop (.exe) |
|--------|----------------------|----------------|
| **Energy Calculator** | Simplified (statistical) | Full OpenMM (MD) |
| **Accuracy** | Approximate | r=0.837, RMSE=0.54 |
| **Speed** | ~Instant | 2-5 sec/mutation |
| **Dependencies** | pip only | Conda (bundled) |
| **Use Case** | Screening, demos | Publication-quality |

---

## 📦 Distribution

### Option 1: Direct Distribution

1. Upload `Mutalyze.exe` to cloud storage (Google Drive, Dropbox)
2. Share download link
3. Users double-click to run (no installation needed)

### Option 2: Installer (Advanced)

Create installer with [Inno Setup](https://jrsoftware.org/isinfo.php):

```iss
[Setup]
AppName=Mutalyze
AppVersion=1.0
DefaultDirName={pf}\Mutalyze
OutputDir=.
OutputBaseFilename=MutalyzeInstaller

[Files]
Source: "dist\Mutalyze.exe"; DestDir: "{app}"

[Icons]
Name: "{commondesktop}\Mutalyze"; Filename: "{app}\Mutalyze.exe"
```

---

## 🎯 Deployment Checklist

Before distributing:

- [ ] Test on clean Windows VM (no Python installed)
- [ ] Verify OpenMM calculations work
- [ ] Check all features (visualization, batch, export)
- [ ] Test with sample PDB files
- [ ] Measure startup time (<30 seconds)
- [ ] Verify file size (<800 MB)
- [ ] Test on Windows 10 and 11
- [ ] Include README with usage instructions

---

## 📝 Notes

1. **First Run:** May take 10-30 seconds to extract and start
2. **Antivirus:** Some may flag PyInstaller exes - add exclusion
3. **Updates:** Rebuild exe for new features/fixes
4. **Size:** Large due to OpenMM + ML models (normal)

---

## 🔗 Resources

- PyInstaller Docs: https://pyinstaller.org/
- OpenMM Installation: https://openmm.org/
- Streamlit Deployment: https://docs.streamlit.io/
- Conda Cheat Sheet: https://docs.conda.io/projects/conda/en/latest/user-guide/cheatsheet.html

---

## ✅ Success Criteria

Your build is successful when:

✅ Executable runs without errors  
✅ OpenMM loads and calculates ΔΔG  
✅ Browser opens automatically  
✅ All tabs functional (mutations, visualization, analysis)  
✅ Can process batch mutations  
✅ Results export to CSV  
✅ 3D visualization works  
✅ No console errors  

---

**🎉 You now have a fully standalone Mutalyze application with publication-quality OpenMM calculations!**
