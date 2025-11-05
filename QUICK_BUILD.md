# 🚀 Quick Start: Build Mutalyze Executable

This guide gets you from zero to a working `.exe` file in 10 minutes.

---

## ⚡ Quick Build (Windows)

### Step 1: Install Conda

Download and install Miniconda:
```
https://docs.conda.io/en/latest/miniconda.html
```

✅ Check "Add Conda to PATH" during installation

---

### Step 2: Run Build Script

Open Command Prompt in the Mutalyze folder and run:

```cmd
build_windows.bat
```

**That's it!** The script will:
1. Create conda environment
2. Install OpenMM + dependencies
3. Build the executable

---

### Step 3: Run Your Executable

```cmd
dist\Mutalyze.exe
```

Browser opens automatically with full OpenMM support! 🎉

---

## 🐧 Quick Build (Linux/Mac)

```bash
# Make script executable
chmod +x build_linux.sh

# Run build
./build_linux.sh

# Run executable
./dist/Mutalyze
```

---

## 🔧 Manual Build (If Scripts Fail)

### 1. Setup Environment

```bash
# Create environment
conda create -n mutalyze_build python=3.10 -y

# Activate
conda activate mutalyze_build

# Install OpenMM
conda install -c conda-forge openmm -y

# Install dependencies
pip install -r requirements_exe.txt
```

### 2. Build Executable

```bash
# Quick build
pyinstaller build_exe_advanced.spec
```

### 3. Find Your Executable

**Windows:** `dist\Mutalyze.exe`  
**Linux/Mac:** `dist/Mutalyze`

---

## ✅ Verify Build

Run the executable and check:

1. ✅ Browser opens automatically
2. ✅ Upload a PDB file (try `data/1CRN.pdb`)
3. ✅ Go to "Single Point Mutation" tab
4. ✅ Select mutation (e.g., A5G)
5. ✅ Click "Calculate"
6. ✅ See ΔΔG value (should NOT say "disabled")

**Success!** You have a working executable with full OpenMM support.

---

## 📊 What You Get

| Feature | Included |
|---------|----------|
| OpenMM (AMBER ff19SB) | ✅ Yes |
| Random Forest v5 | ✅ Yes |
| 3D Visualization | ✅ Yes |
| Batch Processing | ✅ Yes |
| Conservation Analysis | ✅ Yes |
| Interface Analysis | ✅ Yes |
| No Python Required | ✅ Yes |

**File Size:** ~500-800 MB (includes all dependencies)

---

## 🐛 Common Issues

### "Conda not found"
**Fix:** Restart terminal after installing Miniconda

### "OpenMM import failed"
**Fix:** 
```bash
conda activate mutalyze_build
conda install -c conda-forge openmm --force-reinstall
```

### "PyInstaller not found"
**Fix:**
```bash
pip install pyinstaller
```

### Executable too slow to start
**Normal:** First run extracts files (10-30 seconds)

---

## 📦 Distribution

**Share your executable:**

1. **Upload to Google Drive/Dropbox**
   - Just the `.exe` file (no other files needed)
   - Users download and double-click to run

2. **GitHub Release**
   - Create release with `.exe` as asset
   - Users download from releases page

3. **Direct Transfer**
   - Copy `.exe` to USB drive
   - Share via email (may need to zip first)

**Users need:** Just Windows (no Python, no conda, nothing else!)

---

## 🎯 Next Steps

1. Test on clean Windows machine (VM recommended)
2. Create icon.ico for professional look
3. Build installer with Inno Setup (optional)
4. Document for end users

---

**Need help?** See `BUILD_INSTRUCTIONS.md` for detailed guide.

**Ready to build?** Run `build_windows.bat` now! 🚀
