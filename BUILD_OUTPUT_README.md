# 📁 Build and Distribution Files

This directory contains files generated during the build process. These files are created by PyInstaller when you run the build scripts.

## 🗂️ Directory Structure After Build

```
build/                  # Temporary build files (can be deleted)
  └── Mutalyze/        # PyInstaller working directory

dist/                   # Final output directory
  └── Mutalyze.exe     # Your standalone executable! (Windows)
  └── Mutalyze         # Your standalone executable! (Linux/Mac)

*.spec                  # PyInstaller specification files (keep)
```

---

## 📦 What to Distribute

### The Executable
**File:** `dist/Mutalyze.exe` (Windows) or `dist/Mutalyze` (Linux/Mac)  
**Size:** ~700 MB  
**This is the ONLY file users need!**

### Optional Documentation
Include these for users:
- `EXECUTABLE_README.md` - User manual
- `QUICK_START.md` - Quick start guide
- `sample_mutations.csv` - Example data

---

## 🧹 Cleanup

After successful build, you can safely delete:
- `build/` directory (temporary files)
- `__pycache__/` directories

**Keep:**
- `dist/` directory (contains your executable)
- `.spec` files (for future rebuilds)

---

## 🔄 Rebuilding

To rebuild after code changes:

1. **Clean previous build:**
   ```bash
   rm -rf build/ dist/
   ```

2. **Run build script:**
   ```bash
   build_windows.bat    # Windows
   ./build_linux.sh     # Linux/Mac
   ```

3. **New executable will be in:** `dist/`

---

## 📤 Distribution Checklist

Before distributing to users:

- [ ] Executable built successfully
- [ ] Tested on clean machine (no Python installed)
- [ ] All features work (upload PDB, calculate, visualize)
- [ ] OpenMM calculations functional (not "disabled")
- [ ] No console errors
- [ ] Included user documentation
- [ ] Version number updated (if applicable)

---

## 💾 Backup

**Recommended:** Keep a copy of the executable for each version:

```
releases/
  ├── Mutalyze_v1.0.0.exe
  ├── Mutalyze_v1.0.1.exe
  └── Mutalyze_v1.1.0.exe
```

---

## 🚫 .gitignore

These directories are already in `.gitignore`:

```gitignore
build/
dist/
*.spec
```

**Note:** Don't commit executables to git (too large). Use GitHub Releases instead.

---

## 📊 Size Optimization

If executable is too large (>800 MB):

1. **Enable UPX compression** (edit `.spec` file):
   ```python
   upx=True
   ```

2. **Exclude unnecessary packages**:
   ```python
   excludes=['pytest', 'unittest', 'tkinter']
   ```

3. **Use `--onedir` instead of `--onefile`**:
   - Creates folder with exe + dependencies
   - Faster startup, but multiple files

---

**Questions?** See `BUILD_INSTRUCTIONS.md`
