# 🎯 Mutalyze Executable - Quick Reference Card

## 🚀 Build in 3 Steps

### Windows
```cmd
1. Install Miniconda
2. Double-click: build_windows.bat
3. Get: dist\Mutalyze.exe
```

### Linux/Mac
```bash
1. Install Miniconda
2. Run: ./build_linux.sh
3. Get: dist/Mutalyze
```

---

## 📦 What You Get

**Executable Size:** ~700 MB  
**OpenMM:** ✅ Full support  
**Accuracy:** r=0.837, RMSE=0.54 kcal/mol  
**Installation:** None (standalone)  

---

## ✨ Features

- ✅ AMBER ff19SB force field
- ✅ GBSA implicit solvation  
- ✅ Random Forest v5 calibration
- ✅ 3D visualization
- ✅ Batch processing
- ✅ Offline functionality

---

## 👤 User Experience

1. Download `Mutalyze.exe`
2. Double-click to run
3. Browser opens automatically
4. No Python/Conda needed!

---

## 📚 Documentation

- **BUILD_INSTRUCTIONS.md** - Complete build guide
- **QUICK_BUILD.md** - Fast start (5 min)
- **EXECUTABLE_README.md** - User manual
- **EXE_BUILD_SUMMARY.md** - Overview

---

## 🧪 Test Your Build

```bash
python test_exe_build.py
```

All tests should pass ✅

---

## 🆚 Web vs Desktop

| Feature | Web | Desktop |
|---------|-----|---------|
| **Method** | Simplified | Full OpenMM |
| **Accuracy** | Approx | High (0.837) |
| **Speed** | Instant | 2-5 sec |
| **Offline** | ❌ | ✅ |
| **Use** | Screening | Publication |

---

## 🐛 Troubleshooting

**Won't start?**  
→ Run as administrator  
→ Check antivirus  

**Slow startup?**  
→ Normal (10-30s first run)

**OpenMM disabled?**  
→ Rebuild with conda OpenMM

---

## 📊 Performance

- Single mutation: 2-5 sec
- 10 mutations: 15-30 sec
- 100 mutations: 2-5 min

---

## 🔗 Quick Links

- Miniconda: https://docs.conda.io/en/latest/miniconda.html
- PyInstaller: https://pyinstaller.org/
- OpenMM: https://openmm.org/

---

**Ready?** Run `build_windows.bat` now! 🚀
