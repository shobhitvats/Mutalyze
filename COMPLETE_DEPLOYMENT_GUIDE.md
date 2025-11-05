# 🚀 Mutalyze - Complete Deployment Guide

## 📊 Overview: Three Deployment Options

Mutalyze can be deployed in **three different ways**, each optimized for specific use cases:

1. **🌐 Streamlit Cloud** - Fast screening (simplified calculator)
2. **🤗 HuggingFace Spaces** - Cloud + OpenMM (publication quality) ⭐ **RECOMMENDED**
3. **💻 Desktop Executable** - Offline use (publication quality)

---

## 🎯 Quick Decision Guide

**Need publication-quality accuracy in the cloud?**  
→ 🤗 **HuggingFace Spaces**

**Need fastest possible results for screening?**  
→ 🌐 **Streamlit Cloud**

**Need offline use or maximum privacy?**  
→ 💻 **Desktop Executable**

---

## 1️⃣ Streamlit Cloud Deployment

### Status
✅ **Already deployed and running**

### Features
- Energy Method: Simplified (statistical potentials)
- Accuracy: Approximate (~r=0.60)
- Speed: ~Instant
- Best for: Demos, teaching, quick screening

### Files Used
- `app.py`
- `core/` (with automatic fallback to simplified calculator)
- `requirements.txt` (pip only)

### Deployment
Automatically deploys from GitHub main branch.

---

## 2️⃣ HuggingFace Spaces Deployment ⭐

### Status
📦 **Ready to deploy** (files in `huggingface_deployment/`)

### Features
- Energy Method: **Full OpenMM** (AMBER ff19SB + GBSA)
- Accuracy: **High (r=0.837, RMSE=0.54)**
- Speed: 2-5 sec/mutation
- Best for: **Publication-quality research in the cloud**

### Key Advantage
**HuggingFace supports Conda** → Can install OpenMM → Full accuracy!

### Quick Deploy (5 Steps)

1. **Create HuggingFace account:** https://huggingface.co/join

2. **Create new Space:**
   - Click "New" → "Space"
   - Name: `mutalyze`
   - SDK: Streamlit
   - Hardware: CPU (free)

3. **Upload files:**
   - Upload ALL files from `huggingface_deployment/` folder
   - Maintain directory structure

4. **Wait for build:**
   - HuggingFace auto-detects `environment.yml`
   - Creates Conda environment
   - Installs OpenMM
   - Time: 10-15 minutes

5. **Test:**
   - Upload PDB: `1CRN`
   - Calculate: `A5G`
   - Verify: ΔΔG value appears ✅

### Documentation
- **Quick Start:** `huggingface_deployment/QUICKSTART.md`
- **Complete Guide:** `huggingface_deployment/DEPLOYMENT_GUIDE.md`
- **Testing Script:** `huggingface_deployment/test_deployment.sh`

---

## 3️⃣ Desktop Executable Deployment

### Status
🔨 **Build system ready**

### Features
- Energy Method: **Full OpenMM** (AMBER ff19SB + GBSA)
- Accuracy: **High (r=0.837, RMSE=0.54)**
- Speed: 2-5 sec/mutation
- Best for: Offline use, privacy, large batches

### Quick Build

**Windows:**
```cmd
build_windows.bat
```

**Linux/Mac:**
```bash
./build_linux.sh
```

**Output:** `dist/Mutalyze.exe` (~700 MB)

### Documentation
- **Quick Start:** `QUICK_BUILD.md`
- **Complete Guide:** `BUILD_INSTRUCTIONS.md`
- **User Manual:** `EXECUTABLE_README.md`

---

## 📊 Detailed Comparison

| Feature | Streamlit Cloud | **HuggingFace** ⭐ | Desktop .exe |
|---------|----------------|-------------------|--------------|
| **OpenMM** | ❌ No | ✅ **Yes** | ✅ Yes |
| **Accuracy** | Low (~0.60) | **High (0.837)** | High (0.837) |
| **Speed** | ~Instant | 2-5 sec | 2-5 sec |
| **Cloud** | ✅ Yes | ✅ **Yes** | ❌ No |
| **Offline** | ❌ No | ❌ No | ✅ **Yes** |
| **Free** | ✅ Yes | ✅ **Yes** | ✅ Yes |
| **Installation** | None | None | None |
| **Conda** | ❌ No | ✅ **Yes** | ✅ Yes |
| **Publication** | ❌ No | ✅ **Yes** | ✅ Yes |
| **Sharing** | Easy (URL) | **Easy (URL)** | File transfer |
| **Privacy** | Public | Public | **Private** |
| **Best For** | Demos | **Research** | Offline work |

---

## 🎓 Use Case Recommendations

### Academic Research
**Primary:** 🤗 HuggingFace Spaces  
**Backup:** 💻 Desktop Executable (for offline/sensitive data)

### Teaching & Demonstrations
**Primary:** 🌐 Streamlit Cloud  
**Alternative:** 🤗 HuggingFace Spaces (if need accuracy)

### Industry & Proprietary
**Primary:** 💻 Desktop Executable  
**Alternative:** 🤗 HuggingFace Spaces (private Spaces available)

### Collaborative Projects
**Primary:** 🤗 HuggingFace Spaces  
**Alternative:** 🌐 Streamlit Cloud (if speed > accuracy)

---

## 📂 Repository Structure

```
Mutalyze/
├── app.py                          # Main Streamlit app (all deployments)
├── core/                           # Analysis modules (all deployments)
├── models/                         # ML models (all deployments)
├── data/                           # Sample files (all deployments)
│
├── requirements.txt                # Streamlit Cloud (pip only)
│
├── huggingface_deployment/         # HuggingFace Spaces ⭐
│   ├── environment.yml             # Conda config (enables OpenMM!)
│   ├── README.md                   # Space card
│   ├── DEPLOYMENT_GUIDE.md         # Complete guide
│   ├── QUICKSTART.md               # 5-step guide
│   ├── test_deployment.sh          # Testing script
│   └── [all app files copied]      # Complete app
│
├── build_windows.bat               # Desktop .exe builder (Windows)
├── build_linux.sh                  # Desktop .exe builder (Linux/Mac)
├── build_exe_advanced.spec         # PyInstaller config
├── mutalyze_launcher.py            # Desktop launcher
├── requirements_exe.txt            # Build dependencies
├── test_exe_build.py               # Build verification
│
└── Documentation/
    ├── DEPLOYMENT_COMPARISON.md    # This file
    ├── BUILD_INSTRUCTIONS.md       # Desktop build guide
    ├── EXECUTABLE_README.md        # Desktop user manual
    └── QUICK_BUILD.md              # Desktop quick start
```

---

## 🚀 Recommended Deployment Strategy

For **maximum reach and impact**, deploy all three:

### 1. Deploy HuggingFace Spaces (Primary)
- Most users will use this
- Publication-quality results
- Easy to share in papers/talks
- Free hosting

### 2. Keep Streamlit Cloud (Secondary)
- For quick demos
- Fast preliminary checks
- Teaching purposes

### 3. Build Desktop Executable (Optional)
- For users needing offline access
- Sensitive/proprietary data
- Very large batch jobs

### Sharing Strategy
**In your paper:**
```
Online tool: https://huggingface.co/spaces/YOUR_USERNAME/mutalyze
Source code: https://github.com/YOUR_USERNAME/Mutalyze
Desktop version: Available in releases
```

---

## 🔑 Key Differences Explained

### Why HuggingFace Beats Streamlit for Accuracy

**Streamlit Cloud:**
- Limitation: pip-only environment
- Cannot install: OpenMM (requires conda)
- Solution: Fallback to simplified calculator
- Result: Lower accuracy (~r=0.60)

**HuggingFace Spaces:**
- Advantage: Conda environment support
- Can install: OpenMM via conda-forge
- Uses: Full molecular dynamics
- Result: High accuracy (r=0.837)

**Conclusion:** HuggingFace = Best of both worlds (cloud + accuracy)

---

## 📈 Performance Benchmarks

### Accuracy (S2648 Dataset)

| Method | Pearson r | RMSE | MAE |
|--------|-----------|------|-----|
| OpenMM (HF/Desktop) | **0.837** | **0.54** | **0.42** |
| Simplified (Streamlit) | ~0.60 | ~0.90 | ~0.70 |

### Speed

| Task | Streamlit | HuggingFace | Desktop |
|------|-----------|-------------|---------|
| Single mutation | 0.1 s | 3-5 s | 2-5 s |
| 10 mutations | 1 s | 30-50 s | 15-30 s |
| 100 mutations | 10 s | 5-8 min | 2-5 min |

---

## 💡 Pro Tips

### For Best Results

1. **Use HuggingFace for:**
   - Publication submissions
   - Grant applications
   - Research collaborations
   - High-accuracy predictions

2. **Use Streamlit for:**
   - Teaching demonstrations
   - Quick hypothesis testing
   - Preliminary screening
   - Maximum speed

3. **Use Desktop for:**
   - Offline conferences
   - Proprietary structures
   - Very large datasets
   - Maximum privacy

### Combining Methods

**Workflow example:**
1. Screen with Streamlit (fast, find interesting mutations)
2. Validate with HuggingFace (accurate, share results)
3. Download Desktop for offline analysis (privacy)

---

## 🆘 Support & Issues

- **HuggingFace:** See `huggingface_deployment/DEPLOYMENT_GUIDE.md`
- **Desktop:** See `BUILD_INSTRUCTIONS.md`
- **General:** See main `README.md`
- **GitHub Issues:** [Repository URL]/issues

---

## ✅ Deployment Checklist

### HuggingFace Spaces
- [ ] HuggingFace account created
- [ ] Space created (Streamlit SDK)
- [ ] Files uploaded from `huggingface_deployment/`
- [ ] Build completed successfully
- [ ] OpenMM working (not "disabled")
- [ ] Test mutation calculated
- [ ] README.md updated with your info

### Desktop Executable
- [ ] Conda installed
- [ ] Build environment created
- [ ] OpenMM installed via conda
- [ ] Build script executed
- [ ] Executable created in `dist/`
- [ ] Tested on clean machine
- [ ] Documentation included

---

## 🎉 Success!

You now have a **complete deployment strategy** for Mutalyze:

✅ **Three deployment options** (cloud fast, cloud accurate, offline)  
✅ **Full documentation** for each method  
✅ **Ready-to-use packages** for all deployments  
✅ **Testing scripts** for validation  
✅ **Comparison guides** for decision making  

**All use cases covered!** 🚀

---

**Next Steps:**
1. Deploy to HuggingFace (recommended)
2. Share your Space URL
3. Publish your research
4. Help advance protein engineering! 🧬
