# 🤗 HuggingFace Spaces - Quick Start

## ⚡ Deploy in 5 Steps

### 1️⃣ Create HuggingFace Account
Go to: https://huggingface.co/join

### 2️⃣ Create New Space
- Click "New" → "Space"
- **Name:** `mutalyze`
- **SDK:** Streamlit
- **Hardware:** CPU (free)

### 3️⃣ Upload Files
Upload ALL files from `huggingface_deployment/` folder to your Space.

**Critical files:**
- ✅ `environment.yml` (enables OpenMM via Conda)
- ✅ `app.py` (main application)
- ✅ `README.md` (Space description)
- ✅ `core/` folder (analysis modules)
- ✅ `models/` folder (ML model)
- ✅ `data/` folder (sample files)

### 4️⃣ Wait for Build
- HuggingFace auto-detects `environment.yml`
- Creates Conda environment
- Installs OpenMM
- Builds app
- **Time:** 10-15 minutes (first time)

### 5️⃣ Test Your Space
- Open your Space URL
- Upload PDB: `1CRN`
- Calculate: `A5G`
- Verify: ΔΔG value appears ✅

---

## 🎯 What You Get

**Your Space URL:**
```
https://huggingface.co/spaces/YOUR_USERNAME/mutalyze
```

**Features:**
- ✅ Full OpenMM (AMBER ff19SB + GBSA)
- ✅ Publication-quality accuracy (r=0.837)
- ✅ Free hosting
- ✅ Automatic SSL/HTTPS
- ✅ Share with anyone

---

## 🆚 Why HuggingFace?

| Feature | HuggingFace | Streamlit Cloud |
|---------|-------------|-----------------|
| **Conda Support** | ✅ Yes | ❌ No |
| **OpenMM** | ✅ Full | ❌ Simplified |
| **Accuracy** | High (0.837) | Approximate |
| **Free Tier** | ✅ Yes | ✅ Yes |

**Bottom Line:** HuggingFace = Best of both worlds (cloud + accuracy)

---

## 📋 Checklist

Before deploying:
- [ ] All files copied to Space
- [ ] `environment.yml` included
- [ ] Updated README.md with your info
- [ ] Tested locally (run `test_deployment.sh`)

After deploying:
- [ ] Space builds successfully
- [ ] OpenMM loads (check logs)
- [ ] Mutations calculate correctly
- [ ] Share your Space!

---

## 🐛 Troubleshooting

**Space won't build?**
→ Check Logs tab for errors

**OpenMM not working?**
→ Verify `environment.yml` has `openmm>=8.0.0`

**Slow performance?**
→ Upgrade to GPU hardware (paid tier)

---

## 📚 Full Guide

See `DEPLOYMENT_GUIDE.md` for complete instructions.

---

**🚀 Ready? Create your Space now!**

https://huggingface.co/new-space
