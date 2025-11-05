# Mutalyze - HuggingFace Spaces Deployment

## 🎯 About This Folder

This folder contains all files needed to deploy Mutalyze on **HuggingFace Spaces** with **full OpenMM support** using Conda.

## ✨ Key Advantage

**HuggingFace Spaces supports Conda!** This means:
- ✅ Full OpenMM molecular dynamics
- ✅ AMBER ff19SB force field
- ✅ Publication-quality accuracy (r=0.837)
- ✅ Same results as desktop executable
- ✅ Free hosting!

Unlike Streamlit Cloud (pip-only), HuggingFace can run the full high-accuracy version.

---

## 📦 Files in This Folder

```
huggingface_deployment/
├── README.md                 # HuggingFace Space card (this shows on Space page)
├── environment.yml           # Conda environment (CRITICAL for OpenMM)
├── app.py                    # Main Streamlit application
├── core/                     # Core analysis modules
│   ├── pdb_utils.py
│   ├── mutation_builder.py
│   ├── energy_calc.py
│   ├── empirical_correction.py
│   ├── interface_analysis.py
│   ├── conservation.py
│   ├── analysis.py
│   ├── parallel.py
│   ├── visualization.py
│   ├── simple_energy.py     # Fallback (won't be used on HF)
│   └── structure_fixer.py
├── models/
│   └── calibration_v5.pkl    # Random Forest model
├── data/                     # Sample PDB files
└── sample_mutations.csv      # Example mutations
```

---

## 🚀 Deployment Steps

### 1. Create HuggingFace Account
- Go to https://huggingface.co/
- Sign up (free)

### 2. Create New Space
1. Click "New" → "Space"
2. Name: `mutalyze` (or your choice)
3. **SDK:** Select "Streamlit"
4. **Hardware:** CPU (free) or GPU (faster, paid)
5. Click "Create Space"

### 3. Upload Files

**Option A: Via Web Interface**
1. Click "Files" tab
2. Upload all files from this folder
3. Maintain directory structure

**Option B: Via Git (Recommended)**
```bash
# Clone your Space repository
git clone https://huggingface.co/spaces/YOUR_USERNAME/mutalyze
cd mutalyze

# Copy files
cp -r /path/to/Mutalyze/huggingface_deployment/* .

# Commit and push
git add .
git commit -m "Initial deployment with OpenMM support"
git push
```

### 4. Wait for Build
- HuggingFace will detect `environment.yml`
- Conda environment will be created
- OpenMM will be installed
- App will start automatically
- First build: 10-15 minutes

### 5. Test Your Space
- Open your Space URL
- Upload a PDB (try "1CRN")
- Calculate a mutation (e.g., A5G)
- Verify ΔΔG shows a value
- Check: Should NOT say "OpenMM disabled"

---

## ⚙️ Configuration

### environment.yml
This file is **critical** for OpenMM support:

```yaml
name: mutalyze
channels:
  - conda-forge
dependencies:
  - python=3.10
  - openmm>=8.0.0          # Requires Conda!
  - biopython>=1.81
  - scikit-learn>=1.3.0
  # ... other packages
```

**Do NOT use requirements.txt** - it won't install OpenMM correctly.

### README.md (Space Card)
The README.md is displayed on your Space page. It includes:
- Title and description
- Features overview
- Quick start guide
- Accuracy metrics
- Citation information

Edit it to add:
- Your name/institution
- Your Space URL
- Contact information

---

## 🎨 Customization

### Update Space Metadata
Edit the YAML frontmatter in `README.md`:

```yaml
---
title: Your Custom Title
emoji: 🧬  # Change emoji
colorFrom: blue  # Change colors
colorTo: green
---
```

### Add Custom Domain (Pro feature)
HuggingFace Pro allows custom domains like `mutalyze.yourlab.edu`

---

## 📊 Comparison: Deployment Options

| Feature | HuggingFace Space | Streamlit Cloud | Desktop .exe |
|---------|-------------------|-----------------|--------------|
| **OpenMM** | ✅ Yes (Conda) | ❌ No (pip only) | ✅ Yes (bundled) |
| **Accuracy** | High (r=0.837) | Approximate | High (r=0.837) |
| **Speed** | 2-5 sec/mutation | ~Instant | 2-5 sec/mutation |
| **Cost** | Free | Free | One-time build |
| **Hosting** | HuggingFace | Streamlit | User's computer |
| **Best For** | Publication research | Quick screening | Offline work |

---

## 🔧 Troubleshooting

### Space Won't Start
**Check build logs:**
1. Go to your Space
2. Click "Logs" tab
3. Look for errors

**Common issues:**
- Missing files → Re-upload
- Wrong Python version → Check environment.yml
- Out of memory → Upgrade hardware tier

### OpenMM Not Working
**Symptoms:** Message says "OpenMM not available"

**Solutions:**
1. Verify `environment.yml` includes `openmm>=8.0.0`
2. Check it's under `dependencies:` (not `pip:`)
3. Rebuild: Settings → Factory Reboot

### Slow Performance
**Options:**
1. Upgrade to GPU hardware (paid)
2. Disable energy minimization (less accurate but faster)
3. Use persistent storage for model caching

---

## 💾 Space Settings

### Hardware Options
- **CPU Basic:** Free, slower (~5 sec/mutation)
- **CPU Upgrade:** Faster (~2 sec/mutation)
- **GPU:** Much faster (paid)

### Persistence
Enable persistent storage to cache:
- Downloaded PDB files
- Preprocessed structures
- Computation results

---

## 🌐 Sharing Your Space

Your Space will be at:
```
https://huggingface.co/spaces/YOUR_USERNAME/mutalyze
```

**Embed in website:**
```html
<iframe
  src="https://YOUR_USERNAME-mutalyze.hf.space"
  frameborder="0"
  width="850"
  height="450"
></iframe>
```

**Share on social media:**
- Twitter/X card preview
- LinkedIn preview
- Automatic thumbnail generation

---

## 📈 Analytics

HuggingFace provides:
- View count
- User analytics
- Error tracking
- Performance metrics

Access via Space Settings → Analytics

---

## 🔄 Updates

To update your Space:

```bash
# Make changes locally
cd /path/to/huggingface_deployment

# Edit files as needed
vim app.py

# Push to HuggingFace
git add .
git commit -m "Update: description of changes"
git push

# Space rebuilds automatically
```

---

## 🎯 Best Practices

### 1. Version Control
Keep this folder in sync with main Mutalyze repo:
```bash
# After updating main app
cp ../app.py .
cp -r ../core .
git add . && git commit -m "Sync with main repo"
```

### 2. Testing
Test locally before pushing:
```bash
# Create conda env from environment.yml
conda env create -f environment.yml
conda activate mutalyze

# Run locally
streamlit run app.py

# Verify OpenMM works
python -c "import openmm; print('OK')"
```

### 3. Documentation
Keep README.md updated with:
- Recent changes
- Known limitations
- Contact information
- Citation details

---

## 🆘 Support

**HuggingFace Docs:**
- Spaces: https://huggingface.co/docs/hub/spaces
- Conda: https://huggingface.co/docs/hub/spaces-dependencies

**Mutalyze Issues:**
- GitHub: [Your repo]/issues
- Email: [Your email]

---

## 🎉 Success Checklist

Before going public:

- [ ] Space builds successfully
- [ ] OpenMM loads (check logs)
- [ ] Can upload PDB files
- [ ] Mutations calculate correctly
- [ ] ΔΔG values are reasonable
- [ ] 3D visualization works
- [ ] Export functions work
- [ ] README.md has your info
- [ ] Space card looks professional
- [ ] Tested on mobile/tablet

---

## 📊 Expected Performance

With CPU (free tier):
- Structure loading: <1 second
- Single mutation: 3-5 seconds
- 10 mutations: 30-50 seconds
- 100 mutations: 5-8 minutes

With GPU (paid):
- Single mutation: 1-2 seconds
- 10 mutations: 10-20 seconds
- 100 mutations: 2-3 minutes

---

## 🌟 Pro Tips

### 1. Add Badge to Main Repo
```markdown
[![Open in Spaces](https://huggingface.co/datasets/huggingface/badges/raw/main/open-in-hf-spaces-sm.svg)](https://huggingface.co/spaces/YOUR_USERNAME/mutalyze)
```

### 2. Enable Community Features
- Discussions tab for user questions
- Comments on your Space
- Duplicate feature for users to clone

### 3. Showcase Your Work
- Add to HuggingFace papers
- Share on community Discord
- Feature in blog posts

---

**🚀 Ready to deploy? Upload these files to your HuggingFace Space!**

**Questions?** See HuggingFace Spaces documentation or open an issue on GitHub.

---

## 📝 Notes

- **Build time:** First deployment takes 10-15 minutes (conda environment)
- **Subsequent builds:** 2-5 minutes (cache used)
- **Storage limit:** 50 GB (free tier)
- **Compute limit:** Fair use policy (free tier)

For production with heavy traffic, consider upgrading to Pro tier.
