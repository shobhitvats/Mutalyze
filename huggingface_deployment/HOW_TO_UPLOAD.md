# 🤗 How to Upload to HuggingFace Spaces

Complete step-by-step guide for deploying Mutalyze with full OpenMM support.

---

## 📋 Prerequisites

1. **HuggingFace Account**
   - Go to: https://huggingface.co/join
   - Sign up (free)
   - Verify your email

2. **Git Installed** (for Option B)
   - Check: `git --version`
   - Install if needed: https://git-scm.com/

3. **Git LFS** (Large File Storage)
   - Required for ML model files
   - Install: https://git-lfs.com/
   - Or: `git lfs install`

---

## 🎯 Method 1: Via Git (Recommended)

### Step 1: Create Your Space on HuggingFace

1. **Go to HuggingFace:**
   - Visit: https://huggingface.co/new-space

2. **Configure your Space:**
   ```
   Owner: YOUR_USERNAME
   Space name: mutalyze
   License: mit (or your choice)
   Select SDK: Streamlit
   Space hardware: CPU basic (free)
   ```

3. **Click "Create Space"**
   - Your Space is created but empty
   - You'll see instructions to clone it

### Step 2: Clone Your HuggingFace Space

Open terminal and run:

```bash
# Make sure git-lfs is installed
git lfs install

# Clone your Space (replace YOUR_USERNAME with your HuggingFace username)
git clone https://huggingface.co/spaces/YOUR_USERNAME/mutalyze
cd mutalyze
```

**Example:**
If your username is `johndoe`, the command would be:
```bash
git clone https://huggingface.co/spaces/johndoe/mutalyze
cd mutalyze
```

### Step 3: Copy Files from huggingface_deployment/

```bash
# Navigate to your Mutalyze project
cd /home/vats/Desktop/Mutalyze

# Copy all files from huggingface_deployment to your Space
cp -r huggingface_deployment/* /path/to/mutalyze/

# Example if Space is in home directory:
cp -r huggingface_deployment/* ~/mutalyze/
```

**Or use full path:**
```bash
# From your Space directory
cd ~/mutalyze  # Or wherever you cloned it

# Copy files
cp -r /home/vats/Desktop/Mutalyze/huggingface_deployment/* .
```

### Step 4: Track Large Files with Git LFS

The ML model file is large, so we need Git LFS:

```bash
cd ~/mutalyze  # Your Space directory

# Track the model file with Git LFS
git lfs track "models/calibration_v5.pkl"
git lfs track "*.pkl"

# Add the LFS tracking file
git add .gitattributes
```

### Step 5: Commit and Push to HuggingFace

```bash
# Check what files will be added
git status

# Add all files
git add .

# Commit with a message
git commit -m "Initial deployment with OpenMM support"

# Push to HuggingFace
git push
```

**Authentication:**
- You'll be prompted for credentials
- **Username:** Your HuggingFace username
- **Password:** Use an **Access Token** (NOT your password)

**Create Access Token:**
1. Go to: https://huggingface.co/settings/tokens
2. Click "New token"
3. Name: `mutalyze-deployment`
4. Role: `write`
5. Copy the token
6. Use it as password when pushing

### Step 6: Wait for Build

1. **Go to your Space:**
   ```
   https://huggingface.co/spaces/YOUR_USERNAME/mutalyze
   ```

2. **Check "Logs" tab:**
   - You'll see the build process
   - Conda environment being created
   - OpenMM being installed
   - Time: 10-15 minutes (first build)

3. **Build output should show:**
   ```
   Collecting package metadata...
   Solving environment...
   Installing openmm...
   Installing streamlit...
   Starting app...
   ```

4. **When ready:**
   - Status changes to "Running"
   - App tab becomes available
   - Your Space is live! 🎉

---

## 🎯 Method 2: Via Web Interface (Simpler but Manual)

### Step 1: Create Space (Same as Method 1)
- https://huggingface.co/new-space
- Configure as described above

### Step 2: Upload Files Manually

1. **Click "Files" tab** in your Space

2. **Click "Add file" → "Upload files"**

3. **Upload these files/folders:**
   ```
   From huggingface_deployment/:
   ✅ environment.yml (CRITICAL!)
   ✅ README.md
   ✅ app.py
   ✅ .gitignore
   ✅ DEPLOYMENT_GUIDE.md
   ✅ QUICKSTART.md
   ✅ test_deployment.sh
   ✅ sample_mutations.csv
   ```

4. **Create folders and upload:**
   
   **Create "core" folder:**
   - Click "Add file" → "Create a new file"
   - Filename: `core/.gitkeep`
   - Commit
   - Now upload all files from `huggingface_deployment/core/`
   
   **Create "models" folder:**
   - Same process
   - Upload `calibration_v5.pkl`
   
   **Create "data" folder:**
   - Upload files from `huggingface_deployment/data/`
   
   **Create "forcefields" folder:**
   - Upload files from `huggingface_deployment/forcefields/`

5. **Commit changes after each upload**

### Step 3: Wait for Build (Same as Method 1)

---

## 🧪 Testing Your Space

### After Build Completes:

1. **Click "App" tab** in your Space

2. **Test basic functionality:**
   - Interface loads ✅
   - No error messages ✅

3. **Test OpenMM:**
   - Click "Single Point Mutation" tab
   - Enter PDB ID: `1CRN`
   - Click "Fetch from RCSB"
   - Select position: `5` (Alanine)
   - Select mutation: `GLY` (A5G)
   - Click "Calculate"
   - **Expected:** ΔΔG value appears (e.g., ~+0.8 kcal/mol)
   - **Check:** Should NOT say "OpenMM disabled"

4. **Check console/logs:**
   - Go to "Logs" tab
   - Look for: `✅ OPENMM_AVAILABLE = True`
   - Should see: `Using AMBER ff19SB force field`

---

## 🔧 Updating Your Space

### Via Git:

```bash
# Navigate to your Space
cd ~/mutalyze

# Pull latest from HuggingFace (in case of any changes)
git pull

# Make changes to files
# ... edit app.py or other files ...

# Commit and push
git add .
git commit -m "Update: description of changes"
git push
```

### Via Web Interface:

1. Go to your Space
2. Click "Files" tab
3. Click on file to edit
4. Click "Edit" button
5. Make changes
6. Commit with description

---

## 🐛 Troubleshooting

### Problem: "OpenMM not available" in app

**Solution:**
1. Check `environment.yml` is in root directory
2. Verify it contains:
   ```yaml
   dependencies:
     - openmm>=8.0.0
   ```
3. Factory reboot: Settings → Factory Reboot

### Problem: Build fails with "environment.yml not found"

**Solution:**
- Make sure `environment.yml` is in the **root** directory (not in a subfolder)
- Check filename is exactly `environment.yml` (not `environment.yaml`)

### Problem: Large file upload fails

**Solution:**
- Use Git LFS for files >10MB
- For web upload: Use Git method instead

### Problem: "Authentication failed"

**Solution:**
- Use Access Token, not password
- Create token: https://huggingface.co/settings/tokens
- Token needs `write` permission

### Problem: Build is slow/hanging

**Solution:**
- First build takes 10-15 minutes (normal)
- Check "Logs" tab for progress
- If stuck >20 minutes: Factory Reboot

---

## 📊 Expected Build Process

```
1. Detecting dependencies...
   ✅ Found environment.yml
   ✅ Using Conda

2. Creating Conda environment...
   ⏳ This takes 5-10 minutes
   ✅ Environment created

3. Installing dependencies...
   ⏳ Installing OpenMM (largest package)
   ⏳ Installing BioPython, scikit-learn, etc.
   ✅ All packages installed

4. Starting Streamlit app...
   ✅ App running on port 7860

5. Space ready!
   🎉 Your Space is live
```

---

## ✅ Verification Checklist

After deployment:

- [ ] Space shows "Running" status
- [ ] App tab loads without errors
- [ ] Can upload PDB file
- [ ] Can fetch PDB from RCSB (try "1CRN")
- [ ] Single Point Mutation tab works
- [ ] Can calculate a mutation
- [ ] ΔΔG value appears (not "disabled")
- [ ] Logs show "OPENMM_AVAILABLE = True"
- [ ] 3D visualization works
- [ ] Can export results (CSV)
- [ ] README displays correctly on Space page

---

## 🎨 Customizing Your Space

### Update Space Card (README.md)

Edit the YAML frontmatter:

```yaml
---
title: Mutalyze - Protein Mutation Analysis
emoji: 🧬
colorFrom: blue
colorTo: green
sdk: streamlit
sdk_version: "1.28.0"
app_file: app.py
pinned: false
---
```

**Change:**
- `title` - Your Space title
- `emoji` - Icon (any emoji)
- `colorFrom/colorTo` - Gradient colors
- `pinned` - Pin to your profile

### Add Space Description

Edit `README.md` content below the frontmatter to describe your tool.

---

## 💡 Pro Tips

### 1. Use Git for Updates
Git method is faster for updates than web interface.

### 2. Test Locally First
Before pushing to HuggingFace:
```bash
cd huggingface_deployment
./test_deployment.sh
```

### 3. Monitor Build Logs
Always check "Logs" tab during first build to catch issues early.

### 4. Use Private Spaces
For development, create private Space (Settings → Visibility)

### 5. Enable Discussions
Settings → Enable discussions for user feedback

---

## 📖 Quick Reference Commands

### Initial Setup:
```bash
# Clone your Space
git clone https://huggingface.co/spaces/YOUR_USERNAME/mutalyze
cd mutalyze

# Copy files
cp -r /home/vats/Desktop/Mutalyze/huggingface_deployment/* .

# Setup Git LFS
git lfs install
git lfs track "*.pkl"

# Push to HuggingFace
git add .
git commit -m "Initial deployment"
git push
```

### Updates:
```bash
cd ~/mutalyze
# ... make changes ...
git add .
git commit -m "Update: description"
git push
```

---

## 🔗 Useful Links

- **Your Space:** https://huggingface.co/spaces/YOUR_USERNAME/mutalyze
- **HuggingFace Docs:** https://huggingface.co/docs/hub/spaces
- **Streamlit on HF:** https://huggingface.co/docs/hub/spaces-sdks-streamlit
- **Git LFS:** https://git-lfs.com/
- **Access Tokens:** https://huggingface.co/settings/tokens

---

## 🎉 Success!

Once deployed, your Space URL will be:
```
https://huggingface.co/spaces/YOUR_USERNAME/mutalyze
```

**Share it in:**
- Research papers
- GitHub README
- Conference presentations
- Social media

**Example badge for GitHub README:**
```markdown
[![Open in Spaces](https://huggingface.co/datasets/huggingface/badges/raw/main/open-in-hf-spaces-sm.svg)](https://huggingface.co/spaces/YOUR_USERNAME/mutalyze)
```

---

**Need help?** Check the troubleshooting section or the detailed guide in `DEPLOYMENT_GUIDE.md`

**Ready to deploy?** Follow Method 1 (Git) for best results! 🚀
