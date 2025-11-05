---
title: Mutalyze - Protein Mutation Stability Analysis
emoji: 🧬
colorFrom: blue
colorTo: green
sdk: streamlit
sdk_version: "1.28.0"
app_file: app.py
pinned: false
license: mit
python_version: "3.10"
---

# Mutalyze - Protein Mutation Stability Analysis Platform

**Live Demo:** [Your HuggingFace Space URL]

## 🎯 Overview

Mutalyze is a comprehensive platform for predicting the stability effects of protein mutations using:

- **Full OpenMM** molecular dynamics (AMBER ff19SB + GBSA)
- **Random Forest v5** ML calibration (r=0.837, RMSE=0.54 kcal/mol)
- **Publication-quality** accuracy

## ✨ Features

### 🔬 Mutation Analysis
- Single point mutations
- Systematic scanning (enumerate all mutations)
- Batch processing from CSV
- Real-time ΔΔG predictions

### 📊 Advanced Analysis
- Interface residue detection
- Conservation scoring
- Structure validation
- Interactive 3D visualization

### 💾 Export Options
- CSV results
- JSON data
- Publication-ready plots

## 🚀 Quick Start

1. **Upload PDB file** or enter PDB ID (e.g., "1CRN")
2. **Select mutation** (or enumerate all)
3. **Calculate ΔΔG** (2-5 seconds per mutation)
4. **Analyze results** and visualize

## 📊 Accuracy

Benchmarked on S2648 dataset:
- **Pearson r:** 0.837
- **RMSE:** 0.54 kcal/mol
- **MAE:** 0.42 kcal/mol

## 🔬 Methodology

### Energy Calculation
1. **Structure Preparation:** PDB cleaning, missing atoms handled
2. **Force Field:** AMBER ff19SB (state-of-the-art)
3. **Solvation:** GBSA implicit solvent
4. **Minimization:** Energy minimization (optional but recommended)
5. **Calibration:** Random Forest v5 model

### Interpretation
- **ΔΔG < 0:** Stabilizing (protein more stable)
- **ΔΔG = 0:** Neutral (no effect)
- **ΔΔG > 0:** Destabilizing (protein less stable)

## 💡 Use Cases

- **Research:** Rational protein design, stability engineering
- **Industry:** Therapeutic antibodies, enzyme optimization
- **Education:** Structure-function relationships
- **Validation:** Compare predictions with experiments

## ⚙️ Technology Stack

- **Backend:** OpenMM (molecular dynamics)
- **ML:** scikit-learn (Random Forest)
- **Bioinformatics:** BioPython
- **Frontend:** Streamlit
- **Visualization:** py3Dmol, Matplotlib, Seaborn

## 🌐 Other Deployment Options

- **Web (Streamlit Cloud):** Fast screening with simplified calculator
- **Desktop (.exe):** Offline use, same accuracy as this Space
- **GitHub:** Full source code and documentation

## 📚 Citation

If you use Mutalyze in your research, please cite:

```
Mutalyze: A comprehensive protein mutation stability analysis platform
[Your Name/Institution]
2025
```

**Methods:**
- OpenMM: DOI:10.1371/journal.pcbi.1005659
- AMBER ff19SB: DOI:10.1021/acs.jctc.9b00591

## 📖 Documentation

Full documentation available in the app's "About" tab.

## 🐛 Issues & Support

- **GitHub:** [Repository URL]
- **Email:** [Your email]

## ⚖️ License

MIT License - See LICENSE file for details

---

**Note:** This HuggingFace Space uses full OpenMM with Conda environment for publication-quality results. For fastest screening, try our Streamlit Cloud deployment with simplified calculator.
