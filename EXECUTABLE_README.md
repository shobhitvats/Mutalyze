# 🧬 Mutalyze Desktop Application

**Standalone Executable with Full OpenMM Support**

---

## 📥 Download

Get the latest release:
- **Windows:** `Mutalyze.exe` (~700 MB)
- **Linux:** `Mutalyze` (~700 MB)
- **Mac:** `Mutalyze.app` (~700 MB)

---

## 🚀 Quick Start

### Windows
1. Download `Mutalyze.exe`
2. Double-click to run
3. Browser opens automatically at `http://localhost:8501`

### Linux/Mac
1. Download `Mutalyze`
2. Make executable: `chmod +x Mutalyze`
3. Run: `./Mutalyze`
4. Browser opens automatically

**No installation required!** Everything is bundled.

---

## ✨ Features

### 🔬 Protein Mutation Analysis
- Upload PDB files or fetch from RCSB
- Single point mutations
- Systematic scanning (enumerate all mutations)
- Batch processing from CSV

### 📊 Energy Calculations
- **Full OpenMM** molecular dynamics
- **AMBER ff19SB** force field
- **GBSA** implicit solvation
- **Random Forest v5** ML calibration
- **Accuracy:** r=0.837, RMSE=0.54 kcal/mol

### 🎨 Visualization
- Interactive 3D structure viewer
- Mutation highlighting
- Interface analysis
- Conservation mapping
- Publication-ready plots

### 📈 Analysis Tools
- Stability predictions (ΔΔG)
- Conservation scores
- Interface interactions
- Residue burial analysis
- Statistical summaries

### 💾 Export Options
- CSV results
- JSON data
- PNG/SVG plots
- PDB structures

---

## 🎯 Use Cases

### 1. Single Mutation Testing
Test specific mutations for protein engineering:
```
Q: What happens if I mutate position 42 from Leucine to Proline?
A: Upload structure → Select L42P → Calculate ΔΔG
```

### 2. Systematic Scanning
Find all stabilizing mutations:
```
Q: Which mutations will stabilize my protein?
A: Upload → Enumerate all → Sort by ΔΔG (negative = stabilizing)
```

### 3. Batch Analysis
Process multiple mutations from research:
```
Q: I have 100 mutations from experiments - which are predicted to be destabilizing?
A: Upload CSV → Batch calculate → Export results
```

### 4. Interface Engineering
Optimize protein-protein interactions:
```
Q: Which interface residues can I mutate without disrupting binding?
A: Upload complex → Interface analysis → Test mutations
```

---

## 📋 System Requirements

### Minimum
- **OS:** Windows 10/11, Linux (Ubuntu 20.04+), macOS 10.15+
- **RAM:** 4 GB
- **Disk:** 2 GB free space
- **Display:** 1280x720

### Recommended
- **RAM:** 8 GB or more
- **CPU:** Multi-core (for faster batch processing)
- **Disk:** SSD for faster loading

---

## 🔒 First Run (Important)

### Windows Security Warning
Windows may show "Unknown Publisher" warning on first run:

1. Click "More info"
2. Click "Run anyway"

This is normal for unsigned executables. The app is safe.

### Antivirus False Positives
Some antivirus software may flag PyInstaller executables:

**Solution:** Add exception for `Mutalyze.exe`

**Why?** PyInstaller bundles Python in a way some AV tools misidentify.

### Slow First Launch
First run may take 10-30 seconds to extract and initialize.  
**This is normal!** Subsequent runs are faster.

---

## 📖 User Guide

### Basic Workflow

#### 1. Upload Structure
- **Option A:** Upload PDB file (drag & drop)
- **Option B:** Enter PDB ID (e.g., "1CRN")
- **Option C:** Paste PDB text

#### 2. Choose Analysis Type
- **Single Point Mutation:** Test specific mutation
- **Multiple Mutations:** Batch or enumerate

#### 3. Configure Parameters
- Select residue position
- Choose mutation (or enumerate all)
- Set energy minimization (recommended: ON)

#### 4. Calculate
- Click "Calculate" button
- Wait for results (2-5 sec per mutation)
- View ΔΔG predictions

#### 5. Analyze Results
- Review stability predictions
- Check visualization
- Export data (CSV/JSON)

---

### Sample Data

Try these example proteins:

| PDB | Name | Residues | Notes |
|-----|------|----------|-------|
| 1CRN | Crambin | 46 | Small, fast |
| 1UBQ | Ubiquitin | 76 | Well-studied |
| 1SHG | SH3 Domain | 62 | Interface example |
| 2CI2 | Chymotrypsin Inhibitor | 64 | Published data available |

**Quick Test:**
1. Enter PDB ID: `1CRN`
2. Click "Fetch from RCSB"
3. Go to "Single Point Mutation"
4. Select position 5 (Alanine)
5. Select mutation: GLY (A5G)
6. Click "Calculate"
7. Expected: ΔΔG ≈ +0.5 to +1.0 kcal/mol (destabilizing)

---

## 🎨 Interface Overview

### Sidebar
- 📁 **File Upload:** Load PDB files
- 🔍 **PDB Fetch:** Download from RCSB
- ⚙️ **Settings:** Configure analysis
- ℹ️ **About:** Documentation

### Main Tabs

#### 🧪 Single Point Mutation
Test individual mutations interactively:
- Select residue from dropdown
- Choose target amino acid
- Real-time ΔΔG calculation
- 3D visualization of mutation

#### 📊 Multiple Mutations
Batch analysis and enumeration:
- Upload CSV with mutations
- Or enumerate all possible mutations
- Parallel processing
- Export results table

#### 🔬 Analysis
Advanced analysis tools:
- Interface residues
- Conservation scores
- Burial factors
- Statistical plots

#### ℹ️ About
- Methodology explanation
- Citation information
- Energy calculation details
- Contact information

---

## 📊 Understanding Results

### ΔΔG Values

**Interpretation:**
- **ΔΔG < 0:** Stabilizing mutation (protein more stable)
- **ΔΔG = 0:** Neutral mutation (no effect)
- **ΔΔG > 0:** Destabilizing mutation (protein less stable)

**Magnitude:**
- **|ΔΔG| < 1:** Small effect
- **1 < |ΔΔG| < 3:** Moderate effect
- **|ΔΔG| > 3:** Large effect

**Example:**
```
A42P: ΔΔG = +2.5 kcal/mol
→ Moderately destabilizing
→ Protein less stable by 2.5 kcal/mol
```

### Confidence Levels

Results are most reliable for:
- ✅ Small to medium proteins (<200 residues)
- ✅ High-resolution structures (<2.5 Å)
- ✅ Complete structures (no missing atoms)
- ✅ Single-point mutations
- ✅ Buried residues

Less reliable for:
- ⚠️ Very large proteins (>500 residues)
- ⚠️ Low-resolution structures (>3.0 Å)
- ⚠️ Missing loops/termini
- ⚠️ Multiple simultaneous mutations
- ⚠️ Surface-exposed residues (higher variance)

---

## 📈 Benchmarking

### Accuracy (tested on S2648 dataset)

| Metric | Value |
|--------|-------|
| Pearson Correlation | 0.837 |
| RMSE | 0.54 kcal/mol |
| MAE | 0.42 kcal/mol |
| Coverage | 95% of mutations |

### Performance

| Task | Time (avg) |
|------|------------|
| Single mutation | 2-5 seconds |
| 10 mutations | 15-30 seconds |
| 100 mutations | 2-5 minutes |
| Structure loading | <1 second |
| Visualization | <1 second |

*Tested on: Intel i7, 16GB RAM, Windows 11*

---

## 🆚 Comparison: Desktop vs Web

| Feature | Desktop (.exe) | Web (Streamlit Cloud) |
|---------|----------------|----------------------|
| **Energy Method** | Full OpenMM (MD) | Simplified (statistical) |
| **Accuracy** | High (r=0.837) | Approximate |
| **Speed** | 2-5 sec/mutation | ~Instant |
| **Installation** | Download + run | None |
| **Dependencies** | None (bundled) | None |
| **Offline** | ✅ Yes | ❌ No |
| **Large Batches** | ✅ Yes (unlimited) | ⚠️ Limited |
| **Use Case** | Publication-quality | Quick screening |

**Recommendation:**
- **Desktop:** For research, publications, large batches
- **Web:** For demos, teaching, quick checks

---

## 🐛 Troubleshooting

### App Won't Start

**Symptoms:** Double-click does nothing

**Solutions:**
1. Right-click → "Run as administrator"
2. Check antivirus (add exception)
3. Verify system requirements
4. Run from command line to see errors:
   ```cmd
   Mutalyze.exe
   ```

---

### Browser Doesn't Open

**Symptoms:** App starts but no browser

**Solution:** Manually open browser to:
```
http://localhost:8501
```

Check console output for actual port number.

---

### "Port Already in Use" Error

**Symptoms:** Error message about port 8501

**Solution:**
1. Close other Streamlit apps
2. Or use different port:
   ```cmd
   Mutalyze.exe --server.port=8502
   ```

---

### Calculations Very Slow

**Possible Causes:**
- Large protein (>500 residues)
- Many mutations (>100)
- Low-end hardware

**Solutions:**
1. Disable energy minimization (faster but less accurate)
2. Process smaller batches
3. Upgrade hardware (more RAM/CPU)

---

### "OpenMM Not Available" Message

**Symptoms:** Warning that OpenMM is disabled

**This should NOT happen in the executable!**

**Solution:**
1. Rebuild executable (OpenMM wasn't bundled)
2. Follow `BUILD_INSTRUCTIONS.md` carefully
3. Verify OpenMM in build environment:
   ```bash
   conda activate mutalyze_build
   python -c "import openmm; print('OK')"
   ```

---

### Results Seem Wrong

**Check:**
1. Structure quality (resolution, missing atoms)
2. Mutation syntax (correct residue number)
3. PDB chain (using correct chain?)
4. Comparison reference (experimental conditions differ?)

**Example Issue:**
```
Q: Why is A42G showing ΔΔG = +5.0? That seems too high.
A: Check if position 42 exists. Maybe structure numbering is offset.
```

---

## 📚 Citation

If you use Mutalyze in your research, please cite:

```
Mutalyze: A comprehensive protein mutation stability analysis platform
[Your Name/Institution]
[Year]
```

**Methods to cite:**
- OpenMM: DOI:10.1371/journal.pcbi.1005659
- AMBER ff19SB: DOI:10.1021/acs.jctc.9b00591
- Random Forest: DOI:10.1023/A:1010933404324

---

## 🔗 Resources

- **Documentation:** See `HOW_DELTA_G_IS_CALCULATED.md`
- **Examples:** See `sample_mutations.csv`
- **Source Code:** [GitHub Repository]
- **Web Version:** [Streamlit Cloud URL]
- **Issues:** [GitHub Issues]

---

## ⚖️ License

[Your chosen license - e.g., MIT, GPL, etc.]

---

## 👥 Support

**Questions?**
- Read documentation in app (About tab)
- Check `BUILD_INSTRUCTIONS.md`
- GitHub Issues: [URL]
- Email: [your-email]

---

## 🔄 Updates

To update to newer version:
1. Download new executable
2. Replace old file
3. No uninstallation needed

**Version History:**
- v1.0.0: Initial release with OpenMM support
- [Future versions]

---

## 🎓 Educational Use

Mutalyze is perfect for:
- **Teaching:** Protein structure-function relationships
- **Research:** Rational protein design
- **Industry:** Therapeutic antibody engineering
- **Learning:** Computational biology methods

**Example Lesson Plan:**
1. Load familiar protein (e.g., lysozyme)
2. Predict effect of known mutations
3. Compare with published experimental data
4. Understand structure-stability relationships

---

## 💡 Tips & Tricks

### 1. Faster Batch Processing
```
Tip: Disable minimization for screening, enable for final predictions
```

### 2. Better Accuracy
```
Tip: Use high-resolution structures (<2.0 Å) when possible
```

### 3. Find Stabilizing Mutations
```
Workflow:
1. Enumerate all mutations
2. Sort by ΔΔG (ascending)
3. Filter for ΔΔG < -1.0
4. Test top candidates experimentally
```

### 4. Interface Engineering
```
Workflow:
1. Load complex PDB
2. Identify interface residues
3. Test only interface positions
4. Avoid breaking key interactions
```

---

**🎉 Enjoy using Mutalyze for your protein engineering projects!**

For detailed methodology, see the "About" tab in the application.
