#!/usr/bin/env python3
"""
Re-optimize calibration parameters using expanded training set
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import pearsonr

# Combined training data: original 5 + new 9 = 14 mutations
training_data = [
    # Original 1CRN position 25 (from benchmark)
    (-8.8, 1.20),    # THR25→ALA
    (-72.9, -0.50),  # THR25→ASP
    (-26.7, 2.30),   # THR25→LYS
    (-6.0, 0.80),    # THR25→GLY
    (-16.5, 1.50),   # THR25→VAL
]

# Now let's get raw values for the new mutations
print("Getting raw GBSA values for new mutations...")
print()

from core.pdb_utils import fetch_pdb, save_structure
from core.mutation_builder import MutationBuilder
from core.energy_calc import EnergyCalculator
import tempfile

new_mutations = [
    ('1crn', 'A', 16, 'ILE', 'ALA', 1.80),
    ('1crn', 'A', 16, 'ILE', 'VAL', 1.20),
    ('1crn', 'A', 7, 'THR', 'SER', 0.30),
    ('1crn', 'A', 7, 'THR', 'ALA', 1.50),
    ('1ubq', 'A', 23, 'LEU', 'ALA', 2.10),
    ('1ubq', 'A', 8, 'THR', 'ALA', 1.30),
    ('1ubq', 'A', 56, 'ILE', 'VAL', 0.80),
    ('2ci2', 'I', 16, 'ALA', 'GLY', 0.90),
    ('2ci2', 'I', 29, 'ILE', 'VAL', 0.70),
]

builder = MutationBuilder()
calc = EnergyCalculator()

for pdb_id, chain, pos, wt_aa, mut_aa, exp_ddg in new_mutations:
    try:
        structure = fetch_pdb(pdb_id)
        
        # WT
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            save_structure(structure, f.name)
            wt_pdb = f.name
        
        # Mutant
        mutant = builder.apply_mutation(structure.copy(), chain, pos, mut_aa)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            save_structure(mutant, f.name)
            mut_pdb = f.name
        
        # Get RAW energy (no correction)
        wt_energy = calc.calculate_gbsa_energy(wt_pdb, minimize=True)
        mut_energy = calc.calculate_gbsa_energy(mut_pdb, minimize=True)
        ddg_raw = mut_energy - wt_energy
        
        training_data.append((ddg_raw, exp_ddg))
        print(f"{pdb_id.upper()}:{chain}:{pos}:{wt_aa}→{mut_aa}: raw={ddg_raw:+7.1f}, exp={exp_ddg:+.2f}")
        
    except Exception as e:
        print(f"ERROR on {pdb_id}:{pos}: {e}")
        continue

print(f"\n✅ Collected {len(training_data)} training examples")
print()

# Optimize calibration
raw_values = np.array([x[0] for x in training_data])
exp_values = np.array([x[1] for x in training_data])

def linear(x, a, b):
    return a * x + b

# Fit
params, _ = curve_fit(linear, raw_values, exp_values)
a_opt, b_opt = params

# Predictions
pred_values = linear(raw_values, a_opt, b_opt)
r_opt, _ = pearsonr(exp_values, pred_values)
rmse_opt = np.sqrt(np.mean((exp_values - pred_values)**2))
mae_opt = np.mean(np.abs(exp_values - pred_values))

print("="*70)
print("OPTIMIZED CALIBRATION")
print("="*70)
print(f"\nTraining set: {len(training_data)} mutations")
print(f"\nOptimal parameters:")
print(f"  a (scale)  = {a_opt:.6f}")
print(f"  b (offset) = {b_opt:.4f}")
print(f"\nFormula:")
print(f"  ΔΔG_pred = {a_opt:.6f} × ΔΔG_raw + {b_opt:.4f}")
print(f"\nPerformance:")
print(f"  Pearson r  = {r_opt:.3f}")
print(f"  RMSE       = {rmse_opt:.2f} kcal/mol")
print(f"  MAE        = {mae_opt:.2f} kcal/mol")

# Compare to old calibration
old_a, old_b = 0.0247, 1.71
old_pred = linear(raw_values, old_a, old_b)
old_r, _ = pearsonr(exp_values, old_pred)
old_rmse = np.sqrt(np.mean((exp_values - old_pred)**2))

print(f"\nOld calibration (a={old_a:.4f}, b={old_b:.2f}):")
print(f"  Pearson r  = {old_r:.3f}")
print(f"  RMSE       = {old_rmse:.2f} kcal/mol")

if r_opt > old_r and rmse_opt < old_rmse:
    print(f"\n✅ NEW CALIBRATION IS BETTER!")
    print(f"   Δr = +{r_opt - old_r:.3f}, ΔRMSE = -{old_rmse - rmse_opt:.2f}")
elif r_opt > old_r or rmse_opt < old_rmse:
    print(f"\n⚠️  NEW CALIBRATION IS MIXED")
    print(f"   Δr = {r_opt - old_r:+.3f}, ΔRMSE = {rmse_opt - old_rmse:+.2f}")
else:
    print(f"\n⚠️  OLD CALIBRATION WAS BETTER - KEEP IT")

print("\n" + "="*70)
print("\nTo update calibration in core/empirical_correction.py:")
print(f"    ddg_corrected = {a_opt:.6f} * ddg_raw + {b_opt:.4f}")
print("="*70)
