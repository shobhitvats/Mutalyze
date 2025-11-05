#!/usr/bin/env python3
"""
Validate calibration on diverse mutations from different proteins
Tests generalizability beyond 1CRN training set
"""

from core.pdb_utils import fetch_pdb, save_structure
from core.mutation_builder import MutationBuilder
from core.energy_calc import EnergyCalculator
import tempfile
import numpy as np
from scipy.stats import pearsonr

print("="*70)
print("GENERALIZATION VALIDATION")
print("Testing on mutations NOT used in calibration")
print("="*70)

# Diverse test set from ProThermDB (different proteins, mutation types)
test_mutations = [
    # 1CRN - different positions from training (position 25)
    ('1crn', 'A', 16, 'ILE', 'ALA', 1.8),   # Hydrophobic deletion
    ('1crn', 'A', 16, 'ILE', 'VAL', 1.2),   # Conservative hydrophobic
    ('1crn', 'A', 7, 'THR', 'SER', 0.3),    # Polar to polar
    ('1crn', 'A', 7, 'THR', 'ALA', 1.5),    # Polar to hydrophobic
    
    # 1UBQ - Ubiquitin (different protein!)
    ('1ubq', 'A', 23, 'LEU', 'ALA', 2.1),   # Hydrophobic deletion
    ('1ubq', 'A', 8, 'THR', 'ALA', 1.3),    # Polar to hydrophobic
    ('1ubq', 'A', 56, 'ILE', 'VAL', 0.8),   # Conservative
    
    # 2CI2 - Chymotrypsin inhibitor
    ('2ci2', 'I', 16, 'ALA', 'GLY', 0.9),   # Small to smaller
    ('2ci2', 'I', 29, 'ILE', 'VAL', 0.7),   # Conservative
]

builder = MutationBuilder()
calc = EnergyCalculator()

results = []
print(f"\nTesting {len(test_mutations)} diverse mutations...\n")

for pdb_id, chain, pos, wt_aa, mut_aa, exp_ddg in test_mutations:
    mutation_str = f"{pdb_id.upper()}:{chain}:{pos}:{wt_aa}→{mut_aa}"
    
    try:
        # Fetch structure
        print(f"[{len(results)+1}/{len(test_mutations)}] {mutation_str}")
        structure = fetch_pdb(pdb_id)
        
        # Save WT
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            save_structure(structure, f.name)
            wt_pdb = f.name
        
        # Build mutant
        mutant = builder.apply_mutation(structure.copy(), chain, pos, mut_aa)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            save_structure(mutant, f.name)
            mut_pdb = f.name
        
        # Calculate ΔΔG
        ddg_pred = calc.calculate_ddg(wt_pdb, mut_pdb, minimize=True, mutant_residue=mut_aa)
        error = abs(ddg_pred - exp_ddg)
        
        print(f"  Exp: {exp_ddg:+.2f}, Pred: {ddg_pred:+.2f}, Error: {error:.2f} kcal/mol")
        
        results.append({
            'mutation': mutation_str,
            'exp': exp_ddg,
            'pred': ddg_pred,
            'error': error
        })
        
    except Exception as e:
        print(f"  ERROR: {e}")
        continue

print("\n" + "="*70)
print("RESULTS SUMMARY")
print("="*70)

if len(results) >= 3:
    exp_values = [r['exp'] for r in results]
    pred_values = [r['pred'] for r in results]
    errors = [r['error'] for r in results]
    
    # Calculate metrics
    r, p_value = pearsonr(exp_values, pred_values)
    rmse = np.sqrt(np.mean([(e-p)**2 for e, p in zip(exp_values, pred_values)]))
    mae = np.mean(errors)
    max_error = max(errors)
    
    print(f"\nPerformance on {len(results)} diverse mutations:")
    print(f"  Pearson r:    {r:.3f} (p={p_value:.3f})")
    print(f"  RMSE:         {rmse:.2f} kcal/mol")
    print(f"  MAE:          {mae:.2f} kcal/mol")
    print(f"  Max Error:    {max_error:.2f} kcal/mol")
    
    print(f"\n🎯 Target (FoldX): r=0.50-0.65, RMSE=1.5-2.0")
    
    # Status
    if r >= 0.50 and rmse <= 2.0:
        print("\n✅ COMPETITIVE PERFORMANCE CONFIRMED!")
        print("   Calibration generalizes across proteins and mutation types!")
    elif r >= 0.40:
        print("\n⚠️  Good but needs improvement")
        print("   May need protein-specific calibration or more training data")
    else:
        print("\n❌ Poor generalization")
        print("   Calibration may be overfitted to 1CRN")
    
    # Per-protein breakdown
    print("\nPer-protein breakdown:")
    proteins = {}
    for r in results:
        pdb = r['mutation'].split(':')[0]
        if pdb not in proteins:
            proteins[pdb] = []
        proteins[pdb].append(r['error'])
    
    for pdb, errors in proteins.items():
        avg_error = np.mean(errors)
        print(f"  {pdb.upper()}: {len(errors)} mutations, avg error = {avg_error:.2f} kcal/mol")

else:
    print(f"\n⚠️  Only {len(results)} successful predictions - need more data")

print("\n" + "="*70)
