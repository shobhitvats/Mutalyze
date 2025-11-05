#!/usr/bin/env python3
"""
Batch Prediction Tool
Efficiently process multiple mutations and export results
"""

import pandas as pd
import json
from pathlib import Path
from core.pdb_utils import fetch_pdb, save_structure
from core.mutation_builder import MutationBuilder
from core.energy_calc import EnergyCalculator
import tempfile
from datetime import datetime

class MutalyzeBatch:
    """Batch processing for multiple mutations"""
    
    def __init__(self):
        self.builder = MutationBuilder()
        self.calc = EnergyCalculator()
        
    def predict_from_csv(self, input_csv: str, use_ensemble: bool = True, 
                        top_n: int = 3) -> pd.DataFrame:
        """
        Process mutations from CSV file
        
        CSV format:
        pdb_id,chain,position,wild_aa,mut_aa,exp_ddg
        1crn,A,25,THR,ASP,-0.50
        
        Args:
            input_csv: Path to input CSV
            use_ensemble: Use ensemble method (slower but more accurate)
            top_n: Number of rotamers if using ensemble
            
        Returns:
            DataFrame with predictions
        """
        print(f"Loading mutations from {input_csv}...")
        df = pd.read_csv(input_csv)
        
        required_cols = ['pdb_id', 'chain', 'position', 'wild_aa', 'mut_aa']
        for col in required_cols:
            if col not in df.columns:
                raise ValueError(f"Missing required column: {col}")
        
        results = []
        total = len(df)
        
        print(f"\nProcessing {total} mutations...")
        print(f"Method: {'Ensemble' if use_ensemble else 'Single'}")
        print(f"{'='*70}\n")
        
        for idx, row in df.iterrows():
            pdb_id = row['pdb_id']
            chain = row['chain']
            position = int(row['position'])
            wild_aa = row['wild_aa']
            mut_aa = row['mut_aa']
            exp_ddg = row.get('exp_ddg', None)
            
            mutation_str = f"{pdb_id.upper()}:{chain}:{position}:{wild_aa}→{mut_aa}"
            print(f"[{idx+1}/{total}] {mutation_str}...", end=' ')
            
            try:
                # Fetch structure
                structure = fetch_pdb(pdb_id)
                
                # Save wildtype
                with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
                    save_structure(structure, f.name)
                    wt_pdb = f.name
                
                if use_ensemble:
                    # Ensemble prediction
                    mutant_ensemble = self.builder.apply_mutation_ensemble(
                        structure.copy(), chain, position, mut_aa, top_n=top_n
                    )
                    
                    ens_results = self.calc.calculate_ddg_ensemble(
                        wt_pdb, mutant_ensemble, minimize=True, mutant_residue=mut_aa
                    )
                    
                    ddg_pred = ens_results['ddg_ensemble']
                    ddg_std = ens_results['ddg_std']
                    ddg_ci = ens_results['ddg_uncertainty']
                    
                else:
                    # Single structure prediction
                    mutant = self.builder.apply_mutation(structure.copy(), chain, position, mut_aa)
                    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
                        save_structure(mutant, f.name)
                        mut_pdb = f.name
                    
                    ddg_pred = self.calc.calculate_ddg(wt_pdb, mut_pdb, minimize=True, 
                                                      mutant_residue=mut_aa)
                    ddg_std = 0.0
                    ddg_ci = 0.0
                
                # Calculate error if experimental value provided
                error = abs(ddg_pred - exp_ddg) if exp_ddg is not None else None
                
                result = {
                    'mutation': mutation_str,
                    'pdb_id': pdb_id,
                    'chain': chain,
                    'position': position,
                    'wild_aa': wild_aa,
                    'mut_aa': mut_aa,
                    'ddg_predicted': ddg_pred,
                    'ddg_std': ddg_std,
                    'ci_95_lower': ddg_pred - ddg_ci,
                    'ci_95_upper': ddg_pred + ddg_ci,
                    'exp_ddg': exp_ddg,
                    'error': error,
                    'method': 'Ensemble' if use_ensemble else 'Single'
                }
                
                results.append(result)
                print(f"✓ ΔΔG = {ddg_pred:+.2f} ± {ddg_std:.2f} kcal/mol")
                
            except Exception as e:
                print(f"✗ ERROR: {e}")
                result = {
                    'mutation': mutation_str,
                    'pdb_id': pdb_id,
                    'chain': chain,
                    'position': position,
                    'wild_aa': wild_aa,
                    'mut_aa': mut_aa,
                    'ddg_predicted': None,
                    'ddg_std': None,
                    'ci_95_lower': None,
                    'ci_95_upper': None,
                    'exp_ddg': exp_ddg,
                    'error': None,
                    'method': 'Ensemble' if use_ensemble else 'Single'
                }
                results.append(result)
        
        return pd.DataFrame(results)
    
    def export_results(self, df: pd.DataFrame, output_prefix: str):
        """
        Export results in multiple formats
        
        Args:
            df: Results DataFrame
            output_prefix: Prefix for output files (e.g., 'results')
        """
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        
        # CSV export
        csv_file = f"{output_prefix}_{timestamp}.csv"
        df.to_csv(csv_file, index=False)
        print(f"\n✓ CSV saved: {csv_file}")
        
        # JSON export
        json_file = f"{output_prefix}_{timestamp}.json"
        df.to_json(json_file, orient='records', indent=2)
        print(f"✓ JSON saved: {json_file}")
        
        # Summary report
        report_file = f"{output_prefix}_{timestamp}_summary.txt"
        with open(report_file, 'w') as f:
            f.write("="*70 + "\n")
            f.write("MUTALYZE BATCH PREDICTION SUMMARY\n")
            f.write("="*70 + "\n\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Total mutations: {len(df)}\n")
            f.write(f"Successful: {df['ddg_predicted'].notna().sum()}\n")
            f.write(f"Failed: {df['ddg_predicted'].isna().sum()}\n")
            f.write(f"Method: {df['method'].iloc[0] if len(df) > 0 else 'N/A'}\n\n")
            
            if 'exp_ddg' in df.columns and df['exp_ddg'].notna().any():
                valid = df[df['error'].notna()]
                if len(valid) > 0:
                    import numpy as np
                    from scipy.stats import pearsonr
                    
                    exp = valid['exp_ddg'].values
                    pred = valid['ddg_predicted'].values
                    
                    r, p = pearsonr(exp, pred)
                    rmse = np.sqrt(np.mean((exp - pred)**2))
                    mae = np.mean(np.abs(exp - pred))
                    
                    f.write("PERFORMANCE METRICS\n")
                    f.write("-"*70 + "\n")
                    f.write(f"Pearson r:  {r:.3f} (p={p:.3f})\n")
                    f.write(f"RMSE:       {rmse:.2f} kcal/mol\n")
                    f.write(f"MAE:        {mae:.2f} kcal/mol\n")
                    f.write(f"Max error:  {valid['error'].max():.2f} kcal/mol\n\n")
            
            f.write("PREDICTIONS\n")
            f.write("-"*70 + "\n")
            for _, row in df.iterrows():
                f.write(f"{row['mutation']:<30} ")
                if pd.notna(row['ddg_predicted']):
                    f.write(f"ΔΔG = {row['ddg_predicted']:+.2f} ± {row['ddg_std']:.2f} kcal/mol")
                    if pd.notna(row['exp_ddg']):
                        f.write(f" (exp: {row['exp_ddg']:+.2f}, err: {row['error']:.2f})")
                else:
                    f.write("FAILED")
                f.write("\n")
        
        print(f"✓ Summary saved: {report_file}")
        print(f"\n{'='*70}")


if __name__ == "__main__":
    # Demo: Create sample input and run batch prediction
    print("="*70)
    print("MUTALYZE BATCH PREDICTION DEMO")
    print("="*70)
    print()
    
    # Create sample input CSV
    sample_mutations = pd.DataFrame([
        {'pdb_id': '1crn', 'chain': 'A', 'position': 25, 'wild_aa': 'THR', 'mut_aa': 'ASP', 'exp_ddg': -0.50},
        {'pdb_id': '1crn', 'chain': 'A', 'position': 16, 'wild_aa': 'ILE', 'mut_aa': 'ALA', 'exp_ddg': 1.80},
        {'pdb_id': '1crn', 'chain': 'A', 'position': 25, 'wild_aa': 'THR', 'mut_aa': 'ALA', 'exp_ddg': 1.20},
        {'pdb_id': '1crn', 'chain': 'A', 'position': 16, 'wild_aa': 'ILE', 'mut_aa': 'VAL', 'exp_ddg': 1.20},
    ])
    
    input_file = 'sample_mutations.csv'
    sample_mutations.to_csv(input_file, index=False)
    print(f"✓ Created sample input: {input_file}")
    print()
    
    # Run batch prediction
    batch = MutalyzeBatch()
    results = batch.predict_from_csv(input_file, use_ensemble=True, top_n=3)
    
    # Export results
    batch.export_results(results, 'batch_results')
    
    print(f"\n{'='*70}")
    print("BATCH PREDICTION COMPLETE")
    print("="*70)
    print("\n✨ Features demonstrated:")
    print("  • CSV input/output")
    print("  • Batch processing")
    print("  • Ensemble predictions with confidence intervals")
    print("  • Automatic performance metrics")
    print("  • Multiple export formats (CSV, JSON, TXT)")
    print("\n🎯 Use this for high-throughput screening!")
    print("="*70)
