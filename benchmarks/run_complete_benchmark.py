"""
Complete Benchmark: Mutalyze vs ProThermDB
Run full validation pipeline with real predictions
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import logging
import pandas as pd
import tempfile
import time

from core.pdb_utils import fetch_pdb, save_structure
from core.mutation_builder import MutationBuilder
from core.energy_calc import EnergyCalculator
from benchmarks.protherm_validation import ProThermValidator, create_sample_protherm_data

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def run_benchmark(test_set: list, use_ensemble: bool = False, 
                  top_n_rotamers: int = 3, use_local_energy: bool = False) -> list:
    """
    Run benchmark on test set
    
    Args:
        test_set: List of mutation dictionaries
        use_ensemble: Use rotamer ensemble vs single structure
        top_n_rotamers: Number of rotamers to generate
        use_local_energy: Use local interaction energy (FoldX approach)
        
    Returns:
        List of result dictionaries
    """
    builder = MutationBuilder()
    calc = EnergyCalculator()
    
    results = []
    
    print("\n" + "="*70)
    print(f"RUNNING BENCHMARK ({'ENSEMBLE' if use_ensemble else 'SINGLE'} MODE)")
    print("="*70)
    
    for idx, row in enumerate(test_set):
        pdb_id = row['pdb_id']
        chain = row['chain']
        position = int(row['position'])
        wild_aa = row['wild_aa']
        mut_aa = row['mut_aa']
        ddg_exp = row['ddg_exp']
        
        mutation_str = f"{pdb_id}:{chain}:{position}:{wild_aa}→{mut_aa}"
        print(f"\n[{idx+1}/{len(test_set)}] {mutation_str}")
        print(f"  Experimental ΔΔG: {ddg_exp:+.2f} kcal/mol")
        
        try:
            # Fetch structure
            structure = fetch_pdb(pdb_id)
            
            # Save wildtype to temp file
            wt_file = tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False)
            save_structure(structure, wt_file.name)
            wt_file.close()
            
            start_time = time.time()
            
            if use_ensemble:
                # Ensemble method
                ensemble = builder.apply_mutation_ensemble(
                    structure, chain, position, mut_aa, top_n=top_n_rotamers
                )
                
                if len(ensemble) > 0:
                    ens_results = calc.calculate_ddg_ensemble(
                        wt_file.name, ensemble, minimize=True, mutant_residue=mut_aa
                    )
                    ddg_pred = ens_results['ddg_ensemble']
                    method = f"Ensemble (n={len(ensemble)})"
                else:
                    logger.error(f"  Failed to generate ensemble")
                    ddg_pred = float('nan')
                    method = "Ensemble (FAILED)"
            else:
                # Single structure method
                mutant = builder.apply_mutation(structure, chain, position, mut_aa)
                mut_file = tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False)
                save_structure(mutant, mut_file.name)
                mut_file.close()
                
                ddg_pred = calc.calculate_ddg(
                    wt_file.name, mut_file.name, minimize=True, 
                    mutant_residue=mut_aa, mutation_site=position,
                    use_local_energy=use_local_energy
                )
                method = "Single (Local)" if use_local_energy else "Single"
            
            elapsed = time.time() - start_time
            
            error = abs(ddg_pred - ddg_exp)
            
            print(f"  Predicted ΔΔG:    {ddg_pred:+.2f} kcal/mol ({method})")
            print(f"  Error:            {error:.2f} kcal/mol")
            print(f"  Time:             {elapsed:.1f} seconds")
            
            results.append({
                'pdb_id': pdb_id,
                'chain': chain,
                'position': position,
                'wild_aa': wild_aa,
                'mut_aa': mut_aa,
                'ddg_exp': ddg_exp,
                'ddg_pred': ddg_pred,
                'error': error,
                'method': method,
                'time_seconds': elapsed
            })
            
        except Exception as e:
            logger.error(f"  Failed: {e}")
            results.append({
                'pdb_id': pdb_id,
                'chain': chain,
                'position': position,
                'wild_aa': wild_aa,
                'mut_aa': mut_aa,
                'ddg_exp': ddg_exp,
                'ddg_pred': float('nan'),
                'error': float('nan'),
                'method': 'FAILED',
                'time_seconds': 0
            })
    
    return results  # Return list, not DataFrame


def main():
    """Run complete benchmark pipeline"""
    
    print("\n" + "="*70)
    print("MUTALYZE COMPLETE BENCHMARK")
    print("="*70)
    
    # Create sample test set
    # In production, load from actual ProThermDB file
    print("\n1. Creating test dataset...")
    test_data = create_sample_protherm_data()
    
    # Save to file for validator
    test_file = '/tmp/mutalyze_test_set.csv'
    test_data.to_csv(test_file, index=False)
    print(f"   Test set: {len(test_data)} mutations")
    
    # Run benchmark with single structure method
    print("\n2. Running Single Structure Benchmark...")
    results_single_df = run_benchmark(test_data.to_dict('records'), use_ensemble=False, use_local_energy=False)
    results_single = pd.DataFrame(results_single_df)
    
    # Run benchmark with ensemble method  
    print("\n3. Running Rotamer Ensemble Benchmark...")
    results_ensemble_df = run_benchmark(test_data.to_dict('records'), use_ensemble=True, top_n_rotamers=3, use_local_energy=False)
    results_ensemble = pd.DataFrame(results_ensemble_df)
    
    # Initialize validator
    validator = ProThermValidator(test_file)
    
    # Calculate metrics for single method
    print("\n" + "="*70)
    print("SINGLE STRUCTURE RESULTS")
    print("="*70)
    
    valid_single = results_single[~results_single['ddg_pred'].isna()]
    if len(valid_single) > 0:
        metrics_single = validator.calculate_metrics(
            list(valid_single['ddg_pred']),
            list(valid_single['ddg_exp'])
        )
        validator.print_report(metrics_single)
        validator.plot_correlation(metrics_single, 
                                  '/tmp/mutalyze_single_correlation.png',
                                  'Mutalyze (Single Structure)')
    else:
        print("No valid predictions")
    
    # Calculate metrics for ensemble method
    print("\n" + "="*70)
    print("ROTAMER ENSEMBLE RESULTS")
    print("="*70)
    
    valid_ensemble = results_ensemble[~results_ensemble['ddg_pred'].isna()]
    if len(valid_ensemble) > 0:
        metrics_ensemble = validator.calculate_metrics(
            list(valid_ensemble['ddg_pred']),
            list(valid_ensemble['ddg_exp'])
        )
        validator.print_report(metrics_ensemble)
        validator.plot_correlation(metrics_ensemble, 
                                  '/tmp/mutalyze_ensemble_correlation.png',
                                  'Mutalyze (Rotamer Ensemble)')
    else:
        print("No valid predictions")
    
    # Comparison
    print("\n" + "="*70)
    print("COMPARISON: Single vs Ensemble")
    print("="*70)
    
    if len(valid_single) > 0 and len(valid_ensemble) > 0:
        print(f"\n{'Metric':<20} {'Single':>12} {'Ensemble':>12} {'Δ':>10}")
        print("-" * 56)
        print(f"{'Pearson r':<20} {metrics_single.pearson_r:>12.3f} {metrics_ensemble.pearson_r:>12.3f} {metrics_ensemble.pearson_r - metrics_single.pearson_r:>+10.3f}")
        print(f"{'RMSE (kcal/mol)':<20} {metrics_single.rmse:>12.2f} {metrics_ensemble.rmse:>12.2f} {metrics_ensemble.rmse - metrics_single.rmse:>+10.2f}")
        print(f"{'MAE (kcal/mol)':<20} {metrics_single.mae:>12.2f} {metrics_ensemble.mae:>12.2f} {metrics_ensemble.mae - metrics_single.mae:>+10.2f}")
        print(f"{'Accuracy':<20} {metrics_single.accuracy:>11.1%} {metrics_ensemble.accuracy:>11.1%} {metrics_ensemble.accuracy - metrics_single.accuracy:>+9.1%}")
        
        avg_time_single = results_single['time_seconds'].mean()
        avg_time_ensemble = results_ensemble['time_seconds'].mean()
        print(f"{'Avg Time (sec)':<20} {avg_time_single:>12.1f} {avg_time_ensemble:>12.1f} {avg_time_ensemble - avg_time_single:>+10.1f}")
    
    # Save detailed results
    results_file = '/tmp/mutalyze_benchmark_results.csv'
    combined_results = pd.DataFrame({
        'mutation': [f"{r['pdb_id']}:{r['chain']}:{r['position']}:{r['wild_aa']}→{r['mut_aa']}" 
                     for _, r in results_single.iterrows()],
        'ddg_exp': results_single['ddg_exp'],
        'ddg_single': results_single['ddg_pred'],
        'ddg_ensemble': results_ensemble['ddg_pred'],
        'error_single': results_single['error'],
        'error_ensemble': results_ensemble['error']
    })
    combined_results.to_csv(results_file, index=False)
    print(f"\nDetailed results saved to: {results_file}")
    
    print("\n" + "="*70)
    print("BENCHMARK COMPLETE")
    print("="*70)
    print(f"\nPlots saved:")
    print(f"  - /tmp/mutalyze_single_correlation.png")
    print(f"  - /tmp/mutalyze_ensemble_correlation.png")
    print(f"\nResults: {results_file}")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()
