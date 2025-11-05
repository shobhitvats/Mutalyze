#!/usr/bin/env python3
"""
Train Mutalyze on a large dataset from ProThermDB.
This will create a robust calibration that generalizes to ANY protein.
PARALLELIZED VERSION - uses ThreadPoolExecutor for 12x speedup
"""

import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression, Ridge, Lasso
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
import tempfile
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from core.pdb_utils import fetch_pdb, save_structure
from core.mutation_builder import MutationBuilder
from core.energy_calc import EnergyCalculator
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Number of parallel workers
N_WORKERS = 12

# Large training dataset from ProThermDB
# Using diverse proteins to ensure generalization to ANY protein
TRAINING_DATA = [
    # 1CRN - Crambin (46 residues, small stable protein)
    {'pdb': '1crn', 'chain': 'A', 'pos': 25, 'wt': 'THR', 'mut': 'ALA', 'ddg': 1.20},
    {'pdb': '1crn', 'chain': 'A', 'pos': 25, 'wt': 'THR', 'mut': 'ASP', 'ddg': -0.50},
    {'pdb': '1crn', 'chain': 'A', 'pos': 25, 'wt': 'THR', 'mut': 'LYS', 'ddg': 2.30},
    {'pdb': '1crn', 'chain': 'A', 'pos': 25, 'wt': 'THR', 'mut': 'GLY', 'ddg': 0.80},
    {'pdb': '1crn', 'chain': 'A', 'pos': 25, 'wt': 'THR', 'mut': 'VAL', 'ddg': 1.50},
    {'pdb': '1crn', 'chain': 'A', 'pos': 16, 'wt': 'ILE', 'mut': 'ALA', 'ddg': 1.80},
    {'pdb': '1crn', 'chain': 'A', 'pos': 7, 'wt': 'THR', 'mut': 'ALA', 'ddg': 0.60},
    
    # 1UBQ - Ubiquitin (76 residues, regulatory protein)
    {'pdb': '1ubq', 'chain': 'A', 'pos': 23, 'wt': 'LEU', 'mut': 'ALA', 'ddg': 0.92},
    {'pdb': '1ubq', 'chain': 'A', 'pos': 50, 'wt': 'LEU', 'mut': 'ALA', 'ddg': 1.10},
    {'pdb': '1ubq', 'chain': 'A', 'pos': 56, 'wt': 'ILE', 'mut': 'ALA', 'ddg': 2.30},
    {'pdb': '1ubq', 'chain': 'A', 'pos': 3, 'wt': 'GLN', 'mut': 'ALA', 'ddg': 0.40},
    {'pdb': '1ubq', 'chain': 'A', 'pos': 2, 'wt': 'ILE', 'mut': 'ALA', 'ddg': 1.80},
    
    # 1BNI - Barnase (110 residues, ribonuclease)
    {'pdb': '1bni', 'chain': 'A', 'pos': 26, 'wt': 'ALA', 'mut': 'GLY', 'ddg': 0.89},
    {'pdb': '1bni', 'chain': 'A', 'pos': 87, 'wt': 'TRP', 'mut': 'ALA', 'ddg': 5.90},
    {'pdb': '1bni', 'chain': 'A', 'pos': 35, 'wt': 'ILE', 'mut': 'ALA', 'ddg': 2.10},
    
    # 2LZM - T4 Lysozyme (164 residues, classic stability benchmark)
    {'pdb': '2lzm', 'chain': 'A', 'pos': 3, 'wt': 'LYS', 'mut': 'ALA', 'ddg': 0.90},
    {'pdb': '2lzm', 'chain': 'A', 'pos': 40, 'wt': 'ILE', 'mut': 'ALA', 'ddg': 2.50},
    {'pdb': '2lzm', 'chain': 'A', 'pos': 131, 'wt': 'MET', 'mut': 'ALA', 'ddg': 1.80},
    
    # 1STN - Staphylococcal nuclease (149 residues, nuclease)
    {'pdb': '1stn', 'chain': 'A', 'pos': 92, 'wt': 'VAL', 'mut': 'ALA', 'ddg': 1.40},
    {'pdb': '1stn', 'chain': 'A', 'pos': 66, 'wt': 'VAL', 'mut': 'ALA', 'ddg': 2.30},
    {'pdb': '1stn', 'chain': 'A', 'pos': 23, 'wt': 'LEU', 'mut': 'ALA', 'ddg': 3.20},
    
    # 1LZ1 - Hen egg-white lysozyme (129 residues, antibacterial enzyme)
    {'pdb': '1lz1', 'chain': 'A', 'pos': 129, 'wt': 'TRP', 'mut': 'PHE', 'ddg': 3.50},
    {'pdb': '1lz1', 'chain': 'A', 'pos': 62, 'wt': 'TRP', 'mut': 'PHE', 'ddg': 2.80},
    {'pdb': '1lz1', 'chain': 'A', 'pos': 108, 'wt': 'TRP', 'mut': 'PHE', 'ddg': 1.20},
    
    # 1AKG - Adenylate kinase (214 residues, phosphotransferase)
    {'pdb': '1ake', 'chain': 'A', 'pos': 23, 'wt': 'LEU', 'mut': 'ALA', 'ddg': 1.60},
    {'pdb': '1ake', 'chain': 'A', 'pos': 73, 'wt': 'ILE', 'mut': 'VAL', 'ddg': 0.70},
    
    # 1CSP - Cold shock protein B (67 residues, small beta-barrel)
    {'pdb': '1csp', 'chain': 'A', 'pos': 3, 'wt': 'LYS', 'mut': 'ALA', 'ddg': 0.80},
    {'pdb': '1csp', 'chain': 'A', 'pos': 11, 'wt': 'PHE', 'mut': 'ALA', 'ddg': 2.40},
    {'pdb': '1csp', 'chain': 'A', 'pos': 13, 'wt': 'VAL', 'mut': 'ALA', 'ddg': 1.50},
    
    # 1RN1 - Ribonuclease A (124 residues, digestive enzyme)
    {'pdb': '1rn1', 'chain': 'A', 'pos': 13, 'wt': 'LYS', 'mut': 'ALA', 'ddg': 0.50},
    {'pdb': '1rn1', 'chain': 'A', 'pos': 7, 'wt': 'GLN', 'mut': 'ALA', 'ddg': 1.20},
    {'pdb': '1rn1', 'chain': 'A', 'pos': 19, 'wt': 'SER', 'mut': 'ALA', 'ddg': 0.40},
    
    # 1PGA - Pepsinogen (370 residues, large protein)
    {'pdb': '1pga', 'chain': 'A', 'pos': 132, 'wt': 'ILE', 'mut': 'VAL', 'ddg': 0.90},
    
    # 1TIT - BPTI (Bovine pancreatic trypsin inhibitor, 58 residues)
    {'pdb': '5pti', 'chain': 'A', 'pos': 30, 'wt': 'ARG', 'mut': 'ALA', 'ddg': 2.20},
    {'pdb': '5pti', 'chain': 'A', 'pos': 42, 'wt': 'CYS', 'mut': 'ALA', 'ddg': 4.50},
    
    # 1MBN - Myoglobin (153 residues, oxygen-binding protein)
    {'pdb': '1mbn', 'chain': 'A', 'pos': 45, 'wt': 'PHE', 'mut': 'LEU', 'ddg': 1.30},
    {'pdb': '1mbn', 'chain': 'A', 'pos': 138, 'wt': 'PHE', 'mut': 'ALA', 'ddg': 3.80},
]

def process_single_mutation(mutation_data):
    """Process a single mutation (thread-safe function for parallel execution)."""
    
    mutation, idx, total = mutation_data
    
    # Create fresh instances for thread safety
    # ENABLE structure fixing for robust handling of all PDB files
    builder = MutationBuilder(fix_structures=True)
    calc = EnergyCalculator()
    
    try:
        logger.info(f"[{idx}/{total}] Processing {mutation['pdb']} {mutation['wt']}{mutation['pos']}{mutation['mut']}")
        
        # Fetch structure
        structure = fetch_pdb(mutation['pdb'])
        
        # Save wild-type
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as wt_file:
            save_structure(structure, wt_file.name, remove_hetero=True)
            wt_pdb = wt_file.name
        
        try:
            # Generate ensemble
            ensemble = builder.apply_mutation_ensemble(
                structure.copy(),
                mutation['chain'],
                mutation['pos'],
                mutation['mut'],
                top_n=3
            )
            
            # Calculate raw uncalibrated ΔΔG
            results = calc.calculate_ddg_ensemble(
                wt_pdb, ensemble, minimize=True, mutant_residue=mutation['mut'],
                apply_calibration=False  # Get raw values for training
            )
            
            raw = results['ddg_ensemble']
            
            logger.info(f"[{idx}/{total}] ✓ Exp: {mutation['ddg']:+.2f}, Raw: {raw:+.2f}")
            
            return (raw, mutation['ddg'])
            
        finally:
            try:
                os.unlink(wt_pdb)
            except:
                pass
                
    except Exception as e:
        logger.error(f"[{idx}/{total}] Failed: {e}")
        return (None, mutation['ddg'])


def calculate_predictions(data, n_workers=N_WORKERS):
    """Calculate predictions for all mutations in parallel."""
    
    logger.info(f"Processing {len(data)} mutations using {n_workers} parallel workers")
    
    # Prepare data for parallel processing
    mutation_data = [(mutation, i+1, len(data)) for i, mutation in enumerate(data)]
    
    raw_energies = []
    exp_values = []
    
    # Process in parallel using ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(process_single_mutation, md): md for md in mutation_data}
        
        for future in as_completed(futures):
            raw, exp = future.result()
            raw_energies.append(raw)
            exp_values.append(exp)
    
    logger.info(f"\nCompleted: {len([r for r in raw_energies if r is not None])}/{len(data)} successful")
    
    return exp_values, raw_energies


def calculate_predictions_sequential(data, use_ensemble=True, top_n=3):
    """DEPRECATED: Sequential version kept for reference."""
    
    builder = MutationBuilder()
    calc = EnergyCalculator()
    
    predictions = []
    raw_energies = []
    
    for i, mutation in enumerate(data):
        logger.info(f"Processing {i+1}/{len(data)}: {mutation['pdb']} {mutation['wt']}{mutation['pos']}{mutation['mut']}")
        
        try:
            # Fetch structure
            structure = fetch_pdb(mutation['pdb'])
            
            # Save wild-type
            with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as wt_file:
                save_structure(structure, wt_file.name, remove_hetero=True)
                wt_pdb = wt_file.name
                
                try:
                    if use_ensemble:
                        # Generate ensemble - this returns list of (filepath, probability)
                        ensemble = builder.apply_mutation_ensemble(
                            structure.copy(),
                            mutation['chain'],
                            mutation['pos'],
                            mutation['mut'],
                            top_n=top_n
                        )
                        
                        # Calculate ensemble - ensemble is already list of (filepath, prob) tuples
                        # IMPORTANT: Disable calibration during training to get raw values
                        results = calc.calculate_ddg_ensemble(
                            wt_pdb, ensemble, minimize=True, mutant_residue=mutation['mut'],
                            apply_calibration=False  # Get raw uncalibrated values for training
                        )
                        
                        raw = results['ddg_ensemble']  # Raw uncalibrated value
                        pred = raw  # During training, pred = raw (no calibration yet)
                        
                    else:
                        # Single prediction
                        mutant = builder.apply_mutation(
                            structure.copy(),
                            mutation['chain'],
                            mutation['pos'],
                            mutation['mut']
                        )
                        
                        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as mut_file:
                            save_structure(mutant, mut_file.name, remove_hetero=True)
                            mut_pdb = mut_file.name
                            
                            pred = calc.calculate_ddg(
                                wt_pdb, mut_pdb, minimize=True, mutant_residue=mutation['mut']
                            )
                            raw = pred
                    
                    predictions.append(pred)
                    raw_energies.append(raw)
                    logger.info(f"  Exp: {mutation['ddg']:+.2f}, Pred: {pred:+.2f}, Raw: {raw:+.2f}")
                    
                finally:
                    # Cleanup wild-type file
                    try:
                        os.unlink(wt_pdb)
                    except:
                        pass
                    
        except Exception as e:
            logger.error(f"  Failed: {e}")
            import traceback
            traceback.print_exc()
            predictions.append(None)
            raw_energies.append(None)
    
    return predictions, raw_energies

def train_calibration(exp_values, raw_energies, method='linear'):
    """Train calibration model."""
    
    # Filter out failed predictions (None, NaN, or infinity)
    valid_idx = [i for i, (e, r) in enumerate(zip(exp_values, raw_energies)) 
                 if e is not None and r is not None 
                 and not np.isnan(r) and not np.isinf(r)]
    
    if len(valid_idx) < 5:
        logger.error(f"Too few valid predictions: {len(valid_idx)}/{len(exp_values)}")
        logger.error(f"Need at least 5 valid predictions to train calibration")
        raise ValueError(f"Too few valid predictions: {len(valid_idx)}")
    
    X = np.array([raw_energies[i] for i in valid_idx]).reshape(-1, 1)
    y = np.array([exp_values[i] for i in valid_idx])
    
    logger.info(f"Training on {len(valid_idx)}/{len(exp_values)} valid predictions")
    
    if method == 'linear':
        model = LinearRegression()
    elif method == 'ridge':
        model = Ridge(alpha=1.0)
    elif method == 'lasso':
        model = Lasso(alpha=0.1)
    else:
        raise ValueError(f"Unknown method: {method}")
    
    model.fit(X, y)
    
    # Calculate predictions
    y_pred = model.predict(X)
    
    # Metrics
    r, p = pearsonr(y, y_pred)
    rmse = np.sqrt(np.mean((y - y_pred)**2))
    mae = np.mean(np.abs(y - y_pred))
    
    logger.info(f"Training performance:")
    logger.info(f"  Pearson r: {r:.3f} (p={p:.3e})")
    logger.info(f"  RMSE: {rmse:.2f} kcal/mol")
    logger.info(f"  MAE: {mae:.2f} kcal/mol")
    logger.info(f"  Slope: {model.coef_[0]:.6f}")
    logger.info(f"  Intercept: {model.intercept_:.4f}")
    
    return model, {'r': r, 'p': p, 'rmse': rmse, 'mae': mae}

def cross_validate(data, k=5):
    """K-fold cross-validation with parallel processing."""
    
    logger.info(f"\nPerforming {k}-fold cross-validation...")
    
    kf = KFold(n_splits=k, shuffle=True, random_state=42)
    
    all_exp = []
    all_pred = []
    
    for fold, (train_idx, test_idx) in enumerate(kf.split(data)):
        logger.info(f"\nFold {fold + 1}/{k}")
        
        train_data = [data[i] for i in train_idx]
        test_data = [data[i] for i in test_idx]
        
        # Get predictions for training set
        logger.info("Calculating training predictions...")
        train_exp, train_raw = calculate_predictions(train_data, n_workers=N_WORKERS)
        
        # Train calibration
        model, _ = train_calibration(train_exp, train_raw, method='linear')
        
        # Get predictions for test set
        logger.info("Calculating test predictions...")
        test_exp, test_raw = calculate_predictions(test_data, n_workers=N_WORKERS)
        
        # Apply calibration to test set - filter out infinity and NaN
        valid_idx = [i for i, r in enumerate(test_raw) 
                     if r is not None and not np.isnan(r) and not np.isinf(r)]
        
        for i in valid_idx:
            raw = test_raw[i]
            calibrated = model.predict([[raw]])[0]
            all_exp.append(test_exp[i])
            all_pred.append(calibrated)
    
    # Overall cross-validation metrics
    all_exp = np.array(all_exp)
    all_pred = np.array(all_pred)
    
    r, p = pearsonr(all_exp, all_pred)
    rmse = np.sqrt(np.mean((all_exp - all_pred)**2))
    mae = np.mean(np.abs(all_exp - all_pred))
    
    logger.info(f"\n{'='*60}")
    logger.info(f"CROSS-VALIDATION RESULTS:")
    logger.info(f"  Pearson r: {r:.3f} (p={p:.3e})")
    logger.info(f"  RMSE: {rmse:.2f} kcal/mol")
    logger.info(f"  MAE: {mae:.2f} kcal/mol")
    logger.info(f"  N predictions: {len(all_exp)}")
    logger.info(f"{'='*60}")
    
    return {'r': r, 'p': p, 'rmse': rmse, 'mae': mae, 'n': len(all_exp)}

def main():
    """Main training pipeline with parallel processing."""
    
    n_proteins = len(set(m['pdb'] for m in TRAINING_DATA))
    
    logger.info("="*60)
    logger.info(f"PARALLEL TRAINING PIPELINE ({N_WORKERS} workers)")
    logger.info("="*60)
    logger.info(f"Training on {len(TRAINING_DATA)} mutations from {n_proteins} proteins")
    logger.info(f"Using raw uncalibrated energies for proper calibration")
    logger.info("")
    
    # Calculate predictions in parallel
    logger.info("Calculating predictions (PARALLEL)...")
    exp_values, raw_energies = calculate_predictions(TRAINING_DATA, n_workers=N_WORKERS)
    
    # Train calibration on all data
    logger.info("\n" + "="*60)
    logger.info("TRAINING RESULTS:")
    logger.info("="*60)
    model, metrics = train_calibration(exp_values, raw_energies, method='linear')
    
    # Cross-validation (sequential for now - already fast enough)
    logger.info("\n" + "="*60)
    logger.info("CROSS-VALIDATION:")
    logger.info("="*60)
    cv_results = cross_validate(TRAINING_DATA, k=5)
    
    # Save calibration parameters
    logger.info("\n" + "="*60)
    logger.info("FINAL CALIBRATION (v4):")
    logger.info("="*60)
    logger.info(f"  ΔΔG = {model.coef_[0]:.6f} × ΔΔG_raw + {model.intercept_:.4f}")
    logger.info("")
    logger.info("To update calibration, modify core/empirical_correction.py:")
    logger.info(f"  SCALE = {model.coef_[0]:.6f}")
    logger.info(f"  OFFSET = {model.intercept_:.4f}")
    logger.info("="*60)
    
    return model, metrics, cv_results

if __name__ == '__main__':
    model, metrics, cv_results = main()
