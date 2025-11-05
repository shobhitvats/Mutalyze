#!/usr/bin/env python3
"""
Comprehensive testing on diverse proteins to validate v5 calibration
Tests on proteins of different sizes, folds, and characteristics
"""

import logging
import tempfile
import os
import time
from pathlib import Path
from core.pdb_utils import fetch_pdb, save_structure
from core.mutation_builder import MutationBuilder
from core.energy_calc import EnergyCalculator
from core.empirical_correction import EmpiricalCorrection

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# Diverse test set spanning different protein characteristics
TEST_PROTEINS = [
    # Small proteins
    {'pdb': '1crn', 'name': 'Crambin', 'size': 46, 'type': 'Small hydrophobic'},
    {'pdb': '2ci2', 'name': 'Chymotrypsin Inhibitor', 'size': 64, 'type': 'Small disulfide-rich'},
    {'pdb': '1csp', 'name': 'Cold Shock Protein', 'size': 67, 'type': 'Beta-barrel'},
    
    # Medium proteins
    {'pdb': '1ubq', 'name': 'Ubiquitin', 'size': 76, 'type': 'Beta-grasp fold'},
    {'pdb': '1bni', 'name': 'Barnase', 'size': 110, 'type': 'Alpha+beta'},
    {'pdb': '1lz1', 'name': 'Lysozyme', 'size': 129, 'type': 'Alpha+beta'},
    
    # Large proteins
    {'pdb': '2lzm', 'name': 'T4 Lysozyme', 'size': 164, 'type': 'Alpha+beta large'},
    {'pdb': '1ake', 'name': 'Adenylate Kinase', 'size': 214, 'type': 'Multi-domain'},
    {'pdb': '1pga', 'name': 'Pepsinogen', 'size': 370, 'type': 'Very large'},
]

# Test mutations (position adjusted per protein)
TEST_MUTATIONS = [
    {'pos_rel': 0.2, 'wt': None, 'mut': 'ALA', 'desc': 'N-terminal region'},
    {'pos_rel': 0.5, 'wt': None, 'mut': 'ALA', 'desc': 'Core region'},
    {'pos_rel': 0.8, 'wt': None, 'mut': 'ALA', 'desc': 'C-terminal region'},
]


def test_single_protein(pdb_info):
    """Test mutations on a single protein"""
    
    pdb_id = pdb_info['pdb']
    logger.info(f"\n{'='*70}")
    logger.info(f"Testing: {pdb_info['name']} ({pdb_id.upper()})")
    logger.info(f"Size: {pdb_info['size']} residues | Type: {pdb_info['type']}")
    logger.info(f"{'='*70}")
    
    try:
        # Fetch structure
        logger.info(f"Fetching {pdb_id}...")
        structure = fetch_pdb(pdb_id)
        
        # Get residue list
        residues = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] == ' ':  # Standard residues only
                        residues.append(residue)
        
        num_residues = len(residues)
        logger.info(f"Found {num_residues} residues")
        
        # Save wild-type
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as wt_file:
            save_structure(structure, wt_file.name, remove_hetero=True)
            wt_pdb = wt_file.name
        
        builder = MutationBuilder(fix_structures=False)
        calc = EnergyCalculator()
        
        results = []
        
        # Test mutations at different positions
        for mut_template in TEST_MUTATIONS:
            # Calculate actual position
            pos = int(num_residues * mut_template['pos_rel'])
            if pos < 1:
                pos = 1
            if pos > num_residues:
                pos = num_residues
            
            target_res = residues[pos - 1]
            wt_aa = target_res.resname
            mut_aa = mut_template['mut']
            
            if wt_aa == mut_aa:
                logger.info(f"Skipping {wt_aa}{pos}{mut_aa} (same residue)")
                continue
            
            mutation_desc = f"{wt_aa}{pos}{mut_aa}"
            logger.info(f"\nTesting {mutation_desc} ({mut_template['desc']})...")
            
            try:
                start_time = time.time()
                
                # Generate ensemble
                ensemble = builder.apply_mutation_ensemble(
                    structure.copy(), 'A', pos, mut_aa, top_n=3
                )
                
                if not ensemble:
                    logger.warning(f"  ✗ Failed to generate ensemble")
                    results.append({
                        'mutation': mutation_desc,
                        'success': False,
                        'error': 'Ensemble generation failed'
                    })
                    continue
                
                # Calculate ΔΔG
                result = calc.calculate_ddg_ensemble(
                    wt_pdb, ensemble, minimize=True,
                    mutant_residue=mut_aa, apply_calibration=True
                )
                
                elapsed = time.time() - start_time
                ddg = result['ddg_ensemble']
                
                if ddg == float('inf'):
                    logger.warning(f"  ✗ Energy calculation failed")
                    results.append({
                        'mutation': mutation_desc,
                        'success': False,
                        'error': 'Energy calculation failed'
                    })
                else:
                    logger.info(f"  ✓ ΔΔG = {ddg:+.2f} kcal/mol ({elapsed:.1f}s)")
                    results.append({
                        'mutation': mutation_desc,
                        'ddg': ddg,
                        'time': elapsed,
                        'success': True
                    })
                
            except Exception as e:
                logger.error(f"  ✗ Error: {e}")
                results.append({
                    'mutation': mutation_desc,
                    'success': False,
                    'error': str(e)
                })
        
        # Cleanup
        try:
            os.unlink(wt_pdb)
        except:
            pass
        
        # Summary for this protein
        success_count = sum(1 for r in results if r.get('success', False))
        total_count = len(results)
        
        logger.info(f"\n{pdb_info['name']} Summary:")
        logger.info(f"  Success: {success_count}/{total_count} ({100*success_count/total_count if total_count > 0 else 0:.0f}%)")
        
        if success_count > 0:
            avg_time = sum(r['time'] for r in results if r.get('success')) / success_count
            logger.info(f"  Avg time: {avg_time:.1f}s per mutation")
        
        return {
            'pdb': pdb_id,
            'name': pdb_info['name'],
            'size': pdb_info['size'],
            'type': pdb_info['type'],
            'results': results,
            'success_rate': success_count / total_count if total_count > 0 else 0
        }
        
    except Exception as e:
        logger.error(f"Failed to test {pdb_id}: {e}")
        return {
            'pdb': pdb_id,
            'name': pdb_info['name'],
            'error': str(e),
            'success_rate': 0
        }


def main():
    logger.info("="*70)
    logger.info("COMPREHENSIVE PROTEIN TESTING - v5 Calibration")
    logger.info("="*70)
    
    # Show model info
    info = EmpiricalCorrection.get_model_info()
    logger.info(f"\nCalibration: {info['model_type']}")
    logger.info(f"Performance: r={info['training_r']:.3f}, RMSE={info['training_rmse']:.2f}")
    logger.info(f"Training: {info['training_samples']} mutations\n")
    
    all_results = []
    
    # Test each protein
    for pdb_info in TEST_PROTEINS:
        result = test_single_protein(pdb_info)
        all_results.append(result)
    
    # Overall summary
    logger.info("\n" + "="*70)
    logger.info("OVERALL SUMMARY")
    logger.info("="*70)
    
    total_proteins = len(all_results)
    successful_proteins = sum(1 for r in all_results if r.get('success_rate', 0) > 0)
    
    logger.info(f"\nProteins tested: {total_proteins}")
    logger.info(f"Successful: {successful_proteins}/{total_proteins} ({100*successful_proteins/total_proteins:.0f}%)")
    
    # Success rate by protein size
    logger.info(f"\n{'Protein':<25} {'Size':<10} {'Success Rate':<15} {'Type':<20}")
    logger.info("-"*70)
    
    for result in all_results:
        if 'error' in result:
            logger.info(f"{result['name']:<25} {result['size']:<10} {'ERROR':<15} {result.get('type', 'N/A'):<20}")
        else:
            success_pct = result['success_rate'] * 100
            logger.info(f"{result['name']:<25} {result['size']:<10} {success_pct:>6.0f}%{'':<9} {result['type']:<20}")
    
    # Overall statistics
    all_success_rates = [r['success_rate'] for r in all_results if 'error' not in r]
    if all_success_rates:
        avg_success = sum(all_success_rates) / len(all_success_rates) * 100
        logger.info(f"\nAverage success rate: {avg_success:.1f}%")
    
    logger.info("\n" + "="*70)
    logger.info("TESTING COMPLETE")
    logger.info("="*70)


if __name__ == '__main__':
    main()
