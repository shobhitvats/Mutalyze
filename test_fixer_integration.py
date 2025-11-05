#!/usr/bin/env python3
"""
Quick test to verify structure fixer integration with training pipeline
Tests on a few mutations to ensure it works before full training
"""

import logging
import tempfile
import os
from core.pdb_utils import fetch_pdb, save_structure
from core.mutation_builder import MutationBuilder
from core.energy_calc import EnergyCalculator

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# Test cases - includes problematic PDBs from earlier
TEST_MUTATIONS = [
    # Easy case
    {'pdb': '1crn', 'chain': 'A', 'pos': 25, 'wt': 'THR', 'mut': 'ALA', 'ddg': 1.20},
    # Problematic case (had terminal residue issues)
    {'pdb': '5pti', 'chain': 'A', 'pos': 30, 'wt': 'ARG', 'mut': 'ALA', 'ddg': 2.20},
    # Another problematic case
    {'pdb': '1csp', 'chain': 'A', 'pos': 3, 'wt': 'LYS', 'mut': 'ALA', 'ddg': 0.80},
]

def test_mutation(mutation, builder, calc):
    """Test a single mutation"""
    
    logger.info(f"\nTesting {mutation['pdb']} {mutation['wt']}{mutation['pos']}{mutation['mut']}")
    
    try:
        # Fetch structure
        structure = fetch_pdb(mutation['pdb'])
        logger.info(f"✓ Fetched {mutation['pdb']}")
        
        # Save wild-type
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as wt_file:
            save_structure(structure, wt_file.name, remove_hetero=True)
            wt_pdb = wt_file.name
        
        logger.info(f"✓ Saved wild-type: {wt_pdb}")
        
        try:
            # Generate ensemble (with structure fixing enabled)
            ensemble = builder.apply_mutation_ensemble(
                structure.copy(),
                mutation['chain'],
                mutation['pos'],
                mutation['mut'],
                top_n=3
            )
            
            logger.info(f"✓ Generated {len(ensemble)} rotamers")
            
            # Calculate ΔΔG
            results = calc.calculate_ddg_ensemble(
                wt_pdb, ensemble, minimize=True, mutant_residue=mutation['mut'],
                apply_calibration=False
            )
            
            raw = results['ddg_ensemble']
            logger.info(f"✓ Raw ΔΔG: {raw:+.2f} kcal/mol (Expected: {mutation['ddg']:+.2f})")
            
            return True
            
        finally:
            try:
                os.unlink(wt_pdb)
            except:
                pass
                
    except Exception as e:
        logger.error(f"✗ FAILED: {e}")
        return False


if __name__ == '__main__':
    logger.info("="*60)
    logger.info("Structure Fixer Integration Test")
    logger.info("="*60)
    
    # Create instances WITH structure fixing enabled
    logger.info("\nInitializing with fix_structures=True...")
    builder = MutationBuilder(fix_structures=True)
    calc = EnergyCalculator()
    
    # Test each mutation
    results = []
    for i, mutation in enumerate(TEST_MUTATIONS, 1):
        logger.info(f"\n{'='*60}")
        logger.info(f"Test {i}/{len(TEST_MUTATIONS)}")
        logger.info(f"{'='*60}")
        
        success = test_mutation(mutation, builder, calc)
        results.append((mutation['pdb'], success))
    
    # Summary
    logger.info("\n" + "="*60)
    logger.info("SUMMARY")
    logger.info("="*60)
    
    success_count = sum(1 for _, success in results if success)
    logger.info(f"Success rate: {success_count}/{len(results)} ({100*success_count/len(results):.1f}%)")
    
    for pdb, success in results:
        status = "✓" if success else "✗"
        logger.info(f"{status} {pdb}")
    
    if success_count == len(results):
        logger.info("\n✓ ALL TESTS PASSED! Structure fixer integration working correctly.")
        logger.info("✓ Ready for full training with structure fixing enabled.")
    else:
        logger.warning(f"\n⚠ {len(results) - success_count} test(s) failed. Check logs above.")
