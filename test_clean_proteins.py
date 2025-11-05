#!/usr/bin/env python3
"""
Test v5 calibration on clean, well-structured proteins only.
Focuses on proteins known to have proper PDB formatting.
"""

import logging
import time
from pathlib import Path
from core.mutation_builder import MutationBuilder
from core.energy_calc import EnergyCalculator
from core.empirical_correction import EmpiricalCorrection
from core.pdb_utils import PDBFetcher

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Clean proteins with proper PDB structure
TEST_PROTEINS = {
    '1crn': [  # 46 residues - crambin
        ('THR2ALA', 'N-terminal'),
        ('ILE7VAL', 'Core hydrophobic'),
    ],
    '1ubq': [  # 76 residues - ubiquitin
        ('MET1ALA', 'N-terminal'),
        ('LEU8ALA', 'Core hydrophobic'),
        ('VAL70ALA', 'C-terminal region'),
    ],
    '1lz1': [  # 129 residues - lysozyme
        ('LYS1ALA', 'N-terminal'),
        ('TRP62ALA', 'Core aromatic'),
    ],
    '2lzm': [  # 164 residues - lysozyme mutant
        ('LYS1ALA', 'N-terminal'),
        ('TRP63ALA', 'Core aromatic'),
    ],
}

def test_protein(pdb_id: str, mutations: list):
    """Test a protein with specified mutations."""
    logger.info(f"\n{'='*70}")
    logger.info(f"Testing: {pdb_id.upper()}")
    logger.info(f"{'='*70}")
    
    pdb_path = Path('data/pdb_cache') / f'{pdb_id}.pdb'
    
    if not pdb_path.exists():
        logger.error(f"PDB file not found: {pdb_path}")
        return 0, len(mutations)
    
    successful = 0
    total = len(mutations)
    times = []
    
    for mutation_str, description in mutations:
        logger.info(f"\nMutation: {mutation_str} ({description})")
        start_time = time.time()
        
        try:
            # Build mutant structure
            builder = MutationBuilder(str(pdb_path))
            mutant_pdb = builder.build_mutation(mutation_str)
            
            if not mutant_pdb:
                logger.error(f"  ✗ Failed to build mutation")
                continue
            
            # Calculate energies
            calculator = EnergyCalculator()
            wt_energy = calculator.calculate_energy(str(pdb_path))
            mut_energy = calculator.calculate_energy(mutant_pdb)
            
            if wt_energy is None or mut_energy is None:
                logger.error(f"  ✗ Energy calculation failed")
                continue
            
            # Calculate raw and calibrated ΔΔG
            ddg_raw = mut_energy - wt_energy
            ddg_calibrated = EmpiricalCorrection.apply_correction(ddg_raw)
            
            elapsed = time.time() - start_time
            times.append(elapsed)
            
            logger.info(f"  ✓ WT Energy: {wt_energy:.2f} kcal/mol")
            logger.info(f"  ✓ Mut Energy: {mut_energy:.2f} kcal/mol")
            logger.info(f"  ✓ Raw ΔΔG: {ddg_raw:.2f} kcal/mol")
            logger.info(f"  ✓ Calibrated ΔΔG: {ddg_calibrated:.2f} kcal/mol")
            logger.info(f"  ✓ Time: {elapsed:.1f}s")
            successful += 1
        
        except Exception as e:
            elapsed = time.time() - start_time
            logger.error(f"  ✗ Exception: {str(e)}")
            logger.info(f"  ✗ Time: {elapsed:.1f}s")
    
    # Summary
    success_rate = (successful / total * 100) if total > 0 else 0
    avg_time = sum(times) / len(times) if times else 0
    
    logger.info(f"\n{pdb_id.upper()} Summary:")
    logger.info(f"  Success: {successful}/{total} ({success_rate:.0f}%)")
    if times:
        logger.info(f"  Avg time: {avg_time:.1f}s (min: {min(times):.1f}s, max: {max(times):.1f}s)")
    
    return successful, total

def main():
    """Run tests on all clean proteins."""
    logger.info("="*70)
    logger.info("CLEAN PROTEIN TEST SUITE (v5 Calibration)")
    logger.info("="*70)
    logger.info(f"Testing {len(TEST_PROTEINS)} proteins with proper PDB structure")
    logger.info("")
    
    total_successful = 0
    total_tests = 0
    start_time = time.time()
    
    for pdb_id, mutations in TEST_PROTEINS.items():
        successful, total = test_protein(pdb_id, mutations)
        total_successful += successful
        total_tests += total
    
    # Overall summary
    elapsed = time.time() - start_time
    overall_rate = (total_successful / total_tests * 100) if total_tests > 0 else 0
    
    logger.info(f"\n{'='*70}")
    logger.info("OVERALL RESULTS")
    logger.info(f"{'='*70}")
    logger.info(f"Total Success: {total_successful}/{total_tests} ({overall_rate:.0f}%)")
    logger.info(f"Total Time: {elapsed:.1f}s")
    logger.info(f"Avg per mutation: {elapsed/total_tests:.1f}s")
    logger.info(f"{'='*70}")
    
    return 0 if overall_rate >= 90 else 1

if __name__ == '__main__':
    exit(main())
