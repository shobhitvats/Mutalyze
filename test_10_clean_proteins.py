#!/usr/bin/env python3
"""
Test v5 calibration on 10 diverse, clean proteins.
These proteins are selected for:
1. Clean PDB structure (no terminal residue issues)
2. Size diversity (small to large)
3. Different protein families
"""

import logging
import time
import sys
import tempfile
import os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio.PDB import PDBParser, PDBIO, Select
from core.mutation_builder import MutationBuilder
from core.energy_calc import calculate_stability_change
from core.empirical_correction import EmpiricalCorrection
from core.pdb_utils import save_structure

# Use all available cores
MAX_WORKERS = 12  # Your laptop has 12 threads

class ProteinOnlySelect(Select):
    """Select only protein residues (exclude water, ions, ligands)."""
    def accept_residue(self, residue):
        # Only accept standard amino acid residues (hetero flag is ' ')
        return residue.id[0] == ' '

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# 10 diverse clean proteins with mutations
TEST_PROTEINS = {
    '1crn': {  # 46 residues - crambin (very small)
        'name': 'Crambin',
        'mutations': [
            ('THR2ALA', 'N-terminal hydrophilic→hydrophobic'),
            ('ILE7VAL', 'Core hydrophobic (size reduction)'),
        ]
    },
    '1ubq': {  # 76 residues - ubiquitin (small)
        'name': 'Ubiquitin',
        'mutations': [
            ('LEU8ALA', 'Core hydrophobic→small'),
            ('ILE44ALA', 'Binding interface'),
        ]
    },
    '1lz1': {  # 129 residues - lysozyme (medium)
        'name': 'Lysozyme',
        'mutations': [
            ('TRP62ALA', 'Active site aromatic'),
            ('ASP52ALA', 'Active site acidic'),
        ]
    },
    '2lzm': {  # 164 residues - lysozyme mutant (medium)
        'name': 'Lysozyme T4',
        'mutations': [
            ('TRP63ALA', 'Core aromatic'),
            ('LEU99ALA', 'Core hydrophobic'),
        ]
    },
    '1stn': {  # 68 residues - ovomucoid (small)
        'name': 'Ovomucoid',
        'mutations': [
            ('LEU8ALA', 'Core hydrophobic'),
            ('PHE30ALA', 'Core aromatic'),
        ]
    },
    '1aps': {  # 98 residues - protease (small-medium)
        'name': 'Protease',
        'mutations': [
            ('TRP26ALA', 'Active site'),
            ('LEU40ALA', 'Core hydrophobic'),
        ]
    },
    '1mbn': {  # 153 residues - myoglobin (medium)
        'name': 'Myoglobin',
        'mutations': [
            ('LEU29ALA', 'Core hydrophobic'),
            ('PHE43ALA', 'Heme contact'),
        ]
    },
    '1lmb': {  # 165 residues - lambda repressor (medium)
        'name': 'Lambda Repressor',
        'mutations': [
            ('LEU18ALA', 'Core hydrophobic'),
            ('TRP37ALA', 'DNA binding'),
        ]
    },
    '1ris': {  # 97 residues - ribonuclease (small-medium)
        'name': 'Ribonuclease',
        'mutations': [
            ('VAL43ALA', 'Core hydrophobic'),
            ('PHE120ALA', 'Active site'),
        ]
    },
    '1vqb': {  # 123 residues - VQB protein (medium)
        'name': 'VQB Protein',
        'mutations': [
            ('LEU30ALA', 'Core hydrophobic'),
            ('VAL50ALA', 'Core hydrophobic'),
        ]
    },
}

def test_mutation(pdb_path: str, mutation_str: str, description: str) -> dict:
    """Test a single mutation and return results."""
    start_time = time.time()
    result = {
        'mutation': mutation_str,
        'description': description,
        'success': False,
        'time': 0,
        'error': None,
        'ddg_raw': None,
        'ddg_calibrated': None,
        'wt_energy': None,
        'mut_energy': None
    }
    
    try:
        # Parse structure and clean it (remove water, ions, ligands)
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_path)
        
        # Save cleaned structure (protein only) to use as wildtype
        with tempfile.NamedTemporaryFile(mode='w', suffix='_clean.pdb', delete=False) as tmp:
            io = PDBIO()
            io.set_structure(structure)
            io.save(tmp.name, ProteinOnlySelect())
            clean_pdb_path = tmp.name
        
        # Re-parse the cleaned structure
        structure = parser.get_structure('protein', clean_pdb_path)
        
        # Get the first chain (some PDBs use numeric chain IDs)
        chain = list(structure[0])[0]
        chain_id = chain.id
        
        # Parse mutation
        residue_num = int(mutation_str[3:-3])
        wt_res = mutation_str[:3]
        mut_res = mutation_str[-3:]
        
        # Build mutant
        builder = MutationBuilder(fix_structures=True)
        mutant_structure = builder.apply_mutation(
            structure, 
            chain_id,
            residue_num,
            mut_res
        )
        
        # Save mutant to temp file (also protein only)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
            io = PDBIO()
            io.set_structure(mutant_structure)
            io.save(tmp.name, ProteinOnlySelect())
            temp_mutant = tmp.name
        
        # Calculate ΔΔG using cleaned structures
        ddg_raw = calculate_stability_change(clean_pdb_path, temp_mutant, method='gbsa')
        
        # Clean up
        os.unlink(temp_mutant)
        os.unlink(clean_pdb_path)
        
        if ddg_raw is None or ddg_raw == float('inf') or ddg_raw == float('-inf'):
            result['error'] = 'Energy calculation returned invalid value'
            result['error'] = 'Energy calculation returned invalid value'
            return result
        
        # Apply v5 calibration
        ddg_calibrated = EmpiricalCorrection.apply_correction(ddg_raw)
        
        result['success'] = True
        result['ddg_raw'] = ddg_raw
        result['ddg_calibrated'] = ddg_calibrated
        result['time'] = time.time() - start_time
        
    except Exception as e:
        result['error'] = str(e)
        result['time'] = time.time() - start_time
        logger.error(f"    ✗ Error: {str(e)}")
    
    return result

def test_protein(pdb_id: str, protein_info: dict):
    """Test a protein with specified mutations."""
    logger.info(f"\n{'='*80}")
    logger.info(f"Testing: {pdb_id.upper()} - {protein_info['name']}")
    logger.info(f"{'='*80}")
    
    pdb_path = Path('data/pdb_cache') / f'{pdb_id}.pdb'
    
    if not pdb_path.exists():
        logger.error(f"PDB file not found: {pdb_path}")
        return {'successful': 0, 'total': 0, 'results': []}
    
    # Process mutations in parallel
    results = []
    mutations_list = protein_info['mutations']
    
    with ProcessPoolExecutor(max_workers=min(MAX_WORKERS, len(mutations_list))) as executor:
        # Submit all mutations for this protein
        future_to_mutation = {
            executor.submit(test_mutation, str(pdb_path), mut_str, desc): (mut_str, desc)
            for mut_str, desc in mutations_list
        }
        
        # Collect results as they complete
        for future in as_completed(future_to_mutation):
            mut_str, desc = future_to_mutation[future]
            try:
                result = future.result()
                results.append(result)
                
                if result['success']:
                    logger.info(f"  ✓ {mut_str}: {result['ddg_calibrated']:.2f} kcal/mol ({result['time']:.1f}s)")
                else:
                    logger.error(f"  ✗ {mut_str}: {result['error']}")
            except Exception as e:
                logger.error(f"  ✗ {mut_str}: Exception - {str(e)}")
                results.append({
                    'mutation': mut_str,
                    'description': desc,
                    'success': False,
                    'error': str(e)
                })
    
    # Summary
    successful = sum(1 for r in results if r['success'])
    total = len(results)
    success_rate = (successful / total * 100) if total > 0 else 0
    
    logger.info(f"\n  {pdb_id.upper()} Summary: {successful}/{total} ({success_rate:.0f}%) successful")
    
    return {
        'pdb_id': pdb_id,
        'name': protein_info['name'],
        'successful': successful,
        'total': total,
        'results': results
    }

def main():
    """Run tests on 10 clean proteins."""
    logger.info("="*80)
    logger.info("MUTALYZE v5 - 10 CLEAN PROTEIN TEST SUITE")
    logger.info("="*80)
    logger.info(f"Testing {len(TEST_PROTEINS)} proteins with clean PDB structures")
    logger.info("Total mutations: " + str(sum(len(p['mutations']) for p in TEST_PROTEINS.values())))
    
    overall_start = time.time()
    all_results = []
    
    for pdb_id, protein_info in TEST_PROTEINS.items():
        protein_result = test_protein(pdb_id, protein_info)
        all_results.append(protein_result)
    
    # Overall summary
    total_elapsed = time.time() - overall_start
    total_successful = sum(r['successful'] for r in all_results)
    total_tests = sum(r['total'] for r in all_results)
    overall_rate = (total_successful / total_tests * 100) if total_tests > 0 else 0
    
    logger.info(f"\n{'='*80}")
    logger.info("OVERALL RESULTS")
    logger.info(f"{'='*80}")
    logger.info(f"Total Success: {total_successful}/{total_tests} ({overall_rate:.1f}%)")
    logger.info(f"Total Time: {total_elapsed:.1f} seconds")
    logger.info(f"Avg per mutation: {total_elapsed/total_tests:.1f} seconds")
    
    # Per-protein summary
    logger.info(f"\n{'='*80}")
    logger.info("PER-PROTEIN SUMMARY")
    logger.info(f"{'='*80}")
    logger.info(f"{'Protein':<15} {'Name':<25} {'Success':<15} {'Rate':<10}")
    logger.info("-"*80)
    
    for result in all_results:
        rate = (result['successful'] / result['total'] * 100) if result['total'] > 0 else 0
        logger.info(
            f"{result['pdb_id'].upper():<15} "
            f"{result['name']:<25} "
            f"{result['successful']}/{result['total']:<13} "
            f"{rate:>6.0f}%"
        )
    
    logger.info("="*80)
    
    if overall_rate >= 90:
        logger.info("✅ TEST PASSED - Success rate >= 90%")
        return 0
    else:
        logger.warning(f"⚠️  TEST WARNING - Success rate {overall_rate:.1f}% < 90%")
        return 1

if __name__ == '__main__':
    sys.exit(main())
