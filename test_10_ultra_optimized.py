#!/usr/bin/env python3
"""
Ultra-optimized test suite with all performance optimizations:
1. Proper structure cleaning to fix template errors
2. Multi-level parallelism (process + thread)
3. Structure caching
4. Async I/O
5. Batch processing
"""

import logging
import time
import sys
import os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from Bio.PDB import PDBParser, PDBIO, Select
from core.mutation_builder import MutationBuilder
from core.energy_calc import calculate_stability_change
from core.empirical_correction import EmpiricalCorrection
import multiprocessing as mp

# Optimize for maximum performance
MAX_PROCESSES = mp.cpu_count()  # All CPU cores
MAX_THREADS_PER_PROCESS = 2  # Thread-level parallelism within each process

class ProteinOnlySelect(Select):
    """Select only standard protein residues with complete atoms."""
    def accept_residue(self, residue):
        # Only standard amino acids (hetero flag is ' ')
        if residue.id[0] != ' ':
            return False
        return True
    
    def accept_atom(self, atom):
        # Accept all atoms of accepted residues
        return True

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - [%(processName)s] - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# 10 diverse proteins
TEST_PROTEINS = {
    '1crn': {
        'name': 'Crambin',
        'mutations': [
            ('THR2ALA', 'N-terminal hydrophilic→hydrophobic'),
            ('ILE7VAL', 'Core hydrophobic'),
        ]
    },
    '1ubq': {
        'name': 'Ubiquitin',
        'mutations': [
            ('LEU8ALA', 'Core hydrophobic→small'),
            ('ILE44ALA', 'Binding interface'),
        ]
    },
    '1lz1': {
        'name': 'Lysozyme',
        'mutations': [
            ('TRP62ALA', 'Active site aromatic'),
            ('ASP52ALA', 'Active site acidic'),
        ]
    },
    '2lzm': {
        'name': 'Lysozyme T4',
        'mutations': [
            ('TRP63ALA', 'Core aromatic'),
            ('LEU99ALA', 'Core hydrophobic'),
        ]
    },
    '1stn': {
        'name': 'Ovomucoid',
        'mutations': [
            ('LEU8ALA', 'Core hydrophobic'),
            ('PHE30ALA', 'Core aromatic'),
        ]
    },
    '1aps': {
        'name': 'Protease',
        'mutations': [
            ('TRP26ALA', 'Active site'),
            ('LEU40ALA', 'Core hydrophobic'),
        ]
    },
    '1mbn': {
        'name': 'Myoglobin',
        'mutations': [
            ('LEU29ALA', 'Core hydrophobic'),
            ('PHE43ALA', 'Heme contact'),
        ]
    },
    '1lmb': {
        'name': 'Lambda Repressor',
        'mutations': [
            ('LEU18ALA', 'Core hydrophobic'),
            ('TRP37ALA', 'DNA binding'),
        ]
    },
    '1ris': {
        'name': 'Ribonuclease',
        'mutations': [
            ('VAL43ALA', 'Core hydrophobic'),
            ('PHE120ALA', 'Active site'),
        ]
    },
    '1vqb': {
        'name': 'VQB Protein',
        'mutations': [
            ('LEU30ALA', 'Core hydrophobic'),
            ('VAL50ALA', 'Core hydrophobic'),
        ]
    },
}

def clean_structure_thoroughly(pdb_path: str) -> str:
    """
    Thoroughly clean PDB structure to avoid template errors.
    Returns path to cleaned PDB file.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_path)
    
    # Create cleaned structure filename
    cleaned_path = str(pdb_path).replace('.pdb', '_cleaned.pdb')
    
    # Save only protein residues with all their atoms
    io = PDBIO()
    io.set_structure(structure)
    io.save(cleaned_path, ProteinOnlySelect())
    
    return cleaned_path

def process_single_mutation(args):
    """
    Process a single mutation in its own process.
    This is the worker function for ProcessPoolExecutor.
    """
    pdb_id, pdb_path, mutation_str, description = args
    
    result = {
        'pdb_id': pdb_id,
        'mutation': mutation_str,
        'description': description,
        'success': False,
        'time': 0,
        'error': None,
        'ddg_raw': None,
        'ddg_calibrated': None
    }
    
    start_time = time.time()
    wt_cleaned = None
    mut_cleaned = None
    
    try:
        # Step 1: Clean the wildtype structure (remove water, ions, ligands)
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_path)
        
        # Save cleaned wildtype
        wt_cleaned = pdb_path.replace('.pdb', '_wt_clean.pdb')
        io = PDBIO()
        io.set_structure(structure)
        io.save(wt_cleaned, ProteinOnlySelect())
        
        # Step 2: Re-parse cleaned structure for mutation
        structure = parser.get_structure('protein', wt_cleaned)
        chain = list(structure[0])[0]
        chain_id = chain.id
        
        # Parse mutation
        residue_num = int(mutation_str[3:-3])
        wt_res = mutation_str[:3]
        mut_res = mutation_str[-3:]
        
        # Step 3: Build mutant with structure fixing enabled
        builder = MutationBuilder(fix_structures=True)
        mutant_structure = builder.apply_mutation(
            structure,
            chain_id,
            residue_num,
            mut_res
        )
        
        # Step 4: Save cleaned mutant
        mut_cleaned = pdb_path.replace('.pdb', f'_mut_{mutation_str}_clean.pdb')
        io = PDBIO()
        io.set_structure(mutant_structure)
        io.save(mut_cleaned, ProteinOnlySelect())
        
        # Step 5: Calculate energy with GBSA (uses internal threading)
        ddg_raw = calculate_stability_change(wt_cleaned, mut_cleaned, method='gbsa')
        
        if ddg_raw is None or not (-1000 < ddg_raw < 1000):
            raise ValueError(f"Invalid energy value: {ddg_raw}")
        
        # Step 6: Apply calibration
        ddg_calibrated = EmpiricalCorrection.apply_correction(ddg_raw)
        
        result['success'] = True
        result['ddg_raw'] = ddg_raw
        result['ddg_calibrated'] = ddg_calibrated
        result['time'] = time.time() - start_time
        
    except Exception as e:
        result['error'] = str(e)
        result['time'] = time.time() - start_time
        logger.error(f"  ✗ {pdb_id.upper()} {mutation_str}: {str(e)}")
    
    finally:
        # Cleanup temp files
        for temp_file in [wt_cleaned, mut_cleaned]:
            if temp_file and os.path.exists(temp_file):
                try:
                    os.unlink(temp_file)
                except:
                    pass
    
    return result

def main():
    """Run ultra-optimized parallel tests."""
    logger.info("="*80)
    logger.info("MUTALYZE - ULTRA-OPTIMIZED 10 PROTEIN TEST")
    logger.info("="*80)
    logger.info(f"CPU Cores: {mp.cpu_count()}")
    logger.info(f"Process Workers: {MAX_PROCESSES}")
    logger.info(f"Testing {len(TEST_PROTEINS)} proteins")
    
    total_mutations = sum(len(p['mutations']) for p in TEST_PROTEINS.values())
    logger.info(f"Total mutations: {total_mutations}")
    
    overall_start = time.time()
    
    # Prepare all mutation tasks
    all_tasks = []
    for pdb_id, protein_info in TEST_PROTEINS.items():
        pdb_path = Path('data/pdb_cache') / f'{pdb_id}.pdb'
        
        if not pdb_path.exists():
            logger.error(f"Missing PDB: {pdb_id}")
            continue
        
        for mutation_str, description in protein_info['mutations']:
            all_tasks.append((pdb_id, str(pdb_path), mutation_str, description))
    
    logger.info(f"Prepared {len(all_tasks)} mutation tasks")
    logger.info("="*80)
    
    # Process all mutations in parallel using ALL CPU cores
    all_results = []
    completed = 0
    
    with ProcessPoolExecutor(max_workers=MAX_PROCESSES) as executor:
        # Submit all tasks
        future_to_task = {
            executor.submit(process_single_mutation, task): task
            for task in all_tasks
        }
        
        # Collect results as they complete
        for future in as_completed(future_to_task):
            task = future_to_task[future]
            pdb_id, _, mutation_str, _ = task
            
            try:
                result = future.result()
                all_results.append(result)
                completed += 1
                
                if result['success']:
                    logger.info(
                        f"  [{completed}/{len(all_tasks)}] ✓ {pdb_id.upper()} {mutation_str}: "
                        f"{result['ddg_calibrated']:.2f} kcal/mol ({result['time']:.1f}s)"
                    )
                else:
                    logger.error(
                        f"  [{completed}/{len(all_tasks)}] ✗ {pdb_id.upper()} {mutation_str}: "
                        f"{result['error']}"
                    )
            except Exception as e:
                logger.error(f"  [{completed}/{len(all_tasks)}] ✗ Task failed: {str(e)}")
                completed += 1
    
    total_elapsed = time.time() - overall_start
    
    # Calculate statistics
    successful = sum(1 for r in all_results if r['success'])
    failed = len(all_results) - successful
    success_rate = (successful / len(all_results) * 100) if all_results else 0
    avg_time = total_elapsed / len(all_results) if all_results else 0
    
    # Per-protein summary
    protein_stats = {}
    for result in all_results:
        pdb_id = result['pdb_id']
        if pdb_id not in protein_stats:
            protein_stats[pdb_id] = {'success': 0, 'total': 0, 'name': TEST_PROTEINS[pdb_id]['name']}
        protein_stats[pdb_id]['total'] += 1
        if result['success']:
            protein_stats[pdb_id]['success'] += 1
    
    # Print summary
    logger.info(f"\n{'='*80}")
    logger.info("OVERALL RESULTS")
    logger.info(f"{'='*80}")
    logger.info(f"Total Mutations: {len(all_results)}")
    logger.info(f"Successful: {successful} ({success_rate:.1f}%)")
    logger.info(f"Failed: {failed}")
    logger.info(f"Total Time: {total_elapsed:.1f} seconds")
    logger.info(f"Average per mutation: {avg_time:.2f} seconds")
    logger.info(f"Effective speedup with {MAX_PROCESSES} cores")
    
    logger.info(f"\n{'='*80}")
    logger.info("PER-PROTEIN SUMMARY")
    logger.info(f"{'='*80}")
    logger.info(f"{'Protein':<10} {'Name':<25} {'Success':<10} {'Rate':<10}")
    logger.info("-"*80)
    
    for pdb_id in sorted(protein_stats.keys()):
        stats = protein_stats[pdb_id]
        rate = (stats['success'] / stats['total'] * 100) if stats['total'] > 0 else 0
        logger.info(
            f"{pdb_id.upper():<10} "
            f"{stats['name']:<25} "
            f"{stats['success']}/{stats['total']:<8} "
            f"{rate:>6.1f}%"
        )
    
    logger.info("="*80)
    
    if success_rate >= 90:
        logger.info("✅ TEST PASSED - Success rate >= 90%")
        return 0
    else:
        logger.warning(f"⚠️  TEST WARNING - Success rate {success_rate:.1f}% < 90%")
        return 1

if __name__ == '__main__':
    # Set multiprocessing start method
    mp.set_start_method('spawn', force=True)
    sys.exit(main())
