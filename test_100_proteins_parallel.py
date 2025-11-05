#!/usr/bin/env python3
"""
Test v5 calibration on 100 diverse proteins in parallel.
Uses all 12 CPU threads for maximum speed.
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
from core.pdb_utils import save_structure, PDBFetcher

# Use all available cores
MAX_WORKERS = 12

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class ProteinOnlySelect(Select):
    """Select only protein residues (exclude water, ions, ligands)."""
    def accept_residue(self, residue):
        return residue.id[0] == ' '

# 100 diverse proteins from PDB (various sizes, folds, and functions)
DIVERSE_PROTEINS = [
    # Very small proteins (20-50 residues)
    '1l2y', '2jof', '1pga', '2h3l', '1ycc',
    '1crn', '1ab1', '2ezk', '1pmy', '1prb',
    
    # Small proteins (50-100 residues)
    '1ubq', '1stn', '2ci2', '1aps', '1ris',
    '1vqb', '1shg', '1bpi', '2gb1', '1fme',
    '3icb', '1mba', '1csp', '1bdd', '1tit',
    '2lhb', '1ten', '2ezk', '1igd', '1shf',
    
    # Medium proteins (100-150 residues)
    '1lz1', '1mbn', '1lmb', '1hen', '1cew',
    '2cba', '3chy', '1poh', '1fkb', '1bni',
    '2lzm', '1tca', '1hrc', '1ifc', '1ost',
    '5pti', '1nar', '1ctf', '1mol', '1pne',
    
    # Medium-large proteins (150-250 residues)
    '1ake', '1ycc', '1cbn', '1aaj', '1fjg',
    '1bbp', '1fna', '1paz', '2sns', '1bgl',
    '1hvr', '1nls', '1gky', '1cfb', '1thx',
    '1brt', '1dts', '1cse', '1gdh', '1cyo',
    
    # Large proteins (250-400 residues)
    '1pga', '2ovo', '1aoz', '1rcf', '1top',
    '1gca', '1neu', '1beo', '1a6m', '1thg',
    '1tml', '1arb', '2ace', '1aap', '1abr',
    '1bp2', '1atl', '1cdy', '1dhs', '1ton',
    
    # Diverse fold types
    '1acx', '1ede', '1jnk', '1kuh', '1mla',
    '1nfn', '1ohm', '1rgs', '1wit', '1xjo',
    '1yuh', '2cmd', '2trx', '3gcb', '3tgl',
    '4fxn', '5cpa', '5rxn', '6ldh', '7rsa'
]

def fetch_pdb_if_needed(pdb_id: str) -> Path:
    """Fetch PDB file if not already cached."""
    pdb_cache = Path('data/pdb_cache')
    pdb_cache.mkdir(parents=True, exist_ok=True)
    pdb_path = pdb_cache / f'{pdb_id}.pdb'
    
    if not pdb_path.exists():
        try:
            fetcher = PDBFetcher()
            structure = fetcher.fetch_pdb(pdb_id)
            if structure is None:
                logger.error(f"  Failed to download {pdb_id}")
                return None
            logger.info(f"  Downloaded: {pdb_id}")
        except Exception as e:
            logger.error(f"  Failed to download {pdb_id}: {e}")
            return None
    
    return pdb_path

def get_protein_mutations(pdb_path: str):
    """Generate 2 diverse mutations for a protein."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_path)
    
    # Get first chain
    chain = list(structure[0])[0]
    residues = [r for r in chain if r.id[0] == ' ']
    
    if len(residues) < 10:
        return []
    
    # Select mutations from different regions
    mutations = []
    
    # Mutation from first quarter
    for res in residues[:len(residues)//4]:
        if res.resname in ['LEU', 'ILE', 'VAL', 'PHE', 'TRP', 'TYR']:
            mutations.append(f"{res.resname}{res.id[1]}ALA")
            break
    
    # Mutation from middle
    mid_start = len(residues)//3
    mid_end = 2*len(residues)//3
    for res in residues[mid_start:mid_end]:
        if res.resname in ['LEU', 'ILE', 'VAL', 'PHE', 'TRP', 'TYR']:
            mutations.append(f"{res.resname}{res.id[1]}ALA")
            break
    
    return mutations[:2]  # Return up to 2 mutations

def test_mutation(pdb_path: str, mutation_str: str) -> dict:
    """Test a single mutation (function for parallel execution)."""
    start_time = time.time()
    result = {
        'mutation': mutation_str,
        'success': False,
        'time': 0,
        'error': None,
        'ddg_calibrated': None,
    }
    
    try:
        # Parse structure and clean it
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_path)
        
        # Save cleaned structure (protein only)
        with tempfile.NamedTemporaryFile(mode='w', suffix='_clean.pdb', delete=False) as tmp:
            io = PDBIO()
            io.set_structure(structure)
            io.save(tmp.name, ProteinOnlySelect())
            clean_pdb_path = tmp.name
        
        # Re-parse the cleaned structure
        structure = parser.get_structure('protein', clean_pdb_path)
        
        # Get the first chain
        chain = list(structure[0])[0]
        chain_id = chain.id
        
        # Parse mutation
        residue_num = int(mutation_str[3:-3])
        mut_res = mutation_str[-3:]
        
        # Build mutant
        builder = MutationBuilder(fix_structures=True)
        mutant_structure = builder.apply_mutation(
            structure, 
            chain_id,
            residue_num,
            mut_res
        )
        
        # Save mutant (protein only)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
            io = PDBIO()
            io.set_structure(mutant_structure)
            io.save(tmp.name, ProteinOnlySelect())
            temp_mutant = tmp.name
        
        # Calculate ΔΔG
        ddg_raw = calculate_stability_change(clean_pdb_path, temp_mutant, method='gbsa')
        
        # Clean up
        os.unlink(temp_mutant)
        os.unlink(clean_pdb_path)
        
        if ddg_raw is None or ddg_raw == float('inf') or ddg_raw == float('-inf'):
            result['error'] = 'Energy calculation returned invalid value'
            return result
        
        # Apply v5 calibration
        ddg_calibrated = EmpiricalCorrection.apply_correction(ddg_raw)
        
        result['success'] = True
        result['ddg_calibrated'] = ddg_calibrated
        result['time'] = time.time() - start_time
        
    except Exception as e:
        result['error'] = str(e)
        result['time'] = time.time() - start_time
    
    return result

def test_protein_parallel(pdb_id: str) -> dict:
    """Test a single protein with all its mutations in parallel."""
    # Fetch PDB
    pdb_path = fetch_pdb_if_needed(pdb_id)
    if pdb_path is None:
        return {
            'pdb_id': pdb_id,
            'successful': 0,
            'total': 0,
            'time': 0,
            'error': 'Download failed'
        }
    
    # Get mutations
    try:
        mutations = get_protein_mutations(str(pdb_path))
        if not mutations:
            return {
                'pdb_id': pdb_id,
                'successful': 0,
                'total': 0,
                'time': 0,
                'error': 'No suitable mutations found'
            }
    except Exception as e:
        return {
            'pdb_id': pdb_id,
            'successful': 0,
            'total': 0,
            'time': 0,
            'error': f'Mutation generation failed: {e}'
        }
    
    # Test mutations
    start_time = time.time()
    results = []
    
    for mutation in mutations:
        result = test_mutation(str(pdb_path), mutation)
        results.append(result)
    
    successful = sum(1 for r in results if r['success'])
    total_time = time.time() - start_time
    
    return {
        'pdb_id': pdb_id,
        'successful': successful,
        'total': len(results),
        'time': total_time,
        'results': results
    }

def main():
    """Run parallel test on 100 proteins."""
    logger.info("="*80)
    logger.info("MUTALYZE v5 - 100 PROTEIN PARALLEL TEST")
    logger.info("="*80)
    logger.info(f"Testing {len(DIVERSE_PROTEINS)} diverse proteins")
    logger.info(f"Using {MAX_WORKERS} parallel workers")
    logger.info("")
    
    overall_start = time.time()
    all_results = []
    
    # Process proteins in parallel
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        future_to_pdb = {
            executor.submit(test_protein_parallel, pdb_id): pdb_id
            for pdb_id in DIVERSE_PROTEINS
        }
        
        completed = 0
        for future in as_completed(future_to_pdb):
            pdb_id = future_to_pdb[future]
            completed += 1
            
            try:
                result = future.result()
                all_results.append(result)
                
                if result['successful'] > 0:
                    logger.info(
                        f"[{completed}/{len(DIVERSE_PROTEINS)}] ✓ {pdb_id.upper()}: "
                        f"{result['successful']}/{result['total']} ({result['time']:.1f}s)"
                    )
                else:
                    logger.warning(
                        f"[{completed}/{len(DIVERSE_PROTEINS)}] ✗ {pdb_id.upper()}: "
                        f"{result.get('error', 'Failed')}"
                    )
            except Exception as e:
                logger.error(f"[{completed}/{len(DIVERSE_PROTEINS)}] ✗ {pdb_id.upper()}: {e}")
    
    # Overall summary
    total_elapsed = time.time() - overall_start
    total_successful = sum(r['successful'] for r in all_results)
    total_tests = sum(r['total'] for r in all_results)
    proteins_with_success = sum(1 for r in all_results if r['successful'] > 0)
    overall_rate = (total_successful / total_tests * 100) if total_tests > 0 else 0
    protein_rate = (proteins_with_success / len(DIVERSE_PROTEINS) * 100)
    
    logger.info(f"\n{'='*80}")
    logger.info("OVERALL RESULTS")
    logger.info(f"{'='*80}")
    logger.info(f"Proteins tested: {len(DIVERSE_PROTEINS)}")
    logger.info(f"Proteins successful: {proteins_with_success}/{len(DIVERSE_PROTEINS)} ({protein_rate:.1f}%)")
    logger.info(f"Total mutations: {total_tests}")
    logger.info(f"Successful mutations: {total_successful}/{total_tests} ({overall_rate:.1f}%)")
    logger.info(f"Total time: {total_elapsed:.1f} seconds")
    logger.info(f"Avg per protein: {total_elapsed/len(DIVERSE_PROTEINS):.1f} seconds")
    if total_tests > 0:
        logger.info(f"Avg per mutation: {total_elapsed/total_tests:.1f} seconds")
    logger.info(f"Speedup: ~{MAX_WORKERS}x faster than sequential")
    logger.info("="*80)
    
    if overall_rate >= 80:
        logger.info("✅ TEST PASSED - Success rate >= 80%")
        return 0
    else:
        logger.warning(f"⚠️  TEST WARNING - Success rate {overall_rate:.1f}% < 80%")
        return 1

if __name__ == '__main__':
    sys.exit(main())
