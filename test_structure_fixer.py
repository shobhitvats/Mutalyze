#!/usr/bin/env python3
"""
Test Structure Fixer on Problematic PDBs
Verifies the fixer handles missing atoms, hydrogens, and terminal caps
"""

import logging
from pathlib import Path
from core.structure_fixer import StructureFixer, fix_pdb_structure

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def test_on_existing_pdbs():
    """Test fixer on existing PDB files in dataset"""
    
    # Find PDB files
    pdb_dir = Path('data/pdb_cache')
    if not pdb_dir.exists():
        logger.error(f"PDB directory not found: {pdb_dir}")
        return
    
    pdb_files = list(pdb_dir.glob('*.pdb'))
    logger.info(f"Found {len(pdb_files)} PDB files to test")
    
    # Test fixer
    fixer = StructureFixer()
    
    results = []
    for pdb_file in pdb_files[:5]:  # Test first 5
        logger.info(f"\nTesting: {pdb_file.name}")
        
        try:
            fixed_path, success = fixer.fix_structure(str(pdb_file))
            
            if success:
                logger.info(f"✓ Fixed successfully: {fixed_path}")
                
                # Compare file sizes
                original_size = pdb_file.stat().st_size
                fixed_size = Path(fixed_path).stat().st_size
                
                logger.info(f"  Original: {original_size} bytes")
                logger.info(f"  Fixed: {fixed_size} bytes")
                logger.info(f"  Change: {fixed_size - original_size:+d} bytes")
                
                results.append({
                    'file': pdb_file.name,
                    'success': True,
                    'size_change': fixed_size - original_size
                })
            else:
                logger.warning(f"✗ Fixing failed: {pdb_file.name}")
                results.append({
                    'file': pdb_file.name,
                    'success': False,
                    'size_change': 0
                })
                
        except Exception as e:
            logger.error(f"✗ Error: {e}")
            results.append({
                'file': pdb_file.name,
                'success': False,
                'error': str(e)
            })
    
    # Summary
    logger.info("\n" + "="*60)
    logger.info("SUMMARY")
    logger.info("="*60)
    
    success_count = sum(1 for r in results if r['success'])
    logger.info(f"Success rate: {success_count}/{len(results)} ({100*success_count/len(results):.1f}%)")
    
    avg_size_change = sum(r.get('size_change', 0) for r in results) / len(results)
    logger.info(f"Average size increase: {avg_size_change:+.0f} bytes")
    
    for result in results:
        status = "✓" if result['success'] else "✗"
        logger.info(f"{status} {result['file']}")


def test_convenience_function():
    """Test the convenience function"""
    logger.info("\nTesting convenience function fix_pdb_structure()")
    
    pdb_dir = Path('data/pdb_cache')
    test_file = next(pdb_dir.glob('*.pdb'), None)
    
    if test_file:
        logger.info(f"Testing on: {test_file.name}")
        fixed_path, success = fix_pdb_structure(str(test_file))
        
        if success:
            logger.info(f"✓ Convenience function works: {fixed_path}")
        else:
            logger.warning("✗ Convenience function failed")
    else:
        logger.warning("No PDB file found for testing")


if __name__ == '__main__':
    logger.info("Structure Fixer Test Suite")
    logger.info("="*60)
    
    test_on_existing_pdbs()
    test_convenience_function()
    
    logger.info("\n" + "="*60)
    logger.info("Testing complete!")
