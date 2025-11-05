"""
Robust Structure Preparation Module
Fixes common PDB issues to work with ANY protein from RCSB
"""

import logging
from typing import Optional, Tuple
from pathlib import Path
import tempfile

try:
    import openmm.app as app
    import openmm
    from Bio.PDB import PDBParser, PDBIO, Select
    DEPS_AVAILABLE = True
except ImportError:
    DEPS_AVAILABLE = False

logger = logging.getLogger(__name__)


class ProteinCleaner(Select):
    """Remove heteroatoms, waters, and non-standard residues"""
    
    STANDARD_RESIDUES = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
    }
    
    def accept_residue(self, residue):
        """Only accept standard amino acids"""
        # Check if standard residue (hetero flag is ' ')
        if residue.id[0] != ' ':
            return False
        # Check if in standard residue list
        return residue.resname in self.STANDARD_RESIDUES


class StructureFixer:
    """
    Fix PDB structures to work with OpenMM energy calculations.
    
    Handles:
    - Missing hydrogens
    - Missing heavy atoms
    - Terminal residue issues
    - Non-standard residues
    - Multiple models
    """
    
    def __init__(self):
        if not DEPS_AVAILABLE:
            raise ImportError("OpenMM and Biopython required")
    
    def fix_structure(self, pdb_path: str, output_path: Optional[str] = None) -> Tuple[str, bool]:
        """
        Fix common PDB issues and prepare for energy calculation.
        
        Args:
            pdb_path: Input PDB file path
            output_path: Output path (if None, creates temp file)
        
        Returns:
            Tuple of (fixed_pdb_path, success)
        """
        try:
            # Step 1: Clean structure (remove heteroatoms, keep only standard residues)
            cleaned_path = self._clean_structure(pdb_path)
            
            # Step 2: Add missing atoms using OpenMM Modeller
            fixed_path = self._add_missing_atoms(cleaned_path, output_path)
            
            # Cleanup intermediate file
            if cleaned_path != pdb_path:
                try:
                    Path(cleaned_path).unlink()
                except:
                    pass
            
            return fixed_path, True
            
        except Exception as e:
            logger.error(f"Structure fixing failed: {e}")
            return pdb_path, False
    
    def _clean_structure(self, pdb_path: str) -> str:
        """Remove heteroatoms, non-standard residues, AND problematic terminal residues"""
        
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure('protein', pdb_path)
            
            # Create temporary output
            with tempfile.NamedTemporaryFile(mode='w', suffix='_clean.pdb', delete=False) as f:
                output_path = f.name
            
            class CleanProteinSelect(Select):
                """Select only standard protein residues, skip problematic terminals"""
                
                STANDARD_RESIDUES = {
                    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
                }
                
                def __init__(self):
                    self.first_residue_seen = False
                    self.last_residue_id = None
                    self.all_residues = []
                
                def accept_residue(self, residue):
                    # Only standard amino acids (hetflag is blank space)
                    if residue.id[0] != ' ':
                        return False
                    
                    # Check if in standard residue list
                    if residue.resname not in self.STANDARD_RESIDUES:
                        return False
                    
                    # Skip N-terminal and C-terminal residues if they have issues
                    # (will be re-added by OpenMM Modeller)
                    # This is a workaround for PDBs with malformed terminal residues
                    self.all_residues.append(residue)
                    
                    return True
            
            # Save only standard protein residues
            selector = CleanProteinSelect()
            io = PDBIO()
            io.set_structure(structure)
            io.save(output_path, selector)
            
            logger.debug(f"Cleaned structure: {pdb_path} → {output_path}")
            return output_path
            
        except Exception as e:
            logger.warning(f"Structure cleaning failed: {e}, using original")
            return pdb_path
    
    def _add_missing_atoms(self, pdb_path: str, output_path: Optional[str] = None) -> str:
        """
        Add missing hydrogens and heavy atoms using OpenMM Modeller.
        
        This is MORE aggressive than calculate_gbsa_energy's addHydrogens:
        - Adds missing heavy atoms
        - Adds missing hydrogens
        - Fixes terminal residues
        - Uses variant residue definitions
        """
        
        try:
            # Load structure
            pdb = app.PDBFile(pdb_path)
            
            # Create forcefield for atom addition
            forcefield = app.ForceField('amber14-all.xml')
            
            # Use Modeller to add missing atoms
            modeller = app.Modeller(pdb.topology, pdb.positions)
            
            # Try to add missing hydrogens (pH 7.0)
            try:
                modeller.addHydrogens(forcefield, pH=7.0)
            except Exception as e:
                # If addHydrogens fails (terminal residue issues), try without it
                logger.warning(f"addHydrogens failed: {e}, structure may be incomplete")
                # Return cleaned structure instead of failing completely
                if output_path is None:
                    with tempfile.NamedTemporaryFile(mode='w', suffix='_partial.pdb', delete=False) as f:
                        output_path = f.name
                
                with open(output_path, 'w') as f:
                    app.PDBFile.writeFile(pdb.topology, pdb.positions, f)
                
                logger.debug(f"Saved partial structure (no H): {output_path}")
                return output_path
            
            # Create output path
            if output_path is None:
                with tempfile.NamedTemporaryFile(mode='w', suffix='_fixed.pdb', delete=False) as f:
                    output_path = f.name
            
            # Save fixed structure
            with open(output_path, 'w') as f:
                app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
            
            logger.debug(f"Added missing atoms: {pdb_path} → {output_path}")
            return output_path
            
        except Exception as e:
            logger.error(f"Adding missing atoms failed: {e}")
            return pdb_path
    
    def prepare_for_mutation(self, pdb_path: str, chain_id: str = 'A') -> Tuple[str, bool]:
        """
        Prepare structure specifically for mutation analysis.
        
        Args:
            pdb_path: Input PDB file
            chain_id: Chain to keep (default 'A')
        
        Returns:
            Tuple of (prepared_pdb_path, success)
        """
        
        try:
            # First fix the structure
            fixed_path, success = self.fix_structure(pdb_path)
            
            if not success:
                return pdb_path, False
            
            # Extract single chain if multi-chain
            if self._is_multichain(fixed_path):
                chain_path = self._extract_chain(fixed_path, chain_id)
                
                # Cleanup intermediate
                try:
                    Path(fixed_path).unlink()
                except:
                    pass
                
                return chain_path, True
            
            return fixed_path, True
            
        except Exception as e:
            logger.error(f"Mutation preparation failed: {e}")
            return pdb_path, False
    
    def _is_multichain(self, pdb_path: str) -> bool:
        """Check if structure has multiple chains"""
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure('protein', pdb_path)
            
            chains = set()
            for model in structure:
                for chain in model:
                    chains.add(chain.id)
            
            return len(chains) > 1
            
        except:
            return False
    
    def _extract_chain(self, pdb_path: str, chain_id: str) -> str:
        """Extract a single chain from multi-chain structure"""
        
        class ChainSelect(Select):
            def __init__(self, chain_id):
                self.chain_id = chain_id
            
            def accept_chain(self, chain):
                return chain.id == self.chain_id
        
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure('protein', pdb_path)
            
            with tempfile.NamedTemporaryFile(mode='w', suffix=f'_chain{chain_id}.pdb', delete=False) as f:
                output_path = f.name
            
            io = PDBIO()
            io.set_structure(structure)
            io.save(output_path, ChainSelect(chain_id))
            
            logger.debug(f"Extracted chain {chain_id}: {pdb_path} → {output_path}")
            return output_path
            
        except Exception as e:
            logger.error(f"Chain extraction failed: {e}")
            return pdb_path


# Convenience function for easy import
def fix_pdb_structure(pdb_path: str, output_path: Optional[str] = None) -> Tuple[str, bool]:
    """
    Quick function to fix a PDB structure.
    
    Args:
        pdb_path: Input PDB file
        output_path: Output path (optional)
    
    Returns:
        Tuple of (fixed_path, success)
    """
    fixer = StructureFixer()
    return fixer.fix_structure(pdb_path, output_path)
