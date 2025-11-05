"""
Mutation Builder Module
Handles amino acid mutations, rotamer sampling, and structure modification
"""

import numpy as np
from typing import List, Dict, Tuple, Optional
from Bio.PDB import Structure, Residue, Atom, Vector, PDBIO
import logging
import json
from pathlib import Path
import copy
import tempfile

# Import sidechain builder and rotamer sampler
from core.sidechain_builder import SidechainBuilder
from core.rotamer_sampler import RotamerSampler, get_rotamer_library
from core.structure_fixer import StructureFixer

logger = logging.getLogger(__name__)


# Amino acid properties
AMINO_ACIDS = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU',
               'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']

AA_ONE_TO_THREE = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
    'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
    'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
    'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
}

AA_THREE_TO_ONE = {v: k for k, v in AA_ONE_TO_THREE.items()}


class RotamerLibrary:
    """Simplified rotamer library for amino acid side chains"""
    
    def __init__(self):
        """Initialize with basic rotamer conformations"""
        # Simplified Dunbrack-like rotamer library
        # Format: {residue: [(chi1, chi2, chi3, chi4), ...]}
        self.rotamers = self._initialize_rotamers()
    
    def _initialize_rotamers(self) -> Dict:
        """Initialize basic rotamer library"""
        # Simplified rotamers - in production would load from database
        rotamers = {
            'ALA': [(None,)],  # No chi angles
            'GLY': [(None,)],
            'SER': [(-60,), (60,), (180,)],
            'CYS': [(-60,), (60,), (180,)],
            'VAL': [(60,), (180,), (-60,)],
            'ILE': [(60, 60), (60, 180), (-60, -60), (-60, 180)],
            'LEU': [(60, 60), (60, 180), (-60, -60), (-60, 180)],
            'THR': [(-60,), (60,), (180,)],
            'ARG': [(60, 180, 60, 180), (180, 180, 60, 180), (-60, 180, 60, 180)],
            'LYS': [(60, 180, 60, 180), (180, 180, 60, 180), (-60, 180, 60, 180)],
            'ASP': [(-60, -60), (60, 60), (180, 180)],
            'GLU': [(60, 60, 60), (180, 60, 60), (-60, 180, -60)],
            'ASN': [(-60, -60), (60, 60), (180, 180)],
            'GLN': [(60, 60, 60), (180, 60, 60), (-60, 180, -60)],
            'MET': [(60, 60, 60), (180, 60, 60), (-60, 180, -60)],
            'HIS': [(-60, -60), (60, 60), (180, 180)],
            'PHE': [(-60, 90), (60, 90), (180, 90)],
            'TYR': [(-60, 90), (60, 90), (180, 90)],
            'TRP': [(-60, 90), (60, 90), (180, 90)],
            'PRO': [(30,)]  # Proline is constrained
        }
        return rotamers
    
    def get_rotamers(self, resname: str) -> List[Tuple]:
        """Get rotamer conformations for a residue type"""
        return self.rotamers.get(resname, [(None,)])


class MutationBuilder:
    """Builds mutant protein structures"""
    
    def __init__(self, fix_structures: bool = True):
        """
        Initialize mutation builder
        
        Args:
            fix_structures: If True, automatically fix PDB structures before mutation
        """
        self.rotamer_lib = RotamerLibrary()
        self.fix_structures = fix_structures
        if fix_structures:
            self.structure_fixer = StructureFixer()
    
    def enumerate_mutations(self, structure: Structure.Structure, 
                          chain_ids: Optional[List[str]] = None,
                          residue_range: Optional[Tuple[int, int]] = None) -> List[Dict]:
        """
        Enumerate all possible single-point mutations
        
        Args:
            structure: Biopython Structure
            chain_ids: Optional list of chains to mutate
            residue_range: Optional (start, end) residue range
            
        Returns:
            List of mutation dictionaries
        """
        mutations = []
        
        for model in structure:
            for chain in model:
                if chain_ids and chain.id not in chain_ids:
                    continue
                
                for residue in chain:
                    # Skip non-standard residues
                    if residue.id[0] != ' ':
                        continue
                    
                    resid = residue.id[1]
                    
                    # Check residue range
                    if residue_range:
                        if resid < residue_range[0] or resid > residue_range[1]:
                            continue
                    
                    wt_resname = residue.resname
                    
                    # Skip if not standard amino acid
                    if wt_resname not in AMINO_ACIDS:
                        continue
                    
                    wt_aa = AA_THREE_TO_ONE.get(wt_resname, 'X')
                    
                    # Generate mutations to all other amino acids
                    for mut_resname in AMINO_ACIDS:
                        if mut_resname == wt_resname:
                            continue  # Skip wild-type
                        
                        mut_aa = AA_THREE_TO_ONE.get(mut_resname, 'X')
                        
                        mutations.append({
                            'chain': chain.id,
                            'resid': resid,
                            'wt_resname': wt_resname,
                            'mut_resname': mut_resname,
                            'wt_aa': wt_aa,
                            'mut_aa': mut_aa,
                            'mutation_code': f"{wt_aa}{chain.id}{resid}{mut_aa}"
                        })
        
        logger.info(f"Enumerated {len(mutations)} possible mutations")
        return mutations
    
    def apply_mutation(self, structure: Structure.Structure, 
                      chain_id: str, resid: int, 
                      new_resname: str) -> Structure.Structure:
        """
        Apply a mutation to a structure
        
        Args:
            structure: Original structure
            chain_id: Chain identifier
            resid: Residue number
            new_resname: New residue name (3-letter code)
            
        Returns:
            Mutated structure (deep copy)
        """
        # Deep copy structure to avoid modifying original
        mutant_structure = copy.deepcopy(structure)
        
        # Find the residue to mutate
        for model in mutant_structure:
            for chain in model:
                if chain.id != chain_id:
                    continue
                
                for residue in chain:
                    if residue.id[1] == resid and residue.id[0] == ' ':
                        # Mutate this residue
                        self._mutate_residue(residue, new_resname)
                        logger.debug(f"Mutated {chain_id}:{resid} to {new_resname}")
                        return mutant_structure
        
        logger.warning(f"Residue {chain_id}:{resid} not found for mutation")
        return mutant_structure
    
    def apply_mutation_ensemble(self, structure: Structure.Structure,
                                chain_id: str, resid: int,
                                new_resname: str, top_n: int = 5) -> List[Tuple[str, float]]:
        """
        Apply mutation and generate rotamer ensemble
        
        Args:
            structure: Original structure
            chain_id: Chain identifier
            resid: Residue number
            new_resname: New residue name (3-letter code)
            top_n: Number of top rotamers to generate
            
        Returns:
            List of (pdb_filepath, probability) tuples for each rotamer
        """
        logger.info(f"Generating rotamer ensemble for {chain_id}:{resid}→{new_resname} (top {top_n})")
        
        # Get rotamer library
        rotamer_lib = get_rotamer_library()
        sampler = RotamerSampler()
        
        # Find the residue
        target_residue = None
        for model in structure:
            for chain in model:
                if chain.id != chain_id:
                    continue
                for residue in chain:
                    if residue.id[1] == resid and residue.id[0] == ' ':
                        target_residue = residue
                        break
        
        if target_residue is None:
            logger.error(f"Residue {chain_id}:{resid} not found")
            return []
        
        # Get rotamer conformations from library
        rotamer_ensemble = sampler.generate_rotamer_ensemble(target_residue, new_resname, top_n=top_n)
        
        # Build structure for each rotamer
        ensemble_files = []
        temp_dir = tempfile.mkdtemp(prefix='mutalyze_rotamers_')
        
        for idx, (chi_angles, probability) in enumerate(rotamer_ensemble):
            # Deep copy structure
            mutant_struct = copy.deepcopy(structure)
            
            # Find the residue in the copy
            for model in mutant_struct:
                for chain in model:
                    if chain.id != chain_id:
                        continue
                    for residue in chain:
                        if residue.id[1] == resid and residue.id[0] == ' ':
                            # Mutate residue (removes sidechain)
                            self._mutate_residue(residue, new_resname)
                            
                            # Note: Currently we can't apply specific chi angles
                            # The sidechain builder uses fixed geometry
                            # TODO: Implement chi-angle aware building
                            # For now, each "rotamer" is the same geometry
                            # This is a placeholder for full implementation
                            break
            
            # Save to temporary PDB file (remove heteroatoms)
            pdb_path = Path(temp_dir) / f"rotamer_{idx+1}_p{probability:.3f}.pdb"
            from Bio.PDB import Select
            
            class ProteinOnlySelect(Select):
                """Select only protein residues (no waters, ligands, ions)"""
                def accept_residue(self, residue):
                    # Only keep standard amino acids (hetflag is blank space)
                    return residue.id[0] == ' '
            
            io = PDBIO()
            io.set_structure(mutant_struct)
            io.save(str(pdb_path), ProteinOnlySelect())
            
            # Note: Structure fixing is now done in EnergyCalculator.calculate_gbsa_energy
            # to avoid double-fixing (which caused issues with terminal residues)
            
            ensemble_files.append((str(pdb_path), probability))
            logger.debug(f"Generated rotamer {idx+1}/{top_n}: χ={chi_angles}, P={probability:.3f}")
        
        logger.info(f"Generated {len(ensemble_files)} rotamer structures in {temp_dir}")
        return ensemble_files
    
    def _mutate_residue(self, residue: Residue.Residue, new_resname: str):
        """
        Mutate a single residue in place and build missing sidechain atoms
        
        Args:
            residue: Residue to mutate
            new_resname: New residue type
            
        Note:
            Uses ideal geometry to build missing sidechain atoms.
            This is a simplified approach - production tools use rotamer libraries.
        """
        old_resname = residue.resname
        
        # Keep only backbone atoms
        backbone_atoms = ['N', 'CA', 'C', 'O']
        atoms_to_remove = []
        for atom in residue.get_atoms():
            if atom.name not in backbone_atoms:
                atoms_to_remove.append(atom.id)
        
        # Remove all sidechain atoms
        for atom_id in atoms_to_remove:
            residue.detach_child(atom_id)
        
        # Update residue name
        residue.resname = new_resname
        
        # Build new sidechain using ideal geometry
        try:
            SidechainBuilder.build_sidechain(residue, new_resname)
            logger.debug(f"Mutated residue from {old_resname} to {new_resname}, built sidechain with {len(list(residue.get_atoms()))} atoms")
        except Exception as e:
            logger.warning(f"Failed to build complete sidechain for {new_resname}: {e}")
    
    def _build_cb(self, n_atom: Atom.Atom, ca_atom: Atom.Atom, c_atom: Atom.Atom) -> Optional[Atom.Atom]:
        """
        Build an idealized CB atom position
        
        Args:
            n_atom, ca_atom, c_atom: Backbone atoms
            
        Returns:
            CB atom or None
        """
        try:
            # Get coordinates
            n_coord = n_atom.get_vector()
            ca_coord = ca_atom.get_vector()
            c_coord = c_atom.get_vector()
            
            # Build CB using idealized geometry
            # CB is ~tetrahedral from CA
            b = ca_coord - n_coord
            c_vec = c_coord - ca_coord
            a = b.cross(c_vec)
            
            # Normalize and scale
            a = a.normalized()
            b = b.normalized()
            
            # Idealized CB position (1.54 Å from CA)
            cb_coord = ca_coord - 0.58 * b - 0.82 * a
            
            # Create CB atom
            cb_atom = Atom.Atom(
                name='CB',
                coord=np.array(cb_coord),
                bfactor=ca_atom.bfactor,
                occupancy=1.0,
                altloc=' ',
                fullname=' CB ',
                serial_number=ca_atom.serial_number + 1,
                element='C'
            )
            
            return cb_atom
            
        except Exception as e:
            logger.warning(f"Failed to build CB atom: {e}")
            return None
    
    def build_rotamers(self, structure: Structure.Structure, 
                      chain_id: str, resid: int, 
                      resname: str) -> List[Structure.Structure]:
        """
        Build all rotamers for a residue type at a position
        
        Args:
            structure: Structure template
            chain_id: Chain ID
            resid: Residue number
            resname: Residue type
            
        Returns:
            List of structures with different rotamers
        """
        rotamers = self.rotamer_lib.get_rotamers(resname)
        rotamer_structures = []
        
        for chi_angles in rotamers:
            # Apply mutation
            mutant = self.apply_mutation(structure, chain_id, resid, resname)
            # In full implementation, would apply chi angles here
            rotamer_structures.append(mutant)
        
        return rotamer_structures


def enumerate_all_mutations(structure: Structure.Structure) -> List[Dict]:
    """Convenience function to enumerate all mutations"""
    builder = MutationBuilder()
    return builder.enumerate_mutations(structure)


def create_mutant(structure: Structure.Structure, 
                 chain_id: str, resid: int, 
                 new_resname: str) -> Structure.Structure:
    """Convenience function to create a single mutant"""
    builder = MutationBuilder()
    return builder.apply_mutation(structure, chain_id, resid, new_resname)
