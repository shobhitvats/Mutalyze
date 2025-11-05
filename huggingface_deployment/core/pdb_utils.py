"""
PDB Utilities Module
Handles fetching, parsing, and caching of PDB structures from RCSB
"""

import os
import requests
import gzip
import pickle
from pathlib import Path
from typing import Optional, Dict, List, Tuple
from Bio.PDB import PDBParser, PDBIO, Structure
import logging

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def is_aa(residue, standard=False):
    """
    Check if residue is an amino acid
    
    Args:
        residue: Biopython Residue object
        standard: If True, only check for standard amino acids
        
    Returns:
        True if residue is an amino acid
    """
    # Standard amino acids
    standard_aa = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU',
                   'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
    
    # Check if it's a standard residue (hetflag is blank space)
    if residue.id[0] != ' ':
        return False
    
    if standard:
        return residue.resname in standard_aa
    else:
        # Also accept modified amino acids
        return len(residue.resname) == 3


class PDBFetcher:
    """Handles downloading and caching PDB structures from RCSB"""
    
    def __init__(self, cache_dir: str = "data/pdb_cache"):
        """
        Initialize PDB fetcher with cache directory
        
        Args:
            cache_dir: Directory to cache downloaded PDB files
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.parser = PDBParser(QUIET=True)
        
    def fetch_pdb(self, pdb_id: str, force_download: bool = False) -> Optional[Structure.Structure]:
        """
        Fetch PDB structure from RCSB or local cache
        
        Args:
            pdb_id: PDB identifier (e.g., "6spo", "1crn")
            force_download: If True, bypass cache and download fresh
            
        Returns:
            Biopython Structure object or None if failed
        """
        pdb_id = pdb_id.lower().strip()
        cache_file = self.cache_dir / f"{pdb_id}.pdb"
        
        # Check cache first
        if cache_file.exists() and not force_download:
            logger.info(f"Loading {pdb_id} from cache")
            try:
                return self.parser.get_structure(pdb_id, str(cache_file))
            except Exception as e:
                logger.warning(f"Cache read failed: {e}. Attempting download.")
        
        # Download from RCSB
        logger.info(f"Downloading {pdb_id} from RCSB PDB")
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            # Save to cache
            with open(cache_file, 'w') as f:
                f.write(response.text)
            
            # Parse and return
            structure = self.parser.get_structure(pdb_id, str(cache_file))
            logger.info(f"Successfully fetched and cached {pdb_id}")
            return structure
            
        except requests.RequestException as e:
            logger.error(f"Failed to download {pdb_id}: {e}")
            return None
        except Exception as e:
            logger.error(f"Error processing {pdb_id}: {e}")
            return None
    
    def get_pdb_path(self, pdb_id: str) -> Optional[Path]:
        """Get path to cached PDB file"""
        pdb_id = pdb_id.lower().strip()
        cache_file = self.cache_dir / f"{pdb_id}.pdb"
        return cache_file if cache_file.exists() else None


class PDBAnalyzer:
    """Analyzes PDB structures and extracts metadata"""
    
    @staticmethod
    def get_chains(structure: Structure.Structure) -> List[str]:
        """Extract all chain IDs from structure"""
        chains = []
        for model in structure:
            for chain in model:
                chains.append(chain.id)
        return sorted(set(chains))
    
    @staticmethod
    def get_residues(structure: Structure.Structure, chain_id: Optional[str] = None) -> List[Dict]:
        """
        Extract residue information from structure
        
        Args:
            structure: Biopython Structure
            chain_id: Optional chain filter
            
        Returns:
            List of dictionaries with residue info
        """
        residues = []
        for model in structure:
            for chain in model:
                if chain_id and chain.id != chain_id:
                    continue
                    
                for residue in chain:
                    if is_aa(residue, standard=True):
                        residues.append({
                            'chain': chain.id,
                            'resname': residue.resname,
                            'resid': residue.id[1],
                            'insertion': residue.id[2].strip(),
                            'residue': residue
                        })
        return residues
    
    @staticmethod
    def get_ligands(structure: Structure.Structure) -> List[Dict]:
        """
        Extract heteroatoms (ligands, ions, etc.)
        
        Returns:
            List of ligand information dictionaries
        """
        ligands = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    # HETATM residues, excluding water
                    if residue.id[0].startswith('H_') or residue.id[0].startswith('W'):
                        if residue.resname not in ['HOH', 'WAT']:
                            ligands.append({
                                'chain': chain.id,
                                'resname': residue.resname,
                                'resid': residue.id[1],
                                'residue': residue
                            })
        return ligands
    
    @staticmethod
    def get_metadata(structure: Structure.Structure) -> Dict:
        """
        Extract comprehensive metadata from structure
        
        Returns:
            Dictionary with structure metadata
        """
        chains = PDBAnalyzer.get_chains(structure)
        all_residues = PDBAnalyzer.get_residues(structure)
        ligands = PDBAnalyzer.get_ligands(structure)
        
        metadata = {
            'pdb_id': structure.id,
            'num_models': len(structure),
            'num_chains': len(chains),
            'chains': chains,
            'num_residues': len(all_residues),
            'num_ligands': len(ligands),
            'ligands': [lig['resname'] for lig in ligands],
            'residue_range': {}
        }
        
        # Get residue range per chain
        for chain_id in chains:
            chain_residues = [r for r in all_residues if r['chain'] == chain_id]
            if chain_residues:
                resids = [r['resid'] for r in chain_residues]
                metadata['residue_range'][chain_id] = (min(resids), max(resids))
        
        return metadata
    
    @staticmethod
    def calculate_distance(residue1, residue2) -> float:
        """
        Calculate minimum distance between two residues
        
        Args:
            residue1, residue2: Biopython Residue objects
            
        Returns:
            Minimum atomic distance in Angstroms
        """
        min_dist = float('inf')
        for atom1 in residue1:
            for atom2 in residue2:
                dist = atom1 - atom2  # Biopython overloads '-' for distance
                if dist < min_dist:
                    min_dist = dist
        return min_dist


def save_structure(structure: Structure.Structure, filepath: str, remove_hetero: bool = True):
    """
    Save structure to PDB file
    
    Args:
        structure: Biopython Structure object
        filepath: Output PDB file path
        remove_hetero: If True, remove heteroatoms (waters, ligands, ions)
    """
    from Bio.PDB import Select
    
    class ProteinOnlySelect(Select):
        """Select only protein residues (no waters, ligands, ions)"""
        def accept_residue(self, residue):
            # Only keep standard amino acids (hetflag is blank space)
            return residue.id[0] == ' '
    
    io = PDBIO()
    io.set_structure(structure)
    
    if remove_hetero:
        io.save(str(filepath), ProteinOnlySelect())
    else:
        io.save(str(filepath))
    
    logger.info(f"Saved structure to {filepath}")


# Convenience functions
def fetch_pdb(pdb_id: str, cache_dir: str = "data/pdb_cache") -> Optional[Structure.Structure]:
    """Quick function to fetch a PDB structure"""
    fetcher = PDBFetcher(cache_dir)
    return fetcher.fetch_pdb(pdb_id)


def get_structure_info(pdb_id: str, cache_dir: str = "data/pdb_cache") -> Optional[Dict]:
    """Fetch structure and return metadata"""
    structure = fetch_pdb(pdb_id, cache_dir)
    if structure:
        return PDBAnalyzer.get_metadata(structure)
    return None
