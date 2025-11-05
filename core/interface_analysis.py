"""
Interface Analysis Module
Detects ligand proximity, chain interfaces, and active sites
"""

import numpy as np
from typing import List, Dict, Tuple, Optional, Set
from Bio.PDB import Structure, NeighborSearch
import logging

logger = logging.getLogger(__name__)


class InterfaceAnalyzer:
    """Analyzes protein interfaces, ligand binding sites, and active sites"""
    
    def __init__(self, contact_threshold: float = 6.0):
        """
        Initialize interface analyzer
        
        Args:
            contact_threshold: Distance threshold for contacts in Angstroms
        """
        self.contact_threshold = contact_threshold
    
    def detect_interface_residues(self, structure: Structure.Structure, 
                                  chain1: str, chain2: str) -> List[Dict]:
        """
        Detect residues at the interface between two chains
        
        Args:
            structure: Biopython Structure
            chain1, chain2: Chain IDs to analyze
            
        Returns:
            List of interface residue dictionaries
        """
        interface_residues = []
        
        for model in structure:
            if chain1 not in model or chain2 not in model:
                continue
            
            ch1 = model[chain1]
            ch2 = model[chain2]
            
            # Get all atoms from each chain
            atoms1 = [atom for residue in ch1 for atom in residue if residue.id[0] == ' ']
            atoms2 = [atom for residue in ch2 for atom in residue if residue.id[0] == ' ']
            
            # Build neighbor search for chain 2
            ns = NeighborSearch(atoms2)
            
            # Find chain 1 residues near chain 2
            interface_res_ids = set()
            for atom in atoms1:
                neighbors = ns.search(atom.coord, self.contact_threshold)
                if neighbors:
                    residue = atom.get_parent()
                    interface_res_ids.add((chain1, residue.id[1]))
            
            # Build neighbor search for chain 1
            ns = NeighborSearch(atoms1)
            
            # Find chain 2 residues near chain 1
            for atom in atoms2:
                neighbors = ns.search(atom.coord, self.contact_threshold)
                if neighbors:
                    residue = atom.get_parent()
                    interface_res_ids.add((chain2, residue.id[1]))
            
            # Convert to list of dicts
            for chain_id, resid in interface_res_ids:
                chain = model[chain_id]
                for residue in chain:
                    if residue.id[1] == resid and residue.id[0] == ' ':
                        interface_residues.append({
                            'chain': chain_id,
                            'resid': resid,
                            'resname': residue.resname,
                            'type': 'interface'
                        })
        
        logger.info(f"Found {len(interface_residues)} interface residues between {chain1} and {chain2}")
        return interface_residues
    
    def detect_all_interfaces(self, structure: Structure.Structure) -> Dict[Tuple[str, str], List[Dict]]:
        """
        Detect all chain-chain interfaces in structure
        
        Returns:
            Dictionary mapping chain pairs to interface residues
        """
        interfaces = {}
        
        # Get all chain IDs
        chains = []
        for model in structure:
            chains = [chain.id for chain in model]
            break
        
        # Analyze all chain pairs
        for i, chain1 in enumerate(chains):
            for chain2 in chains[i+1:]:
                interface_res = self.detect_interface_residues(structure, chain1, chain2)
                if interface_res:
                    interfaces[(chain1, chain2)] = interface_res
        
        logger.info(f"Found {len(interfaces)} chain-chain interfaces")
        return interfaces
    
    def detect_ligand_proximity(self, structure: Structure.Structure, 
                               proximity_threshold: float = 5.0) -> Dict[str, List[Dict]]:
        """
        Detect residues near ligands (HETATM)
        
        Args:
            structure: Biopython Structure
            proximity_threshold: Distance threshold in Angstroms
            
        Returns:
            Dictionary mapping ligand names to nearby residues
        """
        ligand_proximity = {}
        
        for model in structure:
            # Separate protein residues and ligands
            protein_atoms = []
            ligand_groups = {}
            
            for chain in model:
                for residue in chain:
                    # Protein residues
                    if residue.id[0] == ' ':
                        protein_atoms.extend([(atom, chain.id, residue.id[1], residue.resname) 
                                            for atom in residue])
                    # Ligands (exclude water)
                    elif residue.id[0].startswith('H_') and residue.resname not in ['HOH', 'WAT']:
                        ligand_key = f"{residue.resname}_{chain.id}_{residue.id[1]}"
                        ligand_groups[ligand_key] = {
                            'atoms': list(residue.get_atoms()),
                            'resname': residue.resname,
                            'chain': chain.id,
                            'resid': residue.id[1]
                        }
            
            # For each ligand, find nearby protein residues
            for ligand_key, ligand_info in ligand_groups.items():
                nearby_residues = []
                
                # Build neighbor search for protein atoms
                ns = NeighborSearch([atom for atom, _, _, _ in protein_atoms])
                
                # Check each ligand atom
                nearby_res_ids = set()
                for lig_atom in ligand_info['atoms']:
                    neighbors = ns.search(lig_atom.coord, proximity_threshold)
                    for neighbor in neighbors:
                        # Find which residue this atom belongs to
                        for atom, ch_id, res_id, res_name in protein_atoms:
                            if atom == neighbor:
                                nearby_res_ids.add((ch_id, res_id, res_name))
                
                # Convert to list
                for chain_id, resid, resname in nearby_res_ids:
                    # Calculate minimum distance
                    min_dist = float('inf')
                    for lig_atom in ligand_info['atoms']:
                        for atom, ch_id, res_id, _ in protein_atoms:
                            if ch_id == chain_id and res_id == resid:
                                dist = np.linalg.norm(lig_atom.coord - atom.coord)
                                if dist < min_dist:
                                    min_dist = dist
                    
                    nearby_residues.append({
                        'chain': chain_id,
                        'resid': resid,
                        'resname': resname,
                        'distance': min_dist,
                        'ligand': ligand_info['resname']
                    })
                
                if nearby_residues:
                    ligand_proximity[ligand_key] = nearby_residues
        
        logger.info(f"Analyzed proximity for {len(ligand_proximity)} ligands")
        return ligand_proximity
    
    def classify_residue_location(self, structure: Structure.Structure, 
                                 chain_id: str, resid: int) -> str:
        """
        Classify residue as core, surface, or interface
        
        Args:
            structure: Biopython Structure
            chain_id: Chain ID
            resid: Residue number
            
        Returns:
            Classification: 'core', 'surface', 'interface', or 'unknown'
        """
        for model in structure:
            if chain_id not in model:
                return 'unknown'
            
            chain = model[chain_id]
            target_residue = None
            
            # Find target residue
            for residue in chain:
                if residue.id[1] == resid and residue.id[0] == ' ':
                    target_residue = residue
                    break
            
            if not target_residue:
                return 'unknown'
            
            # Count nearby atoms from same vs different chains
            target_atoms = list(target_residue.get_atoms())
            all_atoms = []
            
            for ch in model:
                for res in ch:
                    if res.id[0] == ' ':  # Only protein residues
                        for atom in res:
                            all_atoms.append((atom, ch.id))
            
            ns = NeighborSearch([a for a, _ in all_atoms])
            
            same_chain_contacts = 0
            diff_chain_contacts = 0
            total_contacts = 0
            
            for atom in target_atoms:
                neighbors = ns.search(atom.coord, 6.0)
                for neighbor in neighbors:
                    neighbor_chain = None
                    for a, ch in all_atoms:
                        if a == neighbor:
                            neighbor_chain = ch
                            break
                    
                    if neighbor_chain:
                        if neighbor_chain == chain_id:
                            same_chain_contacts += 1
                        else:
                            diff_chain_contacts += 1
                        total_contacts += 1
            
            # Classify based on contacts
            if diff_chain_contacts > 5:
                return 'interface'
            elif total_contacts > 15:
                return 'core'
            else:
                return 'surface'
        
        return 'unknown'
    
    def annotate_mutations_with_context(self, structure: Structure.Structure, 
                                       mutations: List[Dict]) -> List[Dict]:
        """
        Annotate mutations with structural context (interface, ligand proximity, etc.)
        
        Args:
            structure: Biopython Structure
            mutations: List of mutation dictionaries
            
        Returns:
            Annotated mutation list
        """
        # Detect ligand proximity
        ligand_proximity = self.detect_ligand_proximity(structure)
        
        # Create lookup for ligand-proximal residues
        ligand_proximal = {}
        for ligand_key, residues in ligand_proximity.items():
            for res in residues:
                key = (res['chain'], res['resid'])
                if key not in ligand_proximal:
                    ligand_proximal[key] = []
                ligand_proximal[key].append({
                    'ligand': res['ligand'],
                    'distance': res['distance']
                })
        
        # Annotate each mutation
        annotated = []
        for mutation in mutations:
            chain = mutation['chain']
            resid = mutation['resid']
            
            # Add location classification
            location = self.classify_residue_location(structure, chain, resid)
            mutation['location'] = location
            
            # Add ligand proximity info
            key = (chain, resid)
            if key in ligand_proximal:
                mutation['ligand_proximity'] = ligand_proximal[key]
                mutation['near_ligand'] = True
            else:
                mutation['ligand_proximity'] = []
                mutation['near_ligand'] = False
            
            annotated.append(mutation)
        
        logger.info(f"Annotated {len(annotated)} mutations with structural context")
        return annotated


def detect_binding_site_residues(structure: Structure.Structure, 
                                 proximity_threshold: float = 5.0) -> List[Dict]:
    """Convenience function to detect all ligand binding site residues"""
    analyzer = InterfaceAnalyzer()
    ligand_prox = analyzer.detect_ligand_proximity(structure, proximity_threshold)
    
    all_residues = []
    for ligand_key, residues in ligand_prox.items():
        all_residues.extend(residues)
    
    return all_residues
