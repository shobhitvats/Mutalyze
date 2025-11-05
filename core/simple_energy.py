"""
Lightweight Energy Calculator for Streamlit Cloud
Works without OpenMM - uses empirical statistical potentials
"""

import numpy as np
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

class SimplifiedEnergyCalculator:
    """
    Statistical potential-based ΔΔG calculator
    Works without OpenMM for Streamlit Cloud deployment
    """
    
    # Miyazawa-Jernigan contact potentials (statistical)
    # Units: kcal/mol
    CONTACT_POTENTIALS = {
        ('ALA', 'ALA'): -0.30, ('ALA', 'CYS'): -0.49, ('ALA', 'ASP'): -0.16, ('ALA', 'GLU'): -0.20,
        ('ALA', 'PHE'): -0.56, ('ALA', 'GLY'): -0.27, ('ALA', 'HIS'): -0.27, ('ALA', 'ILE'): -0.50,
        ('ALA', 'LYS'): -0.12, ('ALA', 'LEU'): -0.53, ('ALA', 'MET'): -0.49, ('ALA', 'ASN'): -0.07,
        ('ALA', 'PRO'): -0.25, ('ALA', 'GLN'): -0.15, ('ALA', 'ARG'): -0.17, ('ALA', 'SER'): -0.16,
        ('ALA', 'THR'): -0.22, ('ALA', 'VAL'): -0.46, ('ALA', 'TRP'): -0.40, ('ALA', 'TYR'): -0.31,
        
        ('CYS', 'CYS'): -0.96, ('CYS', 'ASP'): -0.26, ('CYS', 'GLU'): -0.22, ('CYS', 'PHE'): -0.67,
        ('CYS', 'GLY'): -0.30, ('CYS', 'HIS'): -0.41, ('CYS', 'ILE'): -0.59, ('CYS', 'LYS'): -0.05,
        ('CYS', 'LEU'): -0.61, ('CYS', 'MET'): -0.67, ('CYS', 'ASN'): -0.14, ('CYS', 'PRO'): -0.27,
        ('CYS', 'GLN'): -0.18, ('CYS', 'ARG'): -0.13, ('CYS', 'SER'): -0.18, ('CYS', 'THR'): -0.23,
        ('CYS', 'VAL'): -0.54, ('CYS', 'TRP'): -0.54, ('CYS', 'TYR'): -0.42,
        
        ('ASP', 'ASP'): 0.11, ('ASP', 'GLU'): 0.07, ('ASP', 'PHE'): -0.32, ('ASP', 'GLY'): -0.08,
        ('ASP', 'HIS'): -0.19, ('ASP', 'ILE'): -0.30, ('ASP', 'LYS'): 0.02, ('ASP', 'LEU'): -0.30,
        ('ASP', 'MET'): -0.30, ('ASP', 'ASN'): 0.06, ('ASP', 'PRO'): -0.07, ('ASP', 'GLN'): 0.04,
        ('ASP', 'ARG'): 0.05, ('ASP', 'SER'): 0.03, ('ASP', 'THR'): 0.00, ('ASP', 'VAL'): -0.27,
        ('ASP', 'TRP'): -0.23, ('ASP', 'TYR'): -0.17,
        
        # Add more as needed...
    }
    
    # Hydrophobicity scale (Kyte-Doolittle)
    HYDROPHOBICITY = {
        'ALA': 1.8, 'CYS': 2.5, 'ASP': -3.5, 'GLU': -3.5, 'PHE': 2.8,
        'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5, 'LYS': -3.9, 'LEU': 3.8,
        'MET': 1.9, 'ASN': -3.5, 'PRO': -1.6, 'GLN': -3.5, 'ARG': -4.5,
        'SER': -0.8, 'THR': -0.7, 'VAL': 4.2, 'TRP': -0.9, 'TYR': -1.3
    }
    
    # Volume change (Å³)
    VOLUME = {
        'ALA': 88.6, 'CYS': 108.5, 'ASP': 111.1, 'GLU': 138.4, 'PHE': 189.9,
        'GLY': 60.1, 'HIS': 153.2, 'ILE': 166.7, 'LYS': 168.6, 'LEU': 166.7,
        'MET': 162.9, 'ASN': 114.1, 'PRO': 112.7, 'GLN': 143.8, 'ARG': 173.4,
        'SER': 89.0, 'THR': 116.1, 'VAL': 140.0, 'TRP': 227.8, 'TYR': 193.6
    }
    
    def __init__(self):
        logger.info("Using simplified statistical potential calculator (no OpenMM required)")
    
    def calculate_ddg_simple(self, wt_resname: str, mut_resname: str, 
                            burial_factor: float = 0.5) -> float:
        """
        Calculate ΔΔG using statistical potentials
        
        Args:
            wt_resname: Wild-type residue (3-letter code)
            mut_resname: Mutant residue (3-letter code)
            burial_factor: 0.0 (surface) to 1.0 (buried)
            
        Returns:
            ΔΔG in kcal/mol
        """
        
        # Component 1: Hydrophobicity change
        hydro_wt = self.HYDROPHOBICITY.get(wt_resname, 0)
        hydro_mut = self.HYDROPHOBICITY.get(mut_resname, 0)
        delta_hydro = (hydro_mut - hydro_wt) * burial_factor * 0.3
        
        # Component 2: Volume change (steric clash penalty)
        vol_wt = self.VOLUME.get(wt_resname, 100)
        vol_mut = self.VOLUME.get(mut_resname, 100)
        delta_vol = abs(vol_mut - vol_wt) / 100.0 * burial_factor * 0.5
        
        # Component 3: Charge change penalty
        charged_wt = wt_resname in ['ASP', 'GLU', 'LYS', 'ARG', 'HIS']
        charged_mut = mut_resname in ['ASP', 'GLU', 'LYS', 'ARG', 'HIS']
        charge_penalty = 0
        if charged_wt != charged_mut:
            charge_penalty = 1.0 * burial_factor
        
        # Combine components
        ddg_raw = delta_hydro + delta_vol + charge_penalty
        
        return ddg_raw
    
    def estimate_burial(self, residue_id: int, chain_length: int) -> float:
        """
        Rough estimate of burial based on position
        (More accurate with actual structure analysis)
        """
        # Residues near terminals tend to be more exposed
        terminal_distance = min(residue_id, chain_length - residue_id)
        burial = min(1.0, terminal_distance / 10.0)
        return burial


def calculate_stability_change_simple(wt_resname: str, mut_resname: str,
                                     residue_id: int = 50, 
                                     chain_length: int = 100) -> float:
    """
    Simplified ΔΔG calculation without OpenMM
    
    This is a FALLBACK for when OpenMM is not available.
    Less accurate than full molecular dynamics but works on Streamlit Cloud.
    
    Args:
        wt_resname: Wild-type residue 3-letter code
        mut_resname: Mutant residue 3-letter code
        residue_id: Position in chain
        chain_length: Total chain length
    
    Returns:
        Estimated ΔΔG in kcal/mol
    """
    calc = SimplifiedEnergyCalculator()
    burial = calc.estimate_burial(residue_id, chain_length)
    ddg = calc.calculate_ddg_simple(wt_resname, mut_resname, burial)
    
    logger.info(f"Simple ΔΔG calculation: {wt_resname}→{mut_resname} = {ddg:.2f} kcal/mol (burial={burial:.2f})")
    
    return ddg
