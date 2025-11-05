"""
Conservation Analysis Module
Implements BLAST search, MSA building, and conservation scoring
"""

import numpy as np
from typing import List, Dict, Optional, Tuple
from pathlib import Path
import subprocess
import tempfile
import logging
from collections import Counter
import math

logger = logging.getLogger(__name__)


class ConservationAnalyzer:
    """Analyzes sequence conservation using MSA"""
    
    def __init__(self, blast_db: Optional[str] = None):
        """
        Initialize conservation analyzer
        
        Args:
            blast_db: Path to BLAST database (e.g., 'nr', 'swissprot')
        """
        self.blast_db = blast_db or 'nr'
        self.blast_available = self._check_blast()
        self.mafft_available = self._check_mafft()
    
    def _check_blast(self) -> bool:
        """Check if BLAST+ is available"""
        try:
            subprocess.run(['blastp', '-version'], 
                         capture_output=True, check=True)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning("BLAST+ not found. Sequence search disabled.")
            return False
    
    def _check_mafft(self) -> bool:
        """Check if MAFFT is available"""
        try:
            subprocess.run(['mafft', '--version'], 
                         capture_output=True, check=True)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning("MAFFT not found. MSA building disabled.")
            return False
    
    def extract_sequence_from_pdb(self, structure, chain_id: str) -> str:
        """
        Extract amino acid sequence from PDB structure
        
        Args:
            structure: Biopython Structure
            chain_id: Chain identifier
            
        Returns:
            Single-letter amino acid sequence
        """
        from Bio.PDB.Polypeptide import PPBuilder
        
        ppb = PPBuilder()
        sequence = ""
        
        for model in structure:
            if chain_id in model:
                chain = model[chain_id]
                for pp in ppb.build_peptides(chain):
                    sequence += str(pp.get_sequence())
                break
        
        return sequence
    
    def blast_search(self, sequence: str, 
                    num_hits: int = 100,
                    evalue_threshold: float = 1e-5) -> Optional[List[str]]:
        """
        Perform BLAST search to find homologous sequences
        
        Args:
            sequence: Query sequence
            num_hits: Number of hits to retrieve
            evalue_threshold: E-value threshold
            
        Returns:
            List of homologous sequences or None if failed
        """
        if not self.blast_available:
            logger.error("BLAST not available")
            return None
        
        try:
            # Create temporary files
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as query_file:
                query_file.write(f">query\n{sequence}\n")
                query_path = query_file.name
            
            output_path = tempfile.mktemp(suffix='.txt')
            
            # Run BLAST
            cmd = [
                'blastp',
                '-query', query_path,
                '-db', self.blast_db,
                '-out', output_path,
                '-outfmt', '6 sseq',  # Get sequences only
                '-max_target_seqs', str(num_hits),
                '-evalue', str(evalue_threshold)
            ]
            
            logger.info(f"Running BLAST search (this may take a while)...")
            result = subprocess.run(cmd, capture_output=True, timeout=300)
            
            if result.returncode != 0:
                logger.error(f"BLAST failed: {result.stderr.decode()}")
                return None
            
            # Parse results
            sequences = [sequence]  # Include query
            with open(output_path, 'r') as f:
                for line in f:
                    seq = line.strip()
                    if seq and seq not in sequences:
                        sequences.append(seq)
            
            # Cleanup
            Path(query_path).unlink(missing_ok=True)
            Path(output_path).unlink(missing_ok=True)
            
            logger.info(f"Found {len(sequences)-1} homologous sequences")
            return sequences
            
        except Exception as e:
            logger.error(f"BLAST search failed: {e}")
            return None
    
    def build_msa(self, sequences: List[str], output_file: Optional[str] = None) -> Optional[List[str]]:
        """
        Build multiple sequence alignment using MAFFT
        
        Args:
            sequences: List of sequences to align
            output_file: Optional output file path
            
        Returns:
            List of aligned sequences or None if failed
        """
        if not self.mafft_available:
            logger.warning("MAFFT not available, using simple alignment")
            return sequences  # Return unaligned
        
        try:
            # Create input FASTA file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as input_file:
                for i, seq in enumerate(sequences):
                    input_file.write(f">seq{i}\n{seq}\n")
                input_path = input_file.name
            
            # Run MAFFT
            output_path = output_file or tempfile.mktemp(suffix='.aln')
            
            cmd = ['mafft', '--auto', '--quiet', input_path]
            
            logger.info("Building multiple sequence alignment with MAFFT...")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            
            if result.returncode != 0:
                logger.error(f"MAFFT failed: {result.stderr}")
                return None
            
            # Save output
            with open(output_path, 'w') as f:
                f.write(result.stdout)
            
            # Parse aligned sequences
            aligned_seqs = []
            current_seq = ""
            for line in result.stdout.split('\n'):
                if line.startswith('>'):
                    if current_seq:
                        aligned_seqs.append(current_seq)
                        current_seq = ""
                else:
                    current_seq += line.strip()
            if current_seq:
                aligned_seqs.append(current_seq)
            
            # Cleanup
            Path(input_path).unlink(missing_ok=True)
            
            logger.info(f"MSA built with {len(aligned_seqs)} sequences")
            return aligned_seqs
            
        except Exception as e:
            logger.error(f"MSA building failed: {e}")
            return None
    
    def calculate_conservation_scores(self, msa: List[str]) -> List[float]:
        """
        Calculate per-position conservation scores using Shannon entropy
        
        Args:
            msa: Multiple sequence alignment (list of aligned sequences)
            
        Returns:
            List of conservation scores (0-1, higher = more conserved)
        """
        if not msa:
            return []
        
        alignment_length = len(msa[0])
        num_sequences = len(msa)
        
        conservation_scores = []
        
        for pos in range(alignment_length):
            # Get residues at this position
            residues = [seq[pos] for seq in msa if pos < len(seq)]
            
            # Remove gaps
            residues = [r for r in residues if r not in ['-', 'X']]
            
            if not residues:
                conservation_scores.append(0.0)
                continue
            
            # Calculate Shannon entropy
            counter = Counter(residues)
            total = len(residues)
            
            entropy = 0.0
            for count in counter.values():
                if count > 0:
                    freq = count / total
                    entropy -= freq * math.log2(freq)
            
            # Normalize entropy to 0-1 scale
            # Maximum entropy for 20 amino acids is log2(20) ≈ 4.32
            max_entropy = math.log2(20)
            normalized_entropy = entropy / max_entropy if max_entropy > 0 else 0
            
            # Convert to conservation score (1 - entropy)
            conservation = 1.0 - normalized_entropy
            conservation_scores.append(conservation)
        
        return conservation_scores
    
    def map_conservation_to_structure(self, structure, chain_id: str,
                                     conservation_scores: List[float]) -> Dict[int, float]:
        """
        Map conservation scores to structure residue numbers
        
        Args:
            structure: Biopython Structure
            chain_id: Chain ID
            conservation_scores: Per-position conservation scores
            
        Returns:
            Dictionary mapping residue numbers to conservation scores
        """
        from Bio.PDB.Polypeptide import PPBuilder
        
        residue_conservation = {}
        
        ppb = PPBuilder()
        for model in structure:
            if chain_id not in model:
                continue
            
            chain = model[chain_id]
            peptides = ppb.build_peptides(chain)
            
            seq_index = 0
            for peptide in peptides:
                for residue in peptide:
                    resid = residue.id[1]
                    if seq_index < len(conservation_scores):
                        residue_conservation[resid] = conservation_scores[seq_index]
                        seq_index += 1
            break
        
        return residue_conservation
    
    def analyze_structure_conservation(self, structure, chain_id: str,
                                      num_homologs: int = 50) -> Optional[Dict[int, float]]:
        """
        Complete conservation analysis pipeline for a structure
        
        Args:
            structure: Biopython Structure
            chain_id: Chain to analyze
            num_homologs: Number of homologous sequences to retrieve
            
        Returns:
            Dictionary mapping residue numbers to conservation scores
        """
        # Extract sequence
        sequence = self.extract_sequence_from_pdb(structure, chain_id)
        if not sequence:
            logger.error(f"Could not extract sequence for chain {chain_id}")
            return None
        
        logger.info(f"Extracted sequence of length {len(sequence)}")
        
        # BLAST search
        homologs = self.blast_search(sequence, num_hits=num_homologs)
        if not homologs:
            logger.warning("BLAST search failed, using single sequence")
            homologs = [sequence]
        
        # Build MSA
        msa = self.build_msa(homologs)
        if not msa:
            logger.error("MSA building failed")
            return None
        
        # Calculate conservation
        conservation_scores = self.calculate_conservation_scores(msa)
        
        # Map to structure
        residue_conservation = self.map_conservation_to_structure(
            structure, chain_id, conservation_scores
        )
        
        logger.info(f"Calculated conservation for {len(residue_conservation)} residues")
        return residue_conservation


def calculate_simple_conservation(sequences: List[str]) -> List[float]:
    """
    Simple conservation calculation without external tools
    
    Args:
        sequences: List of sequences (assumed aligned)
        
    Returns:
        Conservation scores per position
    """
    analyzer = ConservationAnalyzer()
    return analyzer.calculate_conservation_scores(sequences)
