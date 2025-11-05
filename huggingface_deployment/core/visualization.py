"""
Visualization Module
3D rendering and plotting for protein structures and mutation analysis
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict, Optional, Tuple
import logging

logger = logging.getLogger(__name__)

# Set plotting style
sns.set_style("whitegrid")


class StructureVisualizer:
    """3D visualization of protein structures"""
    
    @staticmethod
    def generate_py3dmol_view(pdb_file: str, 
                              width: int = 800, 
                              height: int = 600,
                              style: str = 'cartoon') -> str:
        """
        Generate Py3Dmol HTML view of structure
        
        Args:
            pdb_file: Path to PDB file
            width, height: Viewer dimensions
            style: Visualization style (cartoon, stick, sphere, etc.)
            
        Returns:
            HTML string for embedding
        """
        try:
            import py3Dmol
            
            with open(pdb_file, 'r') as f:
                pdb_data = f.read()
            
            view = py3Dmol.view(width=width, height=height)
            view.addModel(pdb_data, 'pdb')
            view.setStyle({style: {'color': 'spectrum'}})
            view.zoomTo()
            
            return view._make_html()
            
        except ImportError:
            logger.error("py3Dmol not available")
            return "<p>3D visualization not available</p>"
        except Exception as e:
            logger.error(f"Visualization failed: {e}")
            return f"<p>Error: {e}</p>"
    
    @staticmethod
    def color_by_conservation(pdb_file: str,
                            conservation_scores: Dict[int, float],
                            chain_id: str = 'A') -> str:
        """
        Color structure by conservation scores
        
        Args:
            pdb_file: PDB file path
            conservation_scores: Dict mapping resid to conservation score
            chain_id: Chain to color
            
        Returns:
            HTML visualization
        """
        try:
            import py3Dmol
            
            with open(pdb_file, 'r') as f:
                pdb_data = f.read()
            
            view = py3Dmol.view(width=800, height=600)
            view.addModel(pdb_data, 'pdb')
            
            # Color by conservation
            for resid, score in conservation_scores.items():
                # Convert score to color (blue = low, red = high conservation)
                color = plt.cm.RdYlBu_r(score)
                hex_color = '#{:02x}{:02x}{:02x}'.format(
                    int(color[0]*255), int(color[1]*255), int(color[2]*255)
                )
                
                view.setStyle({'chain': chain_id, 'resi': resid},
                            {'cartoon': {'color': hex_color}})
            
            view.zoomTo()
            return view._make_html()
            
        except Exception as e:
            logger.error(f"Conservation coloring failed: {e}")
            return "<p>Error in visualization</p>"
    
    @staticmethod
    def color_by_ddg(pdb_file: str,
                    ddg_data: List[Dict],
                    chain_id: str = 'A') -> str:
        """
        Color structure by ΔΔG values
        
        Args:
            pdb_file: PDB file path
            ddg_data: List of mutations with ΔΔG values
            chain_id: Chain to visualize
            
        Returns:
            HTML visualization
        """
        try:
            import py3Dmol
            
            with open(pdb_file, 'r') as f:
                pdb_data = f.read()
            
            view = py3Dmol.view(width=800, height=600)
            view.addModel(pdb_data, 'pdb')
            
            # Get max ΔΔG per residue
            residue_ddg = {}
            for mutation in ddg_data:
                if mutation.get('chain') == chain_id:
                    resid = mutation.get('resid')
                    ddg = abs(mutation.get('ddg', 0))
                    
                    if resid not in residue_ddg:
                        residue_ddg[resid] = ddg
                    else:
                        residue_ddg[resid] = max(residue_ddg[resid], ddg)
            
            # Normalize and color
            if residue_ddg:
                max_ddg = max(residue_ddg.values())
                
                for resid, ddg in residue_ddg.items():
                    # Normalize to 0-1
                    normalized = min(ddg / max(max_ddg, 1.0), 1.0)
                    
                    # Blue (stable) to red (unstable)
                    color = plt.cm.RdYlBu_r(normalized)
                    hex_color = '#{:02x}{:02x}{:02x}'.format(
                        int(color[0]*255), int(color[1]*255), int(color[2]*255)
                    )
                    
                    view.setStyle({'chain': chain_id, 'resi': resid},
                                {'cartoon': {'color': hex_color}})
            
            view.zoomTo()
            return view._make_html()
            
        except Exception as e:
            logger.error(f"ΔΔG coloring failed: {e}")
            return "<p>Error in visualization</p>"


class PlotGenerator:
    """Generate analysis plots"""
    
    @staticmethod
    def plot_ddg_distribution(ddg_values: List[float], 
                             save_path: Optional[str] = None) -> plt.Figure:
        """
        Plot ΔΔG distribution histogram
        
        Args:
            ddg_values: List of ΔΔG values
            save_path: Optional path to save figure
            
        Returns:
            Matplotlib figure
        """
        # Filter out inf and nan values
        import numpy as np
        ddg_values = [x for x in ddg_values if np.isfinite(x)]
        
        if not ddg_values:
            # Create empty plot with message
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.text(0.5, 0.5, 'No valid ΔΔG values to plot',
                   ha='center', va='center', fontsize=14)
            ax.set_xlabel('ΔΔG (kcal/mol)', fontsize=12)
            ax.set_ylabel('Count', fontsize=12)
            ax.set_title('Distribution of ΔΔG Values', fontsize=14, fontweight='bold')
            return fig
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        ax.hist(ddg_values, bins=50, edgecolor='black', alpha=0.7)
        ax.axvline(x=0, color='red', linestyle='--', label='Neutral')
        ax.axvline(x=-1, color='blue', linestyle='--', alpha=0.5, label='Stabilizing threshold')
        ax.axvline(x=1, color='orange', linestyle='--', alpha=0.5, label='Destabilizing threshold')
        
        ax.set_xlabel('ΔΔG (kcal/mol)', fontsize=12)
        ax.set_ylabel('Count', fontsize=12)
        ax.set_title(f'Distribution of ΔΔG Values (n={len(ddg_values)})', 
                    fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    @staticmethod
    def plot_conservation_vs_ddg(mutations: List[Dict],
                                conservation_key: str = 'conservation',
                                ddg_key: str = 'ddg',
                                save_path: Optional[str] = None) -> plt.Figure:
        """
        Scatter plot of conservation vs ΔΔG
        
        Args:
            mutations: List of mutations with conservation and ΔΔG
            conservation_key: Key for conservation scores
            ddg_key: Key for ΔΔG values
            save_path: Optional save path
            
        Returns:
            Matplotlib figure
        """
        conservation = [m[conservation_key] for m in mutations 
                       if conservation_key in m and ddg_key in m]
        ddg = [abs(m[ddg_key]) for m in mutations 
              if conservation_key in m and ddg_key in m]
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        scatter = ax.scatter(conservation, ddg, alpha=0.6, s=50)
        
        # Add trend line
        if len(conservation) > 1:
            z = np.polyfit(conservation, ddg, 1)
            p = np.poly1d(z)
            x_line = np.linspace(min(conservation), max(conservation), 100)
            ax.plot(x_line, p(x_line), "r--", alpha=0.8, label='Trend line')
        
        ax.set_xlabel('Conservation Score', fontsize=12)
        ax.set_ylabel('|ΔΔG| (kcal/mol)', fontsize=12)
        ax.set_title('Conservation vs. Mutation Sensitivity', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    @staticmethod
    def plot_mutation_heatmap(mutations: List[Dict],
                             ddg_key: str = 'ddg',
                             save_path: Optional[str] = None) -> plt.Figure:
        """
        Heatmap of mutations by position and amino acid
        
        Args:
            mutations: List of mutations
            ddg_key: Key for ΔΔG values
            save_path: Optional save path
            
        Returns:
            Matplotlib figure
        """
        # Organize data into matrix
        amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                      'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        
        # Get unique positions
        positions = sorted(set(m.get('resid', 0) for m in mutations))
        
        if not positions:
            return plt.figure()
        
        # Sample positions if too many
        if len(positions) > 50:
            step = len(positions) // 50
            positions = positions[::step]
        
        # Create matrix
        matrix = np.zeros((len(amino_acids), len(positions)))
        
        for mutation in mutations:
            if ddg_key not in mutation:
                continue
            
            mut_aa = mutation.get('mut_aa', '')
            resid = mutation.get('resid', 0)
            ddg = mutation[ddg_key]
            
            if mut_aa in amino_acids and resid in positions:
                aa_idx = amino_acids.index(mut_aa)
                pos_idx = positions.index(resid)
                matrix[aa_idx, pos_idx] = ddg
        
        # Plot heatmap
        fig, ax = plt.subplots(figsize=(14, 8))
        
        im = ax.imshow(matrix, cmap='RdYlBu_r', aspect='auto', 
                      vmin=-3, vmax=3)
        
        ax.set_xticks(range(len(positions)))
        ax.set_xticklabels(positions, rotation=90)
        ax.set_yticks(range(len(amino_acids)))
        ax.set_yticklabels(amino_acids)
        
        ax.set_xlabel('Residue Position', fontsize=12)
        ax.set_ylabel('Mutant Amino Acid', fontsize=12)
        ax.set_title('Mutation ΔΔG Heatmap', fontsize=14, fontweight='bold')
        
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('ΔΔG (kcal/mol)', fontsize=12)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    @staticmethod
    def plot_category_pie(mutations: List[Dict],
                         ddg_key: str = 'ddg',
                         save_path: Optional[str] = None) -> plt.Figure:
        """
        Pie chart of mutation categories
        
        Args:
            mutations: List of mutations
            ddg_key: Key for ΔΔG values
            save_path: Optional save path
            
        Returns:
            Matplotlib figure
        """
        # Categorize mutations
        categories = {
            'Stabilizing': 0,
            'Neutral': 0,
            'Destabilizing': 0,
            'Highly Destabilizing': 0
        }
        
        for mutation in mutations:
            if ddg_key not in mutation:
                continue
            
            ddg = mutation[ddg_key]
            
            if ddg < -1.0:
                categories['Stabilizing'] += 1
            elif -1.0 <= ddg <= 1.0:
                categories['Neutral'] += 1
            elif 1.0 < ddg <= 2.0:
                categories['Destabilizing'] += 1
            else:
                categories['Highly Destabilizing'] += 1
        
        # Plot
        fig, ax = plt.subplots(figsize=(8, 8))
        
        colors = ['#3498db', '#95a5a6', '#e67e22', '#e74c3c']
        wedges, texts, autotexts = ax.pie(
            categories.values(),
            labels=categories.keys(),
            autopct='%1.1f%%',
            colors=colors,
            startangle=90
        )
        
        ax.set_title('Mutation Impact Distribution', fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    @staticmethod
    def plot_ddg_by_position(df):
        """
        Plot ΔΔG values by residue position
        
        Args:
            df: DataFrame with 'position' and 'ddg' columns
            
        Returns:
            Matplotlib figure
        """
        import pandas as pd
        
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Group by position and aggregate
        if 'position' in df.columns and 'ddg' in df.columns:
            position_ddg = df.groupby('position')['ddg'].agg(['mean', 'std', 'count']).reset_index()
            
            # Plot
            ax.bar(position_ddg['position'], position_ddg['mean'], 
                   yerr=position_ddg['std'], alpha=0.7, capsize=3)
            ax.axhline(y=0, color='red', linestyle='--', label='Neutral')
            ax.axhline(y=-1, color='blue', linestyle='--', alpha=0.5, label='Stabilizing')
            ax.axhline(y=1, color='orange', linestyle='--', alpha=0.5, label='Destabilizing')
            
            ax.set_xlabel('Residue Position', fontsize=12)
            ax.set_ylabel('Mean ΔΔG (kcal/mol)', fontsize=12)
            ax.set_title('ΔΔG by Residue Position', fontsize=14, fontweight='bold')
            ax.legend()
            ax.grid(True, alpha=0.3, axis='y')
        else:
            ax.text(0.5, 0.5, 'No position data available',
                   ha='center', va='center', fontsize=14)
        
        plt.tight_layout()
        return fig


def create_all_plots(mutations: List[Dict], output_dir: str = "plots"):
    """Generate all standard plots for mutation analysis"""
    from pathlib import Path
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    plotter = PlotGenerator()
    
    # Extract ΔΔG values
    ddg_values = [m['ddg'] for m in mutations if 'ddg' in m]
    
    if ddg_values:
        # Distribution
        plotter.plot_ddg_distribution(ddg_values, 
                                     str(output_path / 'ddg_distribution.png'))
        
        # Pie chart
        plotter.plot_category_pie(mutations,
                                 str(output_path / 'category_pie.png'))
        
        # Heatmap
        plotter.plot_mutation_heatmap(mutations,
                                     str(output_path / 'mutation_heatmap.png'))
        
        logger.info(f"Generated plots in {output_dir}")
