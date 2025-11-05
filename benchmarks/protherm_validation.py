"""
ProThermDB Validation Module
Benchmark ΔΔG predictions against experimental data

This module provides tools to validate Mutalyze predictions against
the ProThermDB database of experimentally measured protein stability changes.
"""

import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import matthews_corrcoef, mean_squared_error, mean_absolute_error
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import logging
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class BenchmarkResults:
    """Results from benchmark validation"""
    pearson_r: float
    spearman_rho: float
    rmse: float
    mae: float
    mcc: float
    accuracy: float
    n_samples: int
    predictions: List[float]
    experimental: List[float]
    mutations: List[str]


class ProThermValidator:
    """Validates ΔΔG predictions against ProThermDB experimental data"""
    
    def __init__(self, data_file: Optional[str] = None):
        """
        Initialize validator
        
        Args:
            data_file: Path to curated ProThermDB CSV file
                      If None, will use default location
        """
        if data_file is None:
            data_file = "data/protherm_curated.csv"
        
        self.data_file = Path(data_file)
        self.data = None
        
        if self.data_file.exists():
            self.load_data()
        else:
            logger.warning(f"ProThermDB file not found: {data_file}")
            logger.info("You can download it from: https://web.iitm.ac.in/bioinfo2/prothermdb/")
    
    def load_data(self):
        """Load ProThermDB data from CSV"""
        try:
            self.data = pd.read_csv(self.data_file)
            logger.info(f"Loaded {len(self.data)} mutations from ProThermDB")
            
            # Expected columns: pdb_id, chain, position, wild_aa, mut_aa, ddg_exp, temperature, pH
            required_cols = ['pdb_id', 'chain', 'position', 'wild_aa', 'mut_aa', 'ddg_exp']
            missing_cols = [col for col in required_cols if col not in self.data.columns]
            
            if missing_cols:
                logger.warning(f"Missing columns in data: {missing_cols}")
            
        except Exception as e:
            logger.error(f"Failed to load ProThermDB data: {e}")
            self.data = None
    
    def create_test_dataset(self, n_samples: int = 100, 
                           pdb_id: Optional[str] = None,
                           random_state: int = 42) -> pd.DataFrame:
        """
        Create a test dataset from ProThermDB
        
        Args:
            n_samples: Number of mutations to sample
            pdb_id: If specified, only use mutations from this PDB
            random_state: Random seed for reproducibility
            
        Returns:
            DataFrame with sampled mutations
        """
        if self.data is None:
            logger.error("No data loaded")
            return pd.DataFrame()
        
        # Filter by PDB if specified
        if pdb_id:
            subset = self.data[self.data['pdb_id'] == pdb_id.lower()]
            logger.info(f"Filtered to {len(subset)} mutations in {pdb_id}")
        else:
            subset = self.data
        
        # Sample
        if len(subset) > n_samples:
            test_set = subset.sample(n=n_samples, random_state=random_state)
        else:
            test_set = subset
        
        logger.info(f"Created test set with {len(test_set)} mutations")
        return test_set
    
    def calculate_metrics(self, predictions: List[float], 
                         experimental: List[float]) -> BenchmarkResults:
        """
        Calculate all validation metrics
        
        Args:
            predictions: Predicted ΔΔG values
            experimental: Experimental ΔΔG values
            
        Returns:
            BenchmarkResults object with all metrics
        """
        pred = np.array(predictions)
        exp = np.array(experimental)
        
        # Correlation metrics
        pearson_r, _ = pearsonr(pred, exp)
        spearman_rho, _ = spearmanr(pred, exp)
        
        # Error metrics
        rmse = np.sqrt(mean_squared_error(exp, pred))
        mae = mean_absolute_error(exp, pred)
        
        # Classification metrics (destabilizing vs stabilizing)
        pred_sign = (pred > 0).astype(int)
        exp_sign = (exp > 0).astype(int)
        mcc = matthews_corrcoef(exp_sign, pred_sign)
        accuracy = np.mean(pred_sign == exp_sign)
        
        results = BenchmarkResults(
            pearson_r=pearson_r,
            spearman_rho=spearman_rho,
            rmse=rmse,
            mae=mae,
            mcc=mcc,
            accuracy=accuracy,
            n_samples=len(predictions),
            predictions=predictions,
            experimental=experimental,
            mutations=[]
        )
        
        return results
    
    def plot_correlation(self, results: BenchmarkResults, 
                        output_file: str = "correlation_plot.png",
                        title: str = "Mutalyze ΔΔG Predictions"):
        """
        Create correlation plot
        
        Args:
            results: BenchmarkResults object
            output_file: Path to save plot
            title: Plot title
        """
        fig, ax = plt.subplots(figsize=(8, 8))
        
        pred = np.array(results.predictions)
        exp = np.array(results.experimental)
        
        # Scatter plot
        ax.scatter(exp, pred, alpha=0.6, s=50, edgecolors='black', linewidth=0.5)
        
        # Diagonal line (perfect prediction)
        lims = [
            min(exp.min(), pred.min()) - 1,
            max(exp.max(), pred.max()) + 1
        ]
        ax.plot(lims, lims, 'r--', alpha=0.75, label='Perfect Prediction')
        
        # Labels and title
        ax.set_xlabel('Experimental ΔΔG (kcal/mol)', fontsize=12)
        ax.set_ylabel('Predicted ΔΔG (kcal/mol)', fontsize=12)
        ax.set_title(title, fontsize=14, fontweight='bold')
        
        # Add statistics text box
        stats_text = f"""N = {results.n_samples}
Pearson r = {results.pearson_r:.3f}
Spearman ρ = {results.spearman_rho:.3f}
RMSE = {results.rmse:.2f} kcal/mol
MAE = {results.mae:.2f} kcal/mol
MCC = {results.mcc:.3f}
Accuracy = {results.accuracy:.1%}"""
        
        ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
                fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved correlation plot to {output_file}")
        plt.close()
    
    def print_report(self, results: BenchmarkResults):
        """Print formatted benchmark report"""
        print("\n" + "="*70)
        print("MUTALYZE BENCHMARK RESULTS")
        print("="*70)
        print(f"\nDataset: {results.n_samples} mutations")
        print("\nCorrelation Metrics:")
        print(f"  Pearson r:      {results.pearson_r:6.3f}")
        print(f"  Spearman ρ:     {results.spearman_rho:6.3f}")
        print("\nError Metrics:")
        print(f"  RMSE:           {results.rmse:6.2f} kcal/mol")
        print(f"  MAE:            {results.mae:6.2f} kcal/mol")
        print("\nClassification Metrics (Stabilizing vs Destabilizing):")
        print(f"  MCC:            {results.mcc:6.3f}")
        print(f"  Accuracy:       {results.accuracy:6.1%}")
        print("\nLiterature Comparison (ProThermDB):")
        print("  FoldX:          r = 0.50-0.65, RMSE = 1.5-2.0 kcal/mol")
        print("  Rosetta:        r = 0.55-0.70, RMSE = 1.3-1.8 kcal/mol")
        print("  GBSA (basic):   r = 0.30-0.45, RMSE = 2.5-4.0 kcal/mol")
        print("="*70 + "\n")


def create_sample_protherm_data():
    """
    Create a sample ProThermDB-style dataset for testing
    
    In production, this would be replaced with actual ProThermDB data
    """
    # Sample mutations with realistic ΔΔG values
    data = {
        'pdb_id': ['1crn', '1crn', '1crn', '1crn', '1crn'],
        'chain': ['A', 'A', 'A', 'A', 'A'],
        'position': [25, 25, 25, 25, 25],
        'wild_aa': ['THR', 'THR', 'THR', 'THR', 'THR'],
        'mut_aa': ['ALA', 'ASP', 'LYS', 'GLY', 'VAL'],
        'ddg_exp': [1.2, -0.5, 2.3, 0.8, 1.5],  # Example experimental values
        'temperature': [298] * 5,
        'pH': [7.0] * 5
    }
    
    df = pd.DataFrame(data)
    return df


if __name__ == '__main__':
    # Demo/test
    logging.basicConfig(level=logging.INFO)
    
    # Create sample data
    sample_data = create_sample_protherm_data()
    sample_file = '/tmp/sample_protherm.csv'
    sample_data.to_csv(sample_file, index=False)
    
    # Initialize validator
    validator = ProThermValidator(sample_file)
    
    # Create test set
    test_set = validator.create_test_dataset(n_samples=5)
    print("\nTest Set:")
    print(test_set)
    
    # Simulate some predictions (in practice, these come from Mutalyze)
    predictions = [1.5, -0.3, 2.1, 0.9, 1.7]  # Example predictions
    experimental = list(test_set['ddg_exp'])
    
    # Calculate metrics
    results = validator.calculate_metrics(predictions, experimental)
    
    # Print report
    validator.print_report(results)
    
    # Create plot
    validator.plot_correlation(results, '/tmp/test_correlation.png')
    print(f"\nPlot saved to /tmp/test_correlation.png")
