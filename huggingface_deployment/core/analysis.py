"""
Analysis Module
Interprets ΔΔG values and provides biological significance
"""

import numpy as np
from typing import List, Dict, Tuple, Optional
import logging

logger = logging.getLogger(__name__)


class StabilityAnalyzer:
    """Analyzes and interprets protein stability changes"""
    
    # ΔΔG interpretation thresholds (kcal/mol)
    THRESHOLDS = {
        'highly_stabilizing': -2.0,
        'moderately_stabilizing': -1.0,
        'neutral_lower': -1.0,
        'neutral_upper': 1.0,
        'moderately_destabilizing': 1.0,
        'highly_destabilizing': 2.0
    }
    
    @staticmethod
    def interpret_ddg(ddg: float) -> Dict[str, str]:
        """
        Interpret ΔΔG value
        
        Args:
            ddg: ΔΔG value in kcal/mol
            
        Returns:
            Dictionary with interpretation
        """
        if ddg < StabilityAnalyzer.THRESHOLDS['highly_stabilizing']:
            category = 'highly_stabilizing'
            impact = 'high'
            description = 'Strongly stabilizing mutation'
            color = 'darkblue'
        elif ddg < StabilityAnalyzer.THRESHOLDS['moderately_stabilizing']:
            category = 'moderately_stabilizing'
            impact = 'moderate'
            description = 'Moderately stabilizing mutation'
            color = 'blue'
        elif ddg <= StabilityAnalyzer.THRESHOLDS['neutral_upper']:
            category = 'neutral'
            impact = 'low'
            description = 'Neutral or minor effect'
            color = 'gray'
        elif ddg < StabilityAnalyzer.THRESHOLDS['highly_destabilizing']:
            category = 'moderately_destabilizing'
            impact = 'moderate'
            description = 'Moderately destabilizing mutation'
            color = 'orange'
        else:
            category = 'highly_destabilizing'
            impact = 'high'
            description = 'Strongly destabilizing mutation'
            color = 'red'
        
        return {
            'category': category,
            'impact': impact,
            'description': description,
            'color': color,
            'ddg': ddg
        }
    
    @staticmethod
    def classify_mutations(mutations: List[Dict], 
                          ddg_key: str = 'ddg') -> Dict[str, List[Dict]]:
        """
        Classify mutations by their ΔΔG impact
        
        Args:
            mutations: List of mutation dictionaries with ΔΔG values
            ddg_key: Key name for ΔΔG value in dictionaries
            
        Returns:
            Dictionary of mutations grouped by category
        """
        classified = {
            'highly_stabilizing': [],
            'moderately_stabilizing': [],
            'neutral': [],
            'moderately_destabilizing': [],
            'highly_destabilizing': []
        }
        
        for mutation in mutations:
            if ddg_key not in mutation:
                continue
            
            ddg = mutation[ddg_key]
            interpretation = StabilityAnalyzer.interpret_ddg(ddg)
            
            # Add interpretation to mutation
            mutation['interpretation'] = interpretation
            
            # Classify
            category = interpretation['category']
            classified[category].append(mutation)
        
        return classified
    
    @staticmethod
    def find_hotspot_residues(mutations: List[Dict], 
                             threshold: float = 2.0,
                             ddg_key: str = 'ddg') -> List[Dict]:
        """
        Find residues where mutations have large stability effects (hotspots)
        
        Args:
            mutations: List of mutations
            threshold: Absolute ΔΔG threshold for hotspots
            ddg_key: Key for ΔΔG values
            
        Returns:
            List of hotspot residues
        """
        hotspots = {}
        
        for mutation in mutations:
            if ddg_key not in mutation:
                continue
            
            ddg = abs(mutation[ddg_key])
            if ddg >= threshold:
                key = (mutation.get('chain', ''), mutation.get('resid', 0))
                
                if key not in hotspots:
                    hotspots[key] = {
                        'chain': mutation.get('chain', ''),
                        'resid': mutation.get('resid', 0),
                        'wt_resname': mutation.get('wt_resname', ''),
                        'max_ddg': ddg,
                        'sensitive_mutations': []
                    }
                else:
                    if ddg > hotspots[key]['max_ddg']:
                        hotspots[key]['max_ddg'] = ddg
                
                hotspots[key]['sensitive_mutations'].append(mutation)
        
        # Convert to list and sort by max ΔΔG
        hotspot_list = sorted(hotspots.values(), 
                            key=lambda x: x['max_ddg'], 
                            reverse=True)
        
        logger.info(f"Found {len(hotspot_list)} hotspot residues")
        return hotspot_list
    
    @staticmethod
    def analyze_conservation_stability_correlation(mutations: List[Dict],
                                                   conservation_key: str = 'conservation',
                                                   ddg_key: str = 'ddg') -> Dict:
        """
        Analyze correlation between conservation and stability
        
        Args:
            mutations: List of mutations with conservation and ΔΔG data
            conservation_key: Key for conservation scores
            ddg_key: Key for ΔΔG values
            
        Returns:
            Analysis results dictionary
        """
        conservation_scores = []
        ddg_values = []
        
        for mutation in mutations:
            if conservation_key in mutation and ddg_key in mutation:
                conservation_scores.append(mutation[conservation_key])
                ddg_values.append(abs(mutation[ddg_key]))
        
        if len(conservation_scores) < 2:
            return {
                'correlation': None,
                'n_points': len(conservation_scores),
                'message': 'Insufficient data for correlation analysis'
            }
        
        # Calculate Pearson correlation
        conservation_arr = np.array(conservation_scores)
        ddg_arr = np.array(ddg_values)
        
        correlation = np.corrcoef(conservation_arr, ddg_arr)[0, 1]
        
        # Interpretation
        if abs(correlation) > 0.7:
            strength = 'strong'
        elif abs(correlation) > 0.4:
            strength = 'moderate'
        else:
            strength = 'weak'
        
        if correlation > 0:
            direction = 'positive'
            interpretation = 'Highly conserved positions tend to be more sensitive to mutations'
        else:
            direction = 'negative'
            interpretation = 'Weak or inverse relationship between conservation and mutation sensitivity'
        
        return {
            'correlation': correlation,
            'strength': strength,
            'direction': direction,
            'interpretation': interpretation,
            'n_points': len(conservation_scores)
        }
    
    @staticmethod
    def generate_summary_statistics(mutations: List[Dict], 
                                   ddg_key: str = 'ddg') -> Dict:
        """
        Generate summary statistics for mutation analysis
        
        Args:
            mutations: List of mutations with ΔΔG values
            ddg_key: Key for ΔΔG values
            
        Returns:
            Summary statistics dictionary
        """
        # Filter to get valid ddg values (not None, not inf, not nan)
        ddg_values = []
        for m in mutations:
            if ddg_key in m and m[ddg_key] is not None:
                val = m[ddg_key]
                if np.isfinite(val):
                    ddg_values.append(val)
        
        if not ddg_values:
            return {
                'total_mutations': len(mutations),
                'analyzed_mutations': 0,
                'valid_mutations': 0,
                'mean_ddg': 0.0,
                'median_ddg': 0.0,
                'std_ddg': 0.0,
                'min_ddg': 0.0,
                'max_ddg': 0.0,
                'stabilizing_count': 0,
                'neutral_count': 0,
                'destabilizing_count': 0,
                'highly_destabilizing_count': 0,
                'stabilizing_percent': 0.0,
                'neutral_percent': 0.0,
                'destabilizing_percent': 0.0,
                'error': 'No valid ΔΔG values found'
            }
        
        ddg_array = np.array(ddg_values)
        
        summary = {
            'total_mutations': len(mutations),
            'analyzed_mutations': len([m for m in mutations if ddg_key in m]),
            'valid_mutations': len(ddg_values),
            'mean_ddg': float(np.mean(ddg_array)),
            'median_ddg': float(np.median(ddg_array)),
            'std_ddg': float(np.std(ddg_array)),
            'min_ddg': float(np.min(ddg_array)),
            'max_ddg': float(np.max(ddg_array)),
            'stabilizing_count': int(np.sum(ddg_array < -1.0)),
            'neutral_count': int(np.sum((ddg_array >= -1.0) & (ddg_array <= 1.0))),
            'destabilizing_count': int(np.sum(ddg_array > 1.0)),
            'highly_destabilizing_count': int(np.sum(ddg_array > 2.0))
        }
        
        # Add percentages
        total = summary['valid_mutations']
        summary['stabilizing_percent'] = (summary['stabilizing_count'] / total * 100) if total > 0 else 0
        summary['neutral_percent'] = (summary['neutral_count'] / total * 100) if total > 0 else 0
        summary['destabilizing_percent'] = (summary['destabilizing_count'] / total * 100) if total > 0 else 0
        
        return summary
    
    @staticmethod
    def generate_report(mutations: List[Dict], 
                       conservation_data: Optional[Dict] = None) -> str:
        """
        Generate a comprehensive analysis report
        
        Args:
            mutations: List of analyzed mutations
            conservation_data: Optional conservation analysis results
            
        Returns:
            Formatted report string
        """
        report = []
        report.append("=" * 60)
        report.append("MUTALYZE ANALYSIS REPORT")
        report.append("=" * 60)
        report.append("")
        
        # Summary statistics
        summary = StabilityAnalyzer.generate_summary_statistics(mutations)
        report.append("SUMMARY STATISTICS")
        report.append("-" * 60)
        report.append(f"Total mutations analyzed: {summary.get('analyzed_mutations', 0)}")
        report.append(f"Mean ΔΔG: {summary.get('mean_ddg', 0):.2f} ± {summary.get('std_ddg', 0):.2f} kcal/mol")
        report.append(f"Median ΔΔG: {summary.get('median_ddg', 0):.2f} kcal/mol")
        report.append(f"Range: {summary.get('min_ddg', 0):.2f} to {summary.get('max_ddg', 0):.2f} kcal/mol")
        report.append("")
        report.append(f"Stabilizing mutations: {summary.get('stabilizing_count', 0)} ({summary.get('stabilizing_percent', 0):.1f}%)")
        report.append(f"Neutral mutations: {summary.get('neutral_count', 0)} ({summary.get('neutral_percent', 0):.1f}%)")
        report.append(f"Destabilizing mutations: {summary.get('destabilizing_count', 0)} ({summary.get('destabilizing_percent', 0):.1f}%)")
        report.append(f"Highly destabilizing: {summary.get('highly_destabilizing_count', 0)}")
        report.append("")
        
        # Hotspots
        hotspots = StabilityAnalyzer.find_hotspot_residues(mutations)
        if hotspots:
            report.append("MUTATION HOTSPOTS (|ΔΔG| > 2.0 kcal/mol)")
            report.append("-" * 60)
            for i, hotspot in enumerate(hotspots[:10], 1):  # Top 10
                report.append(f"{i}. {hotspot['wt_resname']}{hotspot['chain']}{hotspot['resid']} "
                            f"(max |ΔΔG| = {hotspot['max_ddg']:.2f} kcal/mol, "
                            f"{len(hotspot['sensitive_mutations'])} sensitive mutations)")
            report.append("")
        
        # Conservation correlation
        if conservation_data:
            corr_analysis = StabilityAnalyzer.analyze_conservation_stability_correlation(mutations)
            if corr_analysis.get('correlation') is not None:
                report.append("CONSERVATION-STABILITY CORRELATION")
                report.append("-" * 60)
                report.append(f"Correlation coefficient: {corr_analysis['correlation']:.3f}")
                report.append(f"Strength: {corr_analysis['strength']}")
                report.append(f"Interpretation: {corr_analysis['interpretation']}")
                report.append("")
        
        report.append("=" * 60)
        
        return "\n".join(report)


def interpret_mutation(ddg: float) -> str:
    """Quick interpretation of a single ΔΔG value"""
    result = StabilityAnalyzer.interpret_ddg(ddg)
    return result['description']
