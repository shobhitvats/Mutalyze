"""
Example Script: Basic Mutalyze Usage
Demonstrates the core functionality without the Streamlit GUI
"""

import sys
from pathlib import Path

# Add core to path
sys.path.insert(0, str(Path(__file__).parent))

from core.pdb_utils import PDBFetcher, PDBAnalyzer
from core.mutation_builder import MutationBuilder
from core.interface_analysis import InterfaceAnalyzer
from core.analysis import StabilityAnalyzer

def main():
    """Run a simple analysis example"""
    
    print("=" * 60)
    print("MUTALYZE - Example Analysis")
    print("=" * 60)
    print()
    
    # Step 1: Fetch PDB structure
    print("Step 1: Fetching PDB structure...")
    pdb_id = "1crn"  # Crambin - small protein for testing
    
    fetcher = PDBFetcher()
    structure = fetcher.fetch_pdb(pdb_id)
    
    if not structure:
        print(f"Failed to fetch {pdb_id}")
        return
    
    print(f"✓ Successfully loaded {pdb_id.upper()}")
    print()
    
    # Step 2: Analyze structure
    print("Step 2: Analyzing structure...")
    metadata = PDBAnalyzer.get_metadata(structure)
    
    print(f"  Chains: {metadata['num_chains']}")
    print(f"  Residues: {metadata['num_residues']}")
    print(f"  Ligands: {metadata['num_ligands']}")
    print(f"  Chain IDs: {', '.join(metadata['chains'])}")
    print()
    
    # Step 3: Enumerate mutations
    print("Step 3: Enumerating mutations...")
    builder = MutationBuilder()
    
    # Just analyze first 5 residues for demonstration
    mutations = builder.enumerate_mutations(
        structure,
        chain_ids=['A'],
        residue_range=(1, 5)
    )
    
    print(f"✓ Enumerated {len(mutations)} possible mutations")
    print()
    
    # Show first few mutations
    print("First 10 mutations:")
    for i, mutation in enumerate(mutations[:10], 1):
        print(f"  {i}. {mutation['mutation_code']}: "
              f"{mutation['wt_resname']} → {mutation['mut_resname']}")
    print()
    
    # Step 4: Interface analysis
    print("Step 4: Analyzing interfaces and binding sites...")
    analyzer = InterfaceAnalyzer()
    
    # Detect ligand proximity
    ligand_prox = analyzer.detect_ligand_proximity(structure)
    print(f"  Found binding sites for {len(ligand_prox)} ligands")
    
    # Classify residue locations
    print()
    print("Residue classifications (first 5):")
    for i in range(1, 6):
        location = analyzer.classify_residue_location(structure, 'A', i)
        print(f"  Residue {i}: {location}")
    print()
    
    # Step 5: Demo analysis (without actual energy calculations)
    print("Step 5: Generating demo stability analysis...")
    
    # Add fake ΔΔG values for demonstration
    import random
    random.seed(42)
    
    for mutation in mutations:
        mutation['ddg'] = random.gauss(0, 1.5)  # Random ΔΔG values
    
    # Analyze results
    summary = StabilityAnalyzer.generate_summary_statistics(mutations)
    
    print(f"  Total mutations analyzed: {summary['analyzed_mutations']}")
    print(f"  Mean ΔΔG: {summary['mean_ddg']:.2f} ± {summary['std_ddg']:.2f} kcal/mol")
    print(f"  Range: {summary['min_ddg']:.2f} to {summary['max_ddg']:.2f} kcal/mol")
    print()
    print(f"  Stabilizing: {summary['stabilizing_count']} ({summary['stabilizing_percent']:.1f}%)")
    print(f"  Neutral: {summary['neutral_count']} ({summary['neutral_percent']:.1f}%)")
    print(f"  Destabilizing: {summary['destabilizing_count']} ({summary['destabilizing_percent']:.1f}%)")
    print()
    
    # Find hotspots
    hotspots = StabilityAnalyzer.find_hotspot_residues(mutations, threshold=2.0)
    
    if hotspots:
        print("Mutation hotspots (|ΔΔG| > 2.0 kcal/mol):")
        for hotspot in hotspots[:5]:
            print(f"  {hotspot['wt_resname']}{hotspot['chain']}{hotspot['resid']}: "
                  f"max |ΔΔG| = {hotspot['max_ddg']:.2f} kcal/mol")
    else:
        print("No significant hotspots found")
    
    print()
    print("=" * 60)
    print("Analysis complete!")
    print()
    print("Note: This example uses random ΔΔG values for demonstration.")
    print("For real energy calculations, use the Streamlit app with OpenMM installed.")
    print()
    print("To run the full GUI: streamlit run app.py")
    print("=" * 60)


if __name__ == "__main__":
    main()
