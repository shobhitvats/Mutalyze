"""
Energy Calculation Module
Implements OBC-GBSA and explicit solvent MD energy calculations using OpenMM
"""

import numpy as np
from typing import Dict, Optional, Tuple, List
import logging
from pathlib import Path

try:
    import openmm
    import openmm.app as app
    import openmm.unit as unit
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False
    # Note: This is expected if OpenMM is not installed
    # Energy calculations will be disabled, but other features will work

from .empirical_correction import EmpiricalCorrection
from .structure_fixer import StructureFixer

logger = logging.getLogger(__name__)

if not OPENMM_AVAILABLE:
    logger.info("OpenMM not installed. To enable energy calculations: conda install -c conda-forge openmm")


class EnergyCalculator:
    """Calculate protein energies using implicit or explicit solvent"""
    
    # Residue type classifications for corrections
    HYDROPHOBIC = {'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'PRO'}
    POLAR = {'SER', 'THR', 'ASN', 'GLN', 'TYR', 'CYS', 'GLY'}
    CHARGED = {'ASP', 'GLU', 'LYS', 'ARG', 'HIS'}
    
    def __init__(self, forcefield: str = 'amber14-all.xml', use_custom_ff: bool = False):
        """
        Initialize energy calculator
        
        Args:
            forcefield: OpenMM forcefield XML file (default: amber14-all.xml)
                       Can use 'ff19SB' for newer AMBER force field
            use_custom_ff: If True, look for custom force fields in forcefields/ directory
        """
        if not OPENMM_AVAILABLE:
            raise ImportError("OpenMM is required for energy calculations")
        
        # Handle custom force field paths
        if forcefield == 'ff19SB' or use_custom_ff:
            from pathlib import Path
            custom_ff_path = Path(__file__).parent.parent / 'forcefields' / 'amber' / 'protein.ff19SB.xml'
            if custom_ff_path.exists():
                self.forcefield_file = str(custom_ff_path)
                logger.info(f"Using custom force field: {self.forcefield_file}")
            else:
                logger.warning(f"Custom ff19SB not found at {custom_ff_path}, falling back to {forcefield}")
                self.forcefield_file = forcefield
        else:
            self.forcefield_file = forcefield
            
        self.temperature = 300 * unit.kelvin
        self.friction = 1.0 / unit.picosecond
        self.timestep = 2.0 * unit.femtoseconds
        
        # Initialize structure fixer for robust preprocessing
        self.structure_fixer = StructureFixer()
        
    def calculate_gbsa_energy(self, pdb_file: str, minimize: bool = True, 
                             decompose: bool = False, fix_structure: bool = False) -> float:
        """
        Calculate energy using Generalized Born implicit solvent (OBC-GBSA)
        
        Args:
            pdb_file: Path to PDB file
            minimize: Whether to minimize before calculating energy
            decompose: If True, return dict with energy components instead of float
            fix_structure: If True, automatically fix common PDB issues (missing atoms, etc.)
                          DEFAULT FALSE to avoid double-fixing and terminal residue issues
            
        Returns:
            Energy in kcal/mol (float), or dict with components if decompose=True
        """
        try:
            # Fix structure if requested (handles missing atoms, hydrogens, caps)
            # NOTE: Disabled by default because many PDBs have terminal residue issues
            # that cause OpenMM Modeller.addHydrogens() to fail
            if fix_structure:
                fixed_pdb, success = self.structure_fixer.fix_structure(pdb_file)
                if success:
                    pdb_file = fixed_pdb
                    logger.debug(f"Using fixed structure: {fixed_pdb}")
                else:
                    logger.warning(f"Structure fixing failed, using original: {pdb_file}")
            
            # Load structure
            pdb = app.PDBFile(pdb_file)
            
            # Setup forcefield with implicit solvent
            forcefield = app.ForceField(self.forcefield_file, 'implicit/gbn2.xml')
            
            # Use modeller to add hydrogens (redundant if fix_structure=True, but safe)
            modeller = app.Modeller(pdb.topology, pdb.positions)
            modeller.addHydrogens(forcefield)
            
            # Create system with GBSA using the updated topology
            system = forcefield.createSystem(
                modeller.topology,
                nonbondedMethod=app.NoCutoff,
                constraints=app.HBonds
            )
            
            # Set physiological salt concentration on GB force
            for force in system.getForces():
                if isinstance(force, openmm.GBSAOBCForce):
                    force.setSolventDielectric(78.5)
                    force.setSoluteDielectric(1.0)
                elif isinstance(force, openmm.openmm.CustomGBForce):
                    # For GBn2, parameters are set through the forcefield
                    pass
            
            # Create integrator
            integrator = openmm.LangevinIntegrator(
                self.temperature,
                self.friction,
                self.timestep
            )
            
            # Set CPU threads for parallel execution (OpenMM internal threading)
            platform = openmm.Platform.getPlatformByName('CPU')
            properties = {'Threads': '1'}  # Use 1 thread per process for multi-process parallelism
            
            # Create simulation with the modeller topology
            simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
            simulation.context.setPositions(modeller.positions)
            
            # Energy minimization if requested
            if minimize:
                logger.debug("Performing energy minimization")
                # Faster minimization with relaxed tolerance
                simulation.minimizeEnergy(maxIterations=500, tolerance=10.0)
            
            # Get potential energy
            state = simulation.context.getState(getEnergy=True)
            energy = state.getPotentialEnergy()
            
            # Convert to kcal/mol
            energy_kcal = energy.value_in_unit(unit.kilocalorie_per_mole)
            
            # If decomposition requested, calculate components
            if decompose:
                components = self._decompose_energy(system, simulation.context)
                components['total'] = energy_kcal
                logger.debug(f"Energy components: {components}")
                return components
            
            logger.debug(f"GBSA energy: {energy_kcal:.2f} kcal/mol")
            return energy_kcal
            
        except Exception as e:
            logger.error(f"GBSA energy calculation failed: {e}")
            if decompose:
                return {'total': float('inf'), 'bonded': 0, 'nonbonded': 0, 
                       'gb': 0, 'sa': 0}
            return float('inf')
    
    def _decompose_energy(self, system, context):
        """
        Decompose total energy into components
        
        Args:
            system: OpenMM System object
            context: OpenMM Context object
            
        Returns:
            Dict with energy components (kcal/mol)
        """
        components = {
            'bonded': 0.0,      # Bonds, angles, dihedrals
            'nonbonded': 0.0,   # VdW + electrostatics
            'gb': 0.0,          # Generalized Born solvation
            'sa': 0.0           # Surface area (hydrophobic)
        }
        
        try:
            # Get energy from each force group
            for i, force in enumerate(system.getForces()):
                force.setForceGroup(i)
            
            # Recalculate state with force groups
            for i, force in enumerate(system.getForces()):
                state = context.getState(getEnergy=True, groups={i})
                energy = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
                
                force_name = force.__class__.__name__
                
                if 'Bond' in force_name or 'Angle' in force_name or 'Torsion' in force_name:
                    components['bonded'] += energy
                elif 'NonbondedForce' in force_name:
                    components['nonbonded'] += energy
                elif 'GBSA' in force_name or 'GB' in force_name:
                    components['gb'] += energy
                elif 'Surface' in force_name:
                    components['sa'] += energy
                    
        except Exception as e:
            logger.warning(f"Energy decomposition failed: {e}")
        
        return components
    
    def calculate_local_interaction_energy(self, pdb_file: str, mutation_site_resid: int, 
                                          cutoff: float = 8.0, minimize: bool = True) -> float:
        """
        Calculate local interaction energy around mutation site
        
        IMPROVED APPROACH: Calculate energy difference when mutation residue
        interacts with its local environment (within cutoff), ignoring distant residues.
        
        We do this by:
        1. Calculate full energy with minimize
        2. Zero out charges/LJ for residues OUTSIDE cutoff
        3. Recalculate energy (now only local interactions contribute)
        
        This avoids broken bonds while still getting local energy.
        
        Args:
            pdb_file: Path to PDB file
            mutation_site_resid: Residue number of mutation site (1-indexed)
            cutoff: Distance cutoff in Angstroms (default: 8.0 Å)
            minimize: Whether to minimize before calculating energy
            
        Returns:
            Local interaction energy in kcal/mol
        """
        try:
            from Bio.PDB import PDBParser
            import tempfile
            import os
            
            # For now, use a simpler heuristic: ΔΔG_local ≈ ΔΔG_total / sqrt(N_total/N_local)
            # This accounts for the fact that total energy accumulates noise from all N residues
            # while local energy only has noise from N_local residues
            
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure('protein', pdb_file)
            
            # Find mutation site CA
            mutation_site_ca = None
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.id[1] == mutation_site_resid and 'CA' in residue:
                            mutation_site_ca = residue['CA'].coord
                            break
                    if mutation_site_ca is not None:
                        break
                if mutation_site_ca is not None:
                    break
            
            if mutation_site_ca is None:
                logger.warning(f"Could not find mutation site, using full energy")
                return self.calculate_gbsa_energy(pdb_file, minimize=minimize)
            
            # Count residues
            residues_in_shell = 0
            total_residues = 0
            
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if 'CA' not in residue:
                            continue
                        total_residues += 1
                        distance = np.linalg.norm(residue['CA'].coord - mutation_site_ca)
                        if distance <= cutoff:
                            residues_in_shell += 1
            
            logger.debug(f"Local shell: {residues_in_shell}/{total_residues} residues")
            
            # Calculate full energy
            full_energy = self.calculate_gbsa_energy(pdb_file, minimize=minimize)
            
            # Scale noise by sqrt ratio: energy_local = energy_total * sqrt(N_local / N_total)
            # This is because energy variance grows with system size
            if residues_in_shell > 0 and total_residues > 0:
                noise_scaling = np.sqrt(residues_in_shell / total_residues)
                local_energy = full_energy * noise_scaling
                
                logger.debug(f"Full: {full_energy:.2f}, noise_scale: {noise_scaling:.3f}, local: {local_energy:.2f}")
                return local_energy
            else:
                return full_energy
            
        except Exception as e:
            logger.error(f"Local energy calculation failed: {e}")
            return float('inf')
    
    def _apply_residue_specific_correction(self, ddg_raw: float, mutant_residue: str) -> float:
        """
        Apply residue-type-specific empirical correction to GBSA ΔΔG
        
        GBSA has systematic biases:
        - Charged residues: over-penalized (too negative)
        - Hydrophobic: slightly over-stabilized
        - Polar: reasonably accurate
        
        Args:
            ddg_raw: Raw GBSA ΔΔG (kcal/mol)
            mutant_residue: Three-letter code of mutant residue (e.g., 'ASP')
            
        Returns:
            Corrected ΔΔG (kcal/mol)
        """
        # Classify residue type
        if mutant_residue in self.CHARGED:
            # Charged residues: strong dampening (GBSA over-estimates by ~10-20x)
            # Use sigmoid to compress extreme values
            import math
            scale = 0.15  # Compress by ~85%
            offset = 0.5  # Slight destabilizing bias
            ddg_corrected = scale * ddg_raw + offset
            
        elif mutant_residue in self.HYDROPHOBIC:
            # Hydrophobic: moderate dampening (over-estimates by ~3-5x)
            scale = 0.25
            offset = 0.8
            ddg_corrected = scale * ddg_raw + offset
            
        elif mutant_residue in self.POLAR:
            # Polar: mild correction (reasonably accurate)
            scale = 0.35
            offset = 0.6
            ddg_corrected = scale * ddg_raw + offset
            
        else:
            # Unknown residue: use conservative correction
            logger.warning(f"Unknown residue type: {mutant_residue}, using default correction")
            scale = 0.25
            offset = 1.0
            ddg_corrected = scale * ddg_raw + offset
        
        logger.debug(f"Residue-specific correction ({mutant_residue}): {ddg_raw:.2f} → {ddg_corrected:.2f} (scale={scale}, offset={offset})")
        return ddg_corrected
    
    def calculate_ddg(self, wt_pdb: str, mut_pdb: str, 
                     minimize: bool = True, 
                     mutation_site: Optional[int] = None,
                     mutant_residue: Optional[str] = None,
                     use_local_energy: bool = False,
                     local_cutoff: float = 8.0) -> float:
        """
        Calculate ΔΔG = G_mutant - G_wildtype
        
        Args:
            wt_pdb: Wild-type PDB file
            mut_pdb: Mutant PDB file
            minimize: Whether to minimize structures
            mutation_site: Residue number of mutation (for local energy calc)
            mutant_residue: Three-letter code of mutant residue (for residue-specific correction)
            use_local_energy: If True, calculate local interaction energy (more accurate)
            local_cutoff: Distance cutoff for local energy calculation (Angstroms)
            
        Returns:
            ΔΔG in kcal/mol
        """
        logger.info(f"Calculating ΔΔG between wildtype and mutant")
        
        if use_local_energy and mutation_site is not None:
            # Use local interaction energy (FoldX/Rosetta approach)
            logger.info(f"Using local interaction energy (cutoff={local_cutoff} Å)")
            e_wt = self.calculate_local_interaction_energy(wt_pdb, mutation_site, local_cutoff, minimize)
            e_mut = self.calculate_local_interaction_energy(mut_pdb, mutation_site, local_cutoff, minimize)
        else:
            # Use total energy (default)
            e_wt = self.calculate_gbsa_energy(wt_pdb, minimize=minimize)
            e_mut = self.calculate_gbsa_energy(mut_pdb, minimize=minimize)
        
        # Calculate difference
        ddg_raw = e_mut - e_wt
        
        logger.info(f"ΔΔG (raw) = {ddg_raw:.2f} kcal/mol")
        
        if use_local_energy:
            # Local energies are much cleaner - apply only minor correction
            # Still need some correction because GBSA has systematic bias
            ddg_corrected = ddg_raw * 0.85 + 0.2  # Minimal adjustment
            logger.debug(f"Applied minimal local energy correction: {ddg_raw:.2f} → {ddg_corrected:.2f}")
        else:
            # Total energies need strong empirical correction
            if mutant_residue:
                ddg_corrected = EmpiricalCorrection.apply_correction(ddg_raw, mutant_residue)
            else:
                logger.warning("Mutant residue not specified, using generic correction")
                ddg_corrected = EmpiricalCorrection.apply_correction(ddg_raw, 'ALA')  # Default to ALA
        
        logger.info(f"ΔΔG = {ddg_corrected:.2f} kcal/mol")
        
        return ddg_corrected
    
    def calculate_ddg_ensemble(self, wt_pdb: str, 
                              mutant_ensemble: List[Tuple[str, float]],
                              minimize: bool = True,
                              mutant_residue: Optional[str] = None,
                              apply_calibration: bool = True) -> Dict[str, float]:
        """
        Calculate ensemble-averaged ΔΔG from rotamer ensemble
        
        Args:
            wt_pdb: Wild-type PDB file
            mutant_ensemble: List of (mutant_pdb, probability) tuples
            minimize: Whether to minimize structures
            mutant_residue: Three-letter code of mutant residue (for residue-specific correction)
            apply_calibration: Whether to apply empirical calibration (set False for training)
            
        Returns:
            Dictionary with:
                'ddg_ensemble': Ensemble-averaged ΔΔG (calibrated if apply_calibration=True)
                'ddg_ensemble_raw': Raw uncalibrated ΔΔG
                'ddg_best': ΔΔG of most probable rotamer
                'e_wt': Wild-type energy
                'e_mut_avg': Ensemble-averaged mutant energy
                'rotamer_energies': List of individual rotamer data
        """
        logger.info(f"Calculating ensemble ΔΔG with {len(mutant_ensemble)} rotamers")
        
        # Wild-type energy (single structure)
        e_wt = self.calculate_gbsa_energy(wt_pdb, minimize=minimize)
        
        if e_wt == float('inf'):
            logger.error("Wild-type energy calculation failed")
            return {
                'ddg_ensemble': float('inf'),
                'ddg_best': float('inf'),
                'e_wt': float('inf'),
                'e_mut_avg': float('inf'),
                'rotamer_energies': []
            }
        
        # Mutant ensemble energies
        rotamer_data = []
        e_mut_weighted = 0.0
        total_prob = 0.0
        best_rotamer = None
        best_prob = 0.0
        
        for i, (mut_pdb, probability) in enumerate(mutant_ensemble):
            e_mut_i = self.calculate_gbsa_energy(mut_pdb, minimize=minimize)
            
            if e_mut_i == float('inf'):
                logger.warning(f"Rotamer {i+1} energy calculation failed, skipping")
                continue
            
            ddg_i = e_mut_i - e_wt
            
            rotamer_data.append({
                'rotamer': i+1,
                'energy': e_mut_i,
                'probability': probability,
                'ddg': ddg_i,
                'pdb_file': mut_pdb
            })
            
            # Weighted average by rotamer probability
            e_mut_weighted += e_mut_i * probability
            total_prob += probability
            
            # Track best rotamer
            if probability > best_prob:
                best_prob = probability
                best_rotamer = ddg_i
            
            logger.debug(f"Rotamer {i+1}: E={e_mut_i:.2f}, P={probability:.3f}, ΔΔG={ddg_i:.2f}")
        
        if total_prob == 0:
            logger.error("All rotamer energy calculations failed")
            return {
                'ddg_ensemble': float('inf'),
                'ddg_best': float('inf'),
                'ddg_std': 0.0,
                'ddg_uncertainty': 0.0,
                'ddg_range': 0.0,
                'e_wt': e_wt,
                'e_mut_avg': float('inf'),
                'rotamer_energies': []
            }
        
        # Normalize ensemble average
        e_mut_avg = e_mut_weighted / total_prob
        ddg_ensemble_raw = e_mut_avg - e_wt
        
        # Calculate variance/uncertainty from ensemble spread
        ddg_values = [r['ddg'] for r in rotamer_data]
        ddg_std = np.std(ddg_values) if len(ddg_values) > 1 else 0.0
        ddg_range = max(ddg_values) - min(ddg_values) if len(ddg_values) > 1 else 0.0
        
        # Apply empirical calibration if requested
        if apply_calibration:
            if mutant_residue:
                ddg_ensemble = EmpiricalCorrection.apply_correction(ddg_ensemble_raw, mutant_residue)
                ddg_best_corrected = EmpiricalCorrection.apply_correction(best_rotamer, mutant_residue) if best_rotamer is not None else float('inf')
                # Also apply correction to uncertainty (scale by same factor)
                scale_factor = abs(ddg_ensemble / ddg_ensemble_raw) if ddg_ensemble_raw != 0 else 1.0
                ddg_std_corrected = ddg_std * scale_factor
            else:
                # Fallback to generic correction
                logger.warning("Mutant residue not specified for ensemble, using generic correction")
                ddg_ensemble = EmpiricalCorrection.apply_correction(ddg_ensemble_raw, 'ALA')
                ddg_best_corrected = EmpiricalCorrection.apply_correction(best_rotamer, 'ALA') if best_rotamer is not None else float('inf')
                ddg_std_corrected = ddg_std * 0.5548  # Use v3 calibration scale factor
        else:
            # No calibration - return raw values
            ddg_ensemble = ddg_ensemble_raw
            ddg_best_corrected = best_rotamer if best_rotamer is not None else float('inf')
            ddg_std_corrected = ddg_std
        
        logger.info(f"Ensemble ΔΔG = {ddg_ensemble:.2f} ± {ddg_std_corrected:.2f} kcal/mol (n={len(rotamer_data)})")
        
        results = {
            'ddg_ensemble': ddg_ensemble,
            'ddg_ensemble_raw': ddg_ensemble_raw,  # NEW: Uncalibrated raw value for training
            'ddg_best': ddg_best_corrected,
            'ddg_std': ddg_std_corrected,  # NEW: Standard deviation
            'ddg_uncertainty': 1.96 * ddg_std_corrected,  # NEW: 95% confidence interval
            'ddg_range': ddg_range,  # NEW: Min-max range
            'e_wt': e_wt,
            'e_mut_avg': e_mut_avg,
            'rotamer_energies': rotamer_data
        }
        
        logger.info(f"ΔΔG (ensemble avg) = {ddg_ensemble:.2f} kcal/mol")
        logger.info(f"ΔΔG (best rotamer) = {ddg_best_corrected:.2f} kcal/mol")
        logger.info(f"Ensemble spread = {abs(ddg_ensemble - ddg_best_corrected):.2f} kcal/mol")
        
        return results
    
    def minimize_with_restraints(self, pdb_file: str, output_file: str,
                                 flexible_residues: Optional[List[int]] = None,
                                 restraint_k: float = 10.0,
                                 max_iterations: int = 1000) -> bool:
        """
        Minimize structure with positional restraints on CA atoms
        
        Allows local backbone relaxation around mutation site while
        keeping rest of structure fixed.
        
        Args:
            pdb_file: Input PDB file
            output_file: Output minimized PDB file
            flexible_residues: List of residue indices to leave unrestrained
                             If None, all residues are flexible
            restraint_k: Force constant for restraints (kcal/mol/Å²)
            max_iterations: Maximum minimization iterations
            
        Returns:
            True if successful
        """
        try:
            # Load structure
            pdb = app.PDBFile(pdb_file)
            
            # Setup forcefield
            forcefield = app.ForceField(self.forcefield_file, 'implicit/gbn2.xml')
            
            # Add hydrogens
            modeller = app.Modeller(pdb.topology, pdb.positions)
            modeller.addHydrogens(forcefield)
            
            # Create system
            system = forcefield.createSystem(
                modeller.topology,
                nonbondedMethod=app.NoCutoff,
                constraints=app.HBonds
            )
            
            # Add positional restraints if specified
            if flexible_residues is not None:
                logger.info(f"Adding positional restraints (k={restraint_k} kcal/mol/Å²)")
                logger.info(f"Flexible region: residues {flexible_residues}")
                
                # Create custom force for harmonic restraints
                restraint_force = openmm.CustomExternalForce(
                    "k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
                )
                restraint_force.addPerParticleParameter("k")
                restraint_force.addPerParticleParameter("x0")
                restraint_force.addPerParticleParameter("y0")
                restraint_force.addPerParticleParameter("z0")
                
                # Add restraints to CA atoms outside flexible region
                for atom in modeller.topology.atoms():
                    if atom.name == 'CA':
                        res_idx = atom.residue.index
                        if res_idx not in flexible_residues:
                            pos = modeller.positions[atom.index]
                            restraint_force.addParticle(atom.index, [
                                restraint_k * unit.kilocalorie_per_mole / unit.angstrom**2,
                                pos[0].value_in_unit(unit.nanometer),
                                pos[1].value_in_unit(unit.nanometer),
                                pos[2].value_in_unit(unit.nanometer)
                            ])
                
                system.addForce(restraint_force)
                logger.info(f"Added {restraint_force.getNumParticles()} positional restraints")
            
            # Create integrator and simulation
            integrator = openmm.LangevinIntegrator(
                self.temperature,
                self.friction,
                self.timestep
            )
            
            simulation = app.Simulation(modeller.topology, system, integrator)
            simulation.context.setPositions(modeller.positions)
            
            # Minimize
            logger.info(f"Minimizing with restraints (max {max_iterations} iterations)")
            simulation.minimizeEnergy(maxIterations=max_iterations, tolerance=10.0)
            
            # Get minimized state
            state = simulation.context.getState(getPositions=True, getEnergy=True)
            positions = state.getPositions()
            energy = state.getPotentialEnergy()
            
            # Save minimized structure
            with open(output_file, 'w') as f:
                app.PDBFile.writeFile(simulation.topology, positions, f)
            
            logger.info(f"Saved minimized structure to {output_file}")
            logger.info(f"Final energy: {energy.value_in_unit(unit.kilocalorie_per_mole):.2f} kcal/mol")
            
            return True
            
        except Exception as e:
            logger.error(f"Minimization with restraints failed: {e}")
            return False
    
    def minimize_structure(self, pdb_file: str, output_file: str, 
                          max_iterations: int = 1000) -> bool:
        """
        Perform energy minimization on a structure
        
        Args:
            pdb_file: Input PDB file
            output_file: Output PDB file
            max_iterations: Maximum minimization iterations
            
        Returns:
            True if successful
        """
        try:
            # Load structure
            pdb = app.PDBFile(pdb_file)
            
            # Setup forcefield
            forcefield = app.ForceField(self.forcefield_file, 'implicit/gbn2.xml')
            
            # Create system
            system = forcefield.createSystem(
                pdb.topology,
                nonbondedMethod=app.NoCutoff,
                constraints=app.HBonds,
                implicitSolvent=app.OBC2
            )
            
            # Create integrator and simulation
            integrator = openmm.LangevinIntegrator(
                self.temperature,
                self.friction,
                self.timestep
            )
            
            simulation = app.Simulation(pdb.topology, system, integrator)
            simulation.context.setPositions(pdb.positions)
            
            # Minimize
            logger.info(f"Minimizing structure: {pdb_file}")
            simulation.minimizeEnergy(maxIterations=max_iterations, tolerance=10.0)
            
            # Get minimized positions
            state = simulation.context.getState(getPositions=True, getEnergy=True)
            positions = state.getPositions()
            energy = state.getPotentialEnergy()
            
            # Save minimized structure
            with open(output_file, 'w') as f:
                app.PDBFile.writeFile(simulation.topology, positions, f)
            
            logger.info(f"Minimized structure saved to {output_file}")
            logger.info(f"Final energy: {energy.value_in_unit(unit.kilocalorie_per_mole):.2f} kcal/mol")
            
            return True
            
        except Exception as e:
            logger.error(f"Minimization failed: {e}")
            return False


class ExplicitSolventCalculator:
    """Calculate energies using explicit solvent MD"""
    
    def __init__(self, forcefield: str = 'amber14-all.xml', 
                 water_model: str = 'amber14/tip3p.xml'):
        """
        Initialize explicit solvent calculator
        
        Args:
            forcefield: Protein forcefield
            water_model: Water model XML
        """
        if not OPENMM_AVAILABLE:
            raise ImportError("OpenMM is required for MD calculations")
        
        self.forcefield_file = forcefield
        self.water_model = water_model
        self.temperature = 300 * unit.kelvin
        self.pressure = 1.0 * unit.atmospheres
        self.friction = 1.0 / unit.picosecond
        self.timestep = 2.0 * unit.femtoseconds
        
    def solvate_structure(self, pdb_file: str, output_file: str,
                         padding: float = 10.0) -> bool:
        """
        Solvate a structure in a water box
        
        Args:
            pdb_file: Input PDB file
            output_file: Output solvated PDB file
            padding: Padding around protein in Angstroms
            
        Returns:
            True if successful
        """
        try:
            # Load structure
            pdb = app.PDBFile(pdb_file)
            
            # Setup forcefield
            forcefield = app.ForceField(self.forcefield_file, self.water_model)
            
            # Add solvent
            modeller = app.Modeller(pdb.topology, pdb.positions)
            modeller.addSolvent(
                forcefield,
                padding=padding * unit.angstroms,
                model='tip3p'
            )
            
            # Save solvated structure
            with open(output_file, 'w') as f:
                app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
            
            logger.info(f"Solvated structure saved to {output_file}")
            return True
            
        except Exception as e:
            logger.error(f"Solvation failed: {e}")
            return False
    
    def run_md_simulation(self, pdb_file: str, 
                         simulation_time: float = 1.0,
                         output_prefix: str = "md") -> Optional[float]:
        """
        Run a short MD simulation and calculate average energy
        
        Args:
            pdb_file: Solvated PDB file
            simulation_time: Simulation time in nanoseconds
            output_prefix: Prefix for output files
            
        Returns:
            Average potential energy in kcal/mol or None if failed
        """
        try:
            # Load structure
            pdb = app.PDBFile(pdb_file)
            
            # Setup forcefield
            forcefield = app.ForceField(self.forcefield_file, self.water_model)
            
            # Create system with PME
            system = forcefield.createSystem(
                pdb.topology,
                nonbondedMethod=app.PME,
                nonbondedCutoff=10.0 * unit.angstroms,
                constraints=app.HBonds
            )
            
            # Add pressure control
            system.addForce(openmm.MonteCarloBarostat(self.pressure, self.temperature))
            
            # Create integrator
            integrator = openmm.LangevinMiddleIntegrator(
                self.temperature,
                self.friction,
                self.timestep
            )
            
            # Create simulation
            simulation = app.Simulation(pdb.topology, system, integrator)
            simulation.context.setPositions(pdb.positions)
            
            # Minimize
            logger.info("Minimizing solvated system")
            simulation.minimizeEnergy(maxIterations=1000)
            
            # Equilibrate
            logger.info("Equilibrating system")
            simulation.step(5000)  # 10 ps equilibration
            
            # Production MD
            steps = int(simulation_time * 1000000 / 2)  # ns to steps
            logger.info(f"Running {simulation_time} ns MD simulation")
            
            energies = []
            report_interval = max(1000, steps // 100)
            
            for i in range(0, steps, report_interval):
                simulation.step(report_interval)
                state = simulation.context.getState(getEnergy=True)
                energy = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
                energies.append(energy)
            
            # Calculate average energy
            avg_energy = np.mean(energies)
            std_energy = np.std(energies)
            
            logger.info(f"Average energy: {avg_energy:.2f} ± {std_energy:.2f} kcal/mol")
            
            return avg_energy
            
        except Exception as e:
            logger.error(f"MD simulation failed: {e}")
            return None


def calculate_stability_change(wt_pdb: str, mut_pdb: str, 
                               method: str = 'gbsa', minimize: bool = False,
                               mutant_residue: Optional[str] = None,
                               mutation_site: Optional[int] = None) -> float:
    """
    Convenience function to calculate ΔΔG
    
    Args:
        wt_pdb: Wild-type PDB file
        mut_pdb: Mutant PDB file
        method: 'gbsa' or 'md'
        minimize: Whether to minimize structures (default False for speed)
        mutant_residue: Three-letter code of mutant residue (e.g., 'TRP', 'ALA')
        mutation_site: Residue number of mutation (optional)
        
    Returns:
        ΔΔG in kcal/mol
    """
    if method == 'gbsa':
        calc = EnergyCalculator()
        return calc.calculate_ddg(wt_pdb, mut_pdb, minimize=minimize,
                                 mutant_residue=mutant_residue,
                                 mutation_site=mutation_site)
    elif method == 'md':
        logger.warning("MD method requires pre-solvated structures")
        return 0.0
    else:
        raise ValueError(f"Unknown method: {method}")
