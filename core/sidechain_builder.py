"""
Sidechain Builder Module
Builds missing sidechain atoms using ideal geometry
"""

import numpy as np
from Bio.PDB import Atom, Vector
import logging

logger = logging.getLogger(__name__)


def cross_product(v1, v2):
    """Calculate cross product of two Bio.PDB Vectors"""
    v1_np = np.array([v1[0], v1[1], v1[2]])
    v2_np = np.array([v2[0], v2[1], v2[2]])
    result_np = np.cross(v1_np, v2_np)
    return Vector(result_np[0], result_np[1], result_np[2])


class SidechainBuilder:
    """Builds sidechain atoms using ideal geometry"""
    
    # Standard bond lengths (Angstroms)
    BOND_LENGTHS = {
        'C-C': 1.52,
        'C-N': 1.47,
        'C-O': 1.43,
        'C-S': 1.81,
        'N-C': 1.47,
        'S-C': 1.81,
    }
    
    # Standard bond angles (degrees)
    BOND_ANGLES = {
        'tetrahedral': 109.47,
        'planar': 120.0,
    }
    
    # Chi angle defaults (degrees) - use most common rotamer
    CHI_DEFAULTS = {
        'SER': [-60],
        'CYS': [-60],
        'THR': [-60],
        'VAL': [180],
        'ILE': [-60, 180],
        'LEU': [-60, 180],
        'ASP': [-60, -60],
        'ASN': [-60, -60],
        'GLU': [-60, 180, -60],
        'GLN': [-60, 180, -60],
        'MET': [-60, 180, -60],
        'LYS': [-60, 180, -60, 180],
        'ARG': [-60, 180, -60, 180],
        'HIS': [-60, -90],
        'PHE': [-60, 90],
        'TYR': [-60, 90],
        'TRP': [-60, -90],
        'PRO': [30],
    }
    
    @staticmethod
    def build_sidechain(residue, resname):
        """
        Build missing sidechain atoms for a residue
        
        Args:
            residue: Biopython Residue object
            resname: Target residue name (3-letter code)
        """
        # Get backbone atoms
        if not (residue.has_id('N') and residue.has_id('CA') and residue.has_id('C')):
            logger.warning(f"Residue missing backbone atoms, cannot build sidechain")
            return
        
        n = residue['N'].get_vector()
        ca = residue['CA'].get_vector()
        c = residue['C'].get_vector()
        
        # Build CB first (except for GLY)
        if resname != 'GLY' and not residue.has_id('CB'):
            cb_vec = SidechainBuilder._build_cb(n, ca, c)
            if cb_vec is not None:
                # Convert Vector to numpy array
                cb_coord = np.array([cb_vec[0], cb_vec[1], cb_vec[2]])
                cb_atom = Atom.Atom(
                    name='CB',
                    coord=cb_coord,
                    bfactor=residue['CA'].bfactor,
                    occupancy=1.0,
                    altloc=' ',
                    fullname=' CB ',
                    serial_number=residue['CA'].serial_number + 1,
                    element='C'
                )
                residue.add(cb_atom)
        
        # Build rest of sidechain based on residue type
        if resname == 'SER' and not residue.has_id('OG'):
            SidechainBuilder._build_ser(residue, n, ca, c)
        elif resname == 'CYS' and not residue.has_id('SG'):
            SidechainBuilder._build_cys(residue, n, ca, c)
        elif resname == 'THR' and not (residue.has_id('OG1') and residue.has_id('CG2')):
            SidechainBuilder._build_thr(residue, n, ca, c)
        elif resname == 'VAL' and not (residue.has_id('CG1') and residue.has_id('CG2')):
            SidechainBuilder._build_val(residue, n, ca, c)
        elif resname == 'ILE' and not residue.has_id('CD1'):
            SidechainBuilder._build_ile(residue, n, ca, c)
        elif resname == 'LEU' and not (residue.has_id('CD1') and residue.has_id('CD2')):
            SidechainBuilder._build_leu(residue, n, ca, c)
        elif resname == 'ASP' and not (residue.has_id('CG') and residue.has_id('OD1')):
            SidechainBuilder._build_asp(residue, n, ca, c)
        elif resname == 'ASN' and not (residue.has_id('CG') and residue.has_id('OD1')):
            SidechainBuilder._build_asn(residue, n, ca, c)
        elif resname == 'GLU' and not residue.has_id('CD'):
            SidechainBuilder._build_glu(residue, n, ca, c)
        elif resname == 'GLN' and not residue.has_id('CD'):
            SidechainBuilder._build_gln(residue, n, ca, c)
        elif resname == 'MET' and not residue.has_id('SD'):
            SidechainBuilder._build_met(residue, n, ca, c)
        elif resname == 'LYS' and not residue.has_id('NZ'):
            SidechainBuilder._build_lys(residue, n, ca, c)
        elif resname == 'ARG' and not residue.has_id('CZ'):
            SidechainBuilder._build_arg(residue, n, ca, c)
        elif resname == 'HIS' and not residue.has_id('CE1'):
            SidechainBuilder._build_his(residue, n, ca, c)
        elif resname == 'PHE' and not residue.has_id('CZ'):
            SidechainBuilder._build_phe(residue, n, ca, c)
        elif resname == 'TYR' and not residue.has_id('OH'):
            SidechainBuilder._build_tyr(residue, n, ca, c)
        elif resname == 'TRP' and not residue.has_id('CH2'):
            SidechainBuilder._build_trp(residue, n, ca, c)
        elif resname == 'PRO' and not (residue.has_id('CG') and residue.has_id('CD')):
            SidechainBuilder._build_pro(residue, n, ca, c)
    
    @staticmethod
    def _build_cb(n, ca, c):
        """Build CB atom using ideal tetrahedral geometry"""
        try:
            # Convert to numpy arrays for easier vector math
            n_arr = np.array([n[0], n[1], n[2]])
            ca_arr = np.array([ca[0], ca[1], ca[2]])
            c_arr = np.array([c[0], c[1], c[2]])
            
            # Vectors from CA
            v1 = n_arr - ca_arr
            v1 = v1 / np.linalg.norm(v1)
            
            v2 = c_arr - ca_arr
            v2 = v2 / np.linalg.norm(v2)
            
            # CB direction (opposite to average of N and C directions)
            v_avg = v1 + v2
            v_avg = v_avg / np.linalg.norm(v_avg)
            
            # Cross product to get perpendicular vector
            v_perp = np.cross(v1, v2)
            v_perp = v_perp / np.linalg.norm(v_perp)
            
            # Tetrahedral angle from N-CA-C plane
            cb_direction = -0.55 * v_avg + 0.83 * v_perp
            cb_direction = cb_direction / np.linalg.norm(cb_direction)
            cb_arr = ca_arr + 1.54 * cb_direction
            
            # Convert back to Vector
            return Vector(cb_arr[0], cb_arr[1], cb_arr[2])
        except Exception as e:
            logger.debug(f"Failed to build CB: {e}")
            return None
    
    @staticmethod
    def _add_atom(residue, name, coord, element, serial_offset=10):
        """Helper to add an atom to residue"""
        # Ensure coord is a numpy array (convert Vector or list)
        if not isinstance(coord, np.ndarray):
            coord = np.array([coord[0], coord[1], coord[2]])
        
        atom = Atom.Atom(
            name=name,
            coord=coord,  # np.array is fine for Biopython
            bfactor=residue['CA'].bfactor,
            occupancy=1.0,
            altloc=' ',
            fullname=f' {name:<3s}',
            serial_number=residue['CA'].serial_number + serial_offset,
            element=element
        )
        residue.add(atom)
    
    @staticmethod
    def _build_along_bond(origin, direction, length):
        """Build atom along a bond direction"""
        return origin + length * direction.normalized()
    
    @staticmethod
    def _build_ser(residue, n, ca, c):
        """Build SER sidechain: CB-OG"""
        if not residue.has_id('CB'):
            return
        cb = residue['CB'].get_vector()
        ca_vec = residue['CA'].get_vector()
        
        # Convert to numpy
        cb_arr = np.array([cb[0], cb[1], cb[2]])
        ca_arr = np.array([ca_vec[0], ca_vec[1], ca_vec[2]])
        
        # OG along CB-CA direction, bond length 1.43 A
        direction = cb_arr - ca_arr
        direction = direction / np.linalg.norm(direction)
        og_arr = cb_arr + 1.43 * direction
        
        SidechainBuilder._add_atom(residue, 'OG', og_arr, 'O')
    
    @staticmethod
    def _build_cys(residue, n, ca, c):
        """Build CYS sidechain: CB-SG"""
        if not residue.has_id('CB'):
            return
        cb = residue['CB'].get_vector()
        ca_vec = residue['CA'].get_vector()
        
        # Convert to numpy
        cb_arr = np.array([cb[0], cb[1], cb[2]])
        ca_arr = np.array([ca_vec[0], ca_vec[1], ca_vec[2]])
        
        # SG along CB-CA extended, bond length 1.81 A
        direction = cb_arr - ca_arr
        direction = direction / np.linalg.norm(direction)
        sg_arr = cb_arr + 1.81 * direction
        
        SidechainBuilder._add_atom(residue, 'SG', sg_arr, 'S')
    
    @staticmethod
    @staticmethod
    def _build_thr(residue, n, ca, c):
        """Build THR sidechain: CB-OG1, CB-CG2"""
        if not residue.has_id('CB'):
            return
        cb = residue['CB'].get_vector()
        ca_vec = residue['CA'].get_vector()
        c_vec = residue['C'].get_vector()
        
        cb_arr = np.array([cb[0], cb[1], cb[2]])
        ca_arr = np.array([ca_vec[0], ca_vec[1], ca_vec[2]])
        c_arr = np.array([c_vec[0], c_vec[1], c_vec[2]])
        
        direction = cb_arr - ca_arr
        direction = direction / np.linalg.norm(direction)
        
        # OG1
        og1_arr = cb_arr + 1.43 * direction
        SidechainBuilder._add_atom(residue, 'OG1', og1_arr, 'O', 10)
        
        # CG2 - branched from CB
        perp = np.cross(c_arr - ca_arr, cb_arr - ca_arr)
        perp = perp / np.linalg.norm(perp)
        cg2_arr = cb_arr + 1.52 * perp
        SidechainBuilder._add_atom(residue, 'CG2', cg2_arr, 'C', 11)
    
    @staticmethod
    def _build_val(residue, n, ca, c):
        """Build VAL sidechain: CB-CG1, CB-CG2"""
        if not residue.has_id('CB'):
            return
        cb = residue['CB'].get_vector()
        ca_vec = residue['CA'].get_vector()
        c_vec = residue['C'].get_vector()
        
        cb_arr = np.array([cb[0], cb[1], cb[2]])
        ca_arr = np.array([ca_vec[0], ca_vec[1], ca_vec[2]])
        c_arr = np.array([c_vec[0], c_vec[1], c_vec[2]])
        
        v1 = cb_arr - ca_arr
        v1 = v1 / np.linalg.norm(v1)
        
        v2 = np.cross(c_arr - ca_arr, v1)
        v2 = v2 / np.linalg.norm(v2)
        
        # Two methyl groups branching from CB
        cg1_arr = cb_arr + 1.52 * (v1 + 0.5 * v2) / np.linalg.norm(v1 + 0.5 * v2)
        cg2_arr = cb_arr + 1.52 * (v1 - 0.5 * v2) / np.linalg.norm(v1 - 0.5 * v2)
        
        SidechainBuilder._add_atom(residue, 'CG1', cg1_arr, 'C', 10)
        SidechainBuilder._add_atom(residue, 'CG2', cg2_arr, 'C', 11)
    
    @staticmethod
    def _build_asp(residue, n, ca, c):
        """Build ASP sidechain: CB-CG-OD1/OD2"""
        if not residue.has_id('CB'):
            return
        cb = residue['CB'].get_vector()
        ca_vec = residue['CA'].get_vector()
        c_vec = residue['C'].get_vector()
        
        cb_arr = np.array([cb[0], cb[1], cb[2]])
        ca_arr = np.array([ca_vec[0], ca_vec[1], ca_vec[2]])
        c_arr = np.array([c_vec[0], c_vec[1], c_vec[2]])
        
        direction = cb_arr - ca_arr
        direction = direction / np.linalg.norm(direction)
        
        # CG
        cg_arr = cb_arr + 1.52 * direction
        SidechainBuilder._add_atom(residue, 'CG', cg_arr, 'C', 10)
        
        # OD1, OD2 - carboxyl group
        perp = np.cross(c_arr - ca_arr, cb_arr - ca_arr)
        perp = perp / np.linalg.norm(perp)
        
        od1_dir = direction + 0.3 * perp
        od1_arr = cg_arr + 1.25 * od1_dir / np.linalg.norm(od1_dir)
        od2_dir = direction - 0.3 * perp
        od2_arr = cg_arr + 1.25 * od2_dir / np.linalg.norm(od2_dir)
        
        SidechainBuilder._add_atom(residue, 'OD1', od1_arr, 'O', 11)
        SidechainBuilder._add_atom(residue, 'OD2', od2_arr, 'O', 12)
    
    @staticmethod
    def _build_asn(residue, n, ca, c):
        """Build ASN sidechain: CB-CG-OD1/ND2"""
        if not residue.has_id('CB'):
            return
        cb = residue['CB'].get_vector()
        ca_vec = residue['CA'].get_vector()
        c_vec = residue['C'].get_vector()
        
        cb_arr = np.array([cb[0], cb[1], cb[2]])
        ca_arr = np.array([ca_vec[0], ca_vec[1], ca_vec[2]])
        c_arr = np.array([c_vec[0], c_vec[1], c_vec[2]])
        
        direction = cb_arr - ca_arr
        direction = direction / np.linalg.norm(direction)
        
        # CG
        cg_arr = cb_arr + 1.52 * direction
        SidechainBuilder._add_atom(residue, 'CG', cg_arr, 'C', 10)
        
        # OD1, ND2 - amide group
        perp = np.cross(c_arr - ca_arr, cb_arr - ca_arr)
        perp = perp / np.linalg.norm(perp)
        
        od1_dir = direction + 0.3 * perp
        od1_arr = cg_arr + 1.23 * od1_dir / np.linalg.norm(od1_dir)
        nd2_dir = direction - 0.3 * perp
        nd2_arr = cg_arr + 1.33 * nd2_dir / np.linalg.norm(nd2_dir)
        
        SidechainBuilder._add_atom(residue, 'OD1', od1_arr, 'O', 11)
        SidechainBuilder._add_atom(residue, 'ND2', nd2_arr, 'N', 12)
    
    @staticmethod
    def _build_glu(residue, n, ca, c):
        """Build GLU sidechain: CB-CG-CD-OE1/OE2"""
        if not residue.has_id('CB'):
            return
        cb = residue['CB'].get_vector()
        ca_vec = residue['CA'].get_vector()
        c_vec = residue['C'].get_vector()
        
        cb_arr = np.array([cb[0], cb[1], cb[2]])
        ca_arr = np.array([ca_vec[0], ca_vec[1], ca_vec[2]])
        c_arr = np.array([c_vec[0], c_vec[1], c_vec[2]])
        
        direction = cb_arr - ca_arr
        direction = direction / np.linalg.norm(direction)
        
        # CG
        cg_arr = cb_arr + 1.52 * direction
        SidechainBuilder._add_atom(residue, 'CG', cg_arr, 'C', 10)
        
        # CD
        cd_arr = cg_arr + 1.52 * direction
        SidechainBuilder._add_atom(residue, 'CD', cd_arr, 'C', 11)
        
        # OE1, OE2
        perp = np.cross(c_arr - ca_arr, cb_arr - ca_arr)
        perp = perp / np.linalg.norm(perp)
        
        oe1_dir = direction + 0.3 * perp
        oe1_arr = cd_arr + 1.25 * oe1_dir / np.linalg.norm(oe1_dir)
        oe2_dir = direction - 0.3 * perp
        oe2_arr = cd_arr + 1.25 * oe2_dir / np.linalg.norm(oe2_dir)
        
        SidechainBuilder._add_atom(residue, 'OE1', oe1_arr, 'O', 12)
        SidechainBuilder._add_atom(residue, 'OE2', oe2_arr, 'O', 13)
    
    @staticmethod
    def _build_gln(residue, n, ca, c):
        """Build GLN sidechain: CB-CG-CD-OE1/NE2"""
        if not residue.has_id('CB'):
            return
        cb = residue['CB'].get_vector()
        ca_vec = residue['CA'].get_vector()
        c_vec = residue['C'].get_vector()
        
        cb_arr = np.array([cb[0], cb[1], cb[2]])
        ca_arr = np.array([ca_vec[0], ca_vec[1], ca_vec[2]])
        c_arr = np.array([c_vec[0], c_vec[1], c_vec[2]])
        
        direction = cb_arr - ca_arr
        direction = direction / np.linalg.norm(direction)
        
        # CG
        cg_arr = cb_arr + 1.52 * direction
        SidechainBuilder._add_atom(residue, 'CG', cg_arr, 'C', 10)
        
        # CD
        cd_arr = cg_arr + 1.52 * direction
        SidechainBuilder._add_atom(residue, 'CD', cd_arr, 'C', 11)
        
        # OE1, NE2
        perp = np.cross(c_arr - ca_arr, cb_arr - ca_arr)
        perp = perp / np.linalg.norm(perp)
        
        oe1_dir = direction + 0.3 * perp
        oe1_arr = cd_arr + 1.23 * oe1_dir / np.linalg.norm(oe1_dir)
        ne2_dir = direction - 0.3 * perp
        ne2_arr = cd_arr + 1.33 * ne2_dir / np.linalg.norm(ne2_dir)
        
        SidechainBuilder._add_atom(residue, 'OE1', oe1_arr, 'O', 12)
        SidechainBuilder._add_atom(residue, 'NE2', ne2_arr, 'N', 13)
    
    @staticmethod
    def _build_met(residue, n, ca, c):
        """Build MET sidechain: CB-CG-SD-CE"""
        if not residue.has_id('CB'):
            return
        cb = residue['CB'].get_vector()
        ca_vec = residue['CA'].get_vector()
        
        cb_arr = np.array([cb[0], cb[1], cb[2]])
        ca_arr = np.array([ca_vec[0], ca_vec[1], ca_vec[2]])
        
        direction = cb_arr - ca_arr
        direction = direction / np.linalg.norm(direction)
        
        cg_arr = cb_arr + 1.52 * direction
        SidechainBuilder._add_atom(residue, 'CG', cg_arr, 'C', 10)
        
        sd_arr = cg_arr + 1.81 * direction
        SidechainBuilder._add_atom(residue, 'SD', sd_arr, 'S', 11)
        
        ce_arr = sd_arr + 1.81 * direction
        SidechainBuilder._add_atom(residue, 'CE', ce_arr, 'C', 12)
    
    @staticmethod
    def _build_lys(residue, n, ca, c):
        """Build LYS sidechain: CB-CG-CD-CE-NZ"""
        if not residue.has_id('CB'):
            return
        cb = residue['CB'].get_vector()
        ca_vec = residue['CA'].get_vector()
        
        cb_arr = np.array([cb[0], cb[1], cb[2]])
        ca_arr = np.array([ca_vec[0], ca_vec[1], ca_vec[2]])
        
        direction = cb_arr - ca_arr
        direction = direction / np.linalg.norm(direction)
        
        cg_arr = cb_arr + 1.52 * direction
        SidechainBuilder._add_atom(residue, 'CG', cg_arr, 'C', 10)
        
        cd_arr = cg_arr + 1.52 * direction
        SidechainBuilder._add_atom(residue, 'CD', cd_arr, 'C', 11)
        
        ce_arr = cd_arr + 1.52 * direction
        SidechainBuilder._add_atom(residue, 'CE', ce_arr, 'C', 12)
        
        nz_arr = ce_arr + 1.49 * direction
        SidechainBuilder._add_atom(residue, 'NZ', nz_arr, 'N', 13)
    
    @staticmethod
    def _build_arg(residue, n, ca, c):
        """Build ARG sidechain: CB-CG-CD-NE-CZ-NH1/NH2"""
        if not residue.has_id('CB'):
            return
        cb = residue['CB'].get_vector()
        ca_vec = residue['CA'].get_vector()
        c_vec = residue['C'].get_vector()
        
        cb_arr = np.array([cb[0], cb[1], cb[2]])
        ca_arr = np.array([ca_vec[0], ca_vec[1], ca_vec[2]])
        c_arr = np.array([c_vec[0], c_vec[1], c_vec[2]])
        
        direction = cb_arr - ca_arr
        direction = direction / np.linalg.norm(direction)
        
        cg_arr = cb_arr + 1.52 * direction
        SidechainBuilder._add_atom(residue, 'CG', cg_arr, 'C', 10)
        
        cd_arr = cg_arr + 1.52 * direction
        SidechainBuilder._add_atom(residue, 'CD', cd_arr, 'C', 11)
        
        ne_arr = cd_arr + 1.46 * direction
        SidechainBuilder._add_atom(residue, 'NE', ne_arr, 'N', 12)
        
        cz_arr = ne_arr + 1.33 * direction
        SidechainBuilder._add_atom(residue, 'CZ', cz_arr, 'C', 13)
        
        perp = np.cross(c_arr - ca_arr, cb_arr - ca_arr)
        perp = perp / np.linalg.norm(perp)
        
        nh1_dir = direction + 0.3 * perp
        nh1_arr = cz_arr + 1.33 * nh1_dir / np.linalg.norm(nh1_dir)
        nh2_dir = direction - 0.3 * perp
        nh2_arr = cz_arr + 1.33 * nh2_dir / np.linalg.norm(nh2_dir)
        
        SidechainBuilder._add_atom(residue, 'NH1', nh1_arr, 'N', 14)
        SidechainBuilder._add_atom(residue, 'NH2', nh2_arr, 'N', 15)
    
    @staticmethod
    def _build_leu(residue, n, ca, c):
        """Build LEU sidechain: CB-CG-CD1/CD2"""
        if not residue.has_id('CB'):
            return
        cb = residue['CB'].get_vector()
        ca_vec = residue['CA'].get_vector()
        c_vec = residue['C'].get_vector()
        
        cb_arr = np.array([cb[0], cb[1], cb[2]])
        ca_arr = np.array([ca_vec[0], ca_vec[1], ca_vec[2]])
        c_arr = np.array([c_vec[0], c_vec[1], c_vec[2]])
        
        direction = cb_arr - ca_arr
        direction = direction / np.linalg.norm(direction)
        
        cg_arr = cb_arr + 1.52 * direction
        SidechainBuilder._add_atom(residue, 'CG', cg_arr, 'C', 10)
        
        perp = np.cross(c_arr - ca_arr, cb_arr - ca_arr)
        perp = perp / np.linalg.norm(perp)
        
        cd1_dir = direction + 0.5 * perp
        cd1_arr = cg_arr + 1.52 * cd1_dir / np.linalg.norm(cd1_dir)
        cd2_dir = direction - 0.5 * perp
        cd2_arr = cg_arr + 1.52 * cd2_dir / np.linalg.norm(cd2_dir)
        
        SidechainBuilder._add_atom(residue, 'CD1', cd1_arr, 'C', 11)
        SidechainBuilder._add_atom(residue, 'CD2', cd2_arr, 'C', 12)
    
    @staticmethod
    def _build_ile(residue, n, ca, c):
        """Build ILE sidechain: CB-CG1-CD1, CB-CG2"""
        if not residue.has_id('CB'):
            return
        cb = residue['CB'].get_vector()
        ca_vec = residue['CA'].get_vector()
        c_vec = residue['C'].get_vector()
        
        cb_arr = np.array([cb[0], cb[1], cb[2]])
        ca_arr = np.array([ca_vec[0], ca_vec[1], ca_vec[2]])
        c_arr = np.array([c_vec[0], c_vec[1], c_vec[2]])
        
        direction = cb_arr - ca_arr
        direction = direction / np.linalg.norm(direction)
        
        perp = np.cross(c_arr - ca_arr, cb_arr - ca_arr)
        perp = perp / np.linalg.norm(perp)
        
        # CG1
        cg1_arr = cb_arr + 1.52 * direction
        SidechainBuilder._add_atom(residue, 'CG1', cg1_arr, 'C', 10)
        
        # CG2 - branched
        cg2_arr = cb_arr + 1.52 * perp
        SidechainBuilder._add_atom(residue, 'CG2', cg2_arr, 'C', 11)
        
        # CD1
        cd1_arr = cg1_arr + 1.52 * direction
        SidechainBuilder._add_atom(residue, 'CD1', cd1_arr, 'C', 12)
    
    @staticmethod
    def _build_phe(residue, n, ca, c):
        """Build PHE sidechain: CB-CG-CD1/CD2-CE1/CE2-CZ (aromatic ring)"""
        if not residue.has_id('CB'):
            return
        cb = residue['CB'].get_vector()
        ca_vec = residue['CA'].get_vector()
        c_vec = residue['C'].get_vector()
        
        cb_arr = np.array([cb[0], cb[1], cb[2]])
        ca_arr = np.array([ca_vec[0], ca_vec[1], ca_vec[2]])
        c_arr = np.array([c_vec[0], c_vec[1], c_vec[2]])
        
        direction = cb_arr - ca_arr
        direction = direction / np.linalg.norm(direction)
        
        cg_arr = cb_arr + 1.50 * direction
        SidechainBuilder._add_atom(residue, 'CG', cg_arr, 'C', 10)
        
        # Aromatic ring (simplified planar geometry)
        perp = np.cross(c_arr - ca_arr, cb_arr - ca_arr)
        perp = perp / np.linalg.norm(perp)
        
        cd1_dir = direction + 0.5 * perp
        cd1_arr = cg_arr + 1.39 * cd1_dir / np.linalg.norm(cd1_dir)
        cd2_dir = direction - 0.5 * perp
        cd2_arr = cg_arr + 1.39 * cd2_dir / np.linalg.norm(cd2_dir)
        SidechainBuilder._add_atom(residue, 'CD1', cd1_arr, 'C', 11)
        SidechainBuilder._add_atom(residue, 'CD2', cd2_arr, 'C', 12)
        
        ce1_arr = cd1_arr + 1.39 * direction
        ce2_arr = cd2_arr + 1.39 * direction
        SidechainBuilder._add_atom(residue, 'CE1', ce1_arr, 'C', 13)
        SidechainBuilder._add_atom(residue, 'CE2', ce2_arr, 'C', 14)
        
        cz_arr = cg_arr + 2.78 * direction
        SidechainBuilder._add_atom(residue, 'CZ', cz_arr, 'C', 15)
    
    @staticmethod
    def _build_tyr(residue, n, ca, c):
        """Build TYR sidechain: like PHE + OH"""
        SidechainBuilder._build_phe(residue, n, ca, c)
        
        if residue.has_id('CZ'):
            cz = residue['CZ'].get_vector()
            cb = residue['CB'].get_vector()
            
            cz_arr = np.array([cz[0], cz[1], cz[2]])
            cb_arr = np.array([cb[0], cb[1], cb[2]])
            
            direction = cz_arr - cb_arr
            direction = direction / np.linalg.norm(direction)
            oh_arr = cz_arr + 1.36 * direction
            SidechainBuilder._add_atom(residue, 'OH', oh_arr, 'O', 16)
    
    @staticmethod
    def _build_trp(residue, n, ca, c):
        """Build TRP sidechain (indole ring - simplified)"""
        if not residue.has_id('CB'):
            return
        cb = residue['CB'].get_vector()
        ca_vec = residue['CA'].get_vector()
        c_vec = residue['C'].get_vector()
        
        cb_arr = np.array([cb[0], cb[1], cb[2]])
        ca_arr = np.array([ca_vec[0], ca_vec[1], ca_vec[2]])
        c_arr = np.array([c_vec[0], c_vec[1], c_vec[2]])
        
        direction = cb_arr - ca_arr
        direction = direction / np.linalg.norm(direction)
        
        perp = np.cross(c_arr - ca_arr, cb_arr - ca_arr)
        perp = perp / np.linalg.norm(perp)
        
        cg_arr = cb_arr + 1.50 * direction
        SidechainBuilder._add_atom(residue, 'CG', cg_arr, 'C', 10)
        
        cd1_dir = direction + 0.3 * perp
        cd1_arr = cg_arr + 1.37 * cd1_dir / np.linalg.norm(cd1_dir)
        SidechainBuilder._add_atom(residue, 'CD1', cd1_arr, 'C', 11)
        
        cd2_dir = direction - 0.3 * perp
        cd2_arr = cg_arr + 1.43 * cd2_dir / np.linalg.norm(cd2_dir)
        SidechainBuilder._add_atom(residue, 'CD2', cd2_arr, 'C', 12)
        
        ne1_arr = cd1_arr + 1.38 * direction
        SidechainBuilder._add_atom(residue, 'NE1', ne1_arr, 'N', 13)
        
        ce2_arr = cd2_arr + 1.40 * direction
        SidechainBuilder._add_atom(residue, 'CE2', ce2_arr, 'C', 14)
        
        ce3_arr = cd2_arr + 1.40 * (-perp)
        SidechainBuilder._add_atom(residue, 'CE3', ce3_arr, 'C', 15)
        
        cz2_arr = ce2_arr + 1.40 * direction
        SidechainBuilder._add_atom(residue, 'CZ2', cz2_arr, 'C', 16)
        
        cz3_arr = ce3_arr + 1.40 * direction
        SidechainBuilder._add_atom(residue, 'CZ3', cz3_arr, 'C', 17)
        
        ch2_arr = cz2_arr + 1.40 * (-perp)
        SidechainBuilder._add_atom(residue, 'CH2', ch2_arr, 'C', 18)
    
    @staticmethod
    def _build_his(residue, n, ca, c):
        """Build HIS sidechain (imidazole ring)"""
        if not residue.has_id('CB'):
            return
        cb = residue['CB'].get_vector()
        ca_vec = residue['CA'].get_vector()
        c_vec = residue['C'].get_vector()
        
        cb_arr = np.array([cb[0], cb[1], cb[2]])
        ca_arr = np.array([ca_vec[0], ca_vec[1], ca_vec[2]])
        c_arr = np.array([c_vec[0], c_vec[1], c_vec[2]])
        
        direction = cb_arr - ca_arr
        direction = direction / np.linalg.norm(direction)
        
        perp = np.cross(c_arr - ca_arr, cb_arr - ca_arr)
        perp = perp / np.linalg.norm(perp)
        
        cg_arr = cb_arr + 1.50 * direction
        SidechainBuilder._add_atom(residue, 'CG', cg_arr, 'C', 10)
        
        nd1_dir = direction + 0.4 * perp
        nd1_arr = cg_arr + 1.38 * nd1_dir / np.linalg.norm(nd1_dir)
        SidechainBuilder._add_atom(residue, 'ND1', nd1_arr, 'N', 11)
        
        cd2_dir = direction - 0.4 * perp
        cd2_arr = cg_arr + 1.37 * cd2_dir / np.linalg.norm(cd2_dir)
        SidechainBuilder._add_atom(residue, 'CD2', cd2_arr, 'C', 12)
        
        ce1_arr = nd1_arr + 1.32 * direction
        SidechainBuilder._add_atom(residue, 'CE1', ce1_arr, 'C', 13)
        
        ne2_arr = cd2_arr + 1.38 * direction
        SidechainBuilder._add_atom(residue, 'NE2', ne2_arr, 'N', 14)
    
    @staticmethod
    def _build_pro(residue, n, ca, c):
        """Build PRO sidechain: CB-CG-CD (ring back to N)"""
        if not residue.has_id('CB'):
            return
        cb = residue['CB'].get_vector()
        ca_vec = residue['CA'].get_vector()
        
        cb_arr = np.array([cb[0], cb[1], cb[2]])
        ca_arr = np.array([ca_vec[0], ca_vec[1], ca_vec[2]])
        
        direction = cb_arr - ca_arr
        direction = direction / np.linalg.norm(direction)
        
        cg_arr = cb_arr + 1.50 * direction
        SidechainBuilder._add_atom(residue, 'CG', cg_arr, 'C', 10)
        
        # CD connects back to N
        n_vec = residue['N'].get_vector()
        n_arr = np.array([n_vec[0], n_vec[1], n_vec[2]])
        cd_direction = n_arr - cg_arr
        cd_direction = cd_direction / np.linalg.norm(cd_direction)
        cd_arr = cg_arr + 1.51 * cd_direction
        SidechainBuilder._add_atom(residue, 'CD', cd_arr, 'C', 11)
