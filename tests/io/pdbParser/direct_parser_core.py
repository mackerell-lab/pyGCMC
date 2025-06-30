# tests/io/pdbParser/direct_parser_core.py
"""
Core direct PDB parsing functionality.
Provides independent Python-based PDB parsing for comparison testing.
"""

import os
import math
from typing import Dict, List, Tuple, NamedTuple

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


class DirectPDBAtom(NamedTuple):
    """Simple atom structure from direct PDB parsing."""
    record_type: str    # ATOM or HETATM
    atom_num: int       # Atom serial number
    atom_name: str      # Atom name
    alt_loc: str        # Alternate location indicator
    res_name: str       # Residue name
    chain_id: str       # Chain identifier
    res_seq: int        # Residue sequence number
    icode: str          # Insertion code
    x: float            # X coordinate
    y: float            # Y coordinate
    z: float            # Z coordinate
    occupancy: float    # Occupancy
    temp_factor: float  # Temperature factor
    element: str        # Element symbol
    charge: str         # Charge


def parse_pdb_line_direct(line: str) -> DirectPDBAtom:
    """
    Direct PDB line parsing following PDB format specification.
    
    PDB format (fixed-width):
    COLUMNS        DATA TYPE       CONTENTS
    1-6            Record name     "ATOM  " or "HETATM"
    7-11           Integer         Atom serial number
    13-16          Atom            Atom name
    17             Character       Alternate location indicator
    18-20          Residue name    Residue name
    22             Character       Chain identifier
    23-26          Integer         Residue sequence number
    27             AChar           Code for insertion of residues
    31-38          Real(8.3)       Orthogonal coordinates for X
    39-46          Real(8.3)       Orthogonal coordinates for Y
    47-54          Real(8.3)       Orthogonal coordinates for Z
    55-60          Real(6.2)       Occupancy
    61-66          Real(6.2)       Temperature factor
    77-78          LString(2)      Element symbol
    79-80          LString(2)      Charge on the atom
    """
    if len(line) < 54:
        raise ValueError(f"PDB line too short: {len(line)} chars")
    
    try:
        record_type = line[0:6].strip()
        atom_num = int(line[6:11])
        atom_name = line[12:16].strip()
        alt_loc = line[16:17].strip()
        res_name = line[17:20].strip()
        chain_id = line[21:22].strip()
        res_seq = int(line[22:26])
        icode = line[26:27].strip()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        
        # Optional fields
        occupancy = float(line[54:60]) if len(line) > 54 and line[54:60].strip() else 1.0
        temp_factor = float(line[60:66]) if len(line) > 60 and line[60:66].strip() else 0.0
        element = line[76:78].strip() if len(line) > 76 else ""
        charge = line[78:80].strip() if len(line) > 78 else ""
        
        return DirectPDBAtom(
            record_type, atom_num, atom_name, alt_loc, res_name, chain_id,
            res_seq, icode, x, y, z, occupancy, temp_factor, element, charge
        )
    except (ValueError, IndexError) as e:
        raise ValueError(f"Failed to parse PDB line: {line.strip()[:50]}... Error: {e}")


def read_pdb_direct(pdb_path: str) -> Tuple[List[DirectPDBAtom], Dict[str, str]]:
    """
    Read PDB file directly using Python, returning atoms and metadata.
    
    Returns:
        Tuple of (atom_list, metadata_dict)
    """
    atoms = []
    metadata = {}
    
    with open(pdb_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.rstrip('\n\r')
            
            if line.startswith(('ATOM  ', 'HETATM')):
                try:
                    atom = parse_pdb_line_direct(line)
                    atoms.append(atom)
                except ValueError as e:
                    print(f"Warning: Skipping invalid PDB line {line_num}: {e}")
                    continue
            elif line.startswith('CRYST1'):
                metadata['crystal'] = line
            elif line.startswith('TITLE'):
                metadata['title'] = line[6:].strip()
            elif line.startswith('REMARK'):
                if 'remarks' not in metadata:
                    metadata['remarks'] = []
                metadata['remarks'].append(line)
    
    return atoms, metadata


def compare_atom_data(direct_atom: DirectPDBAtom, pygcmc_atom, tolerance: float = 1e-3) -> List[str]:
    """
    Compare direct PDB atom with pygcmc parsed atom.
    
    Returns:
        List of differences found (empty if atoms match)
    """
    differences = []
    
    # Basic identification
    if direct_atom.atom_num != pygcmc_atom.get_bynu():
        differences.append(f"Atom number: direct={direct_atom.atom_num}, pygcmc={pygcmc_atom.get_bynu()}")
    
    if direct_atom.atom_name != pygcmc_atom.get_type():
        differences.append(f"Atom name: direct='{direct_atom.atom_name}', pygcmc='{pygcmc_atom.get_type()}'")
    
    if direct_atom.res_name != pygcmc_atom.get_resname():
        differences.append(f"Residue name: direct='{direct_atom.res_name}', pygcmc='{pygcmc_atom.get_resname()}'")
    
    if direct_atom.res_seq != pygcmc_atom.get_ires():
        differences.append(f"Residue number: direct={direct_atom.res_seq}, pygcmc={pygcmc_atom.get_ires()}")
    
    # Coordinates
    pygcmc_coords = pygcmc_atom.get_coor()
    if not math.isclose(direct_atom.x, pygcmc_coords[0], abs_tol=tolerance):
        differences.append(f"X coordinate: direct={direct_atom.x:.3f}, pygcmc={pygcmc_coords[0]:.3f}")
    
    if not math.isclose(direct_atom.y, pygcmc_coords[1], abs_tol=tolerance):
        differences.append(f"Y coordinate: direct={direct_atom.y:.3f}, pygcmc={pygcmc_coords[1]:.3f}")
    
    if not math.isclose(direct_atom.z, pygcmc_coords[2], abs_tol=tolerance):
        differences.append(f"Z coordinate: direct={direct_atom.z:.3f}, pygcmc={pygcmc_coords[2]:.3f}")
    
    # PDB specific fields
    if hasattr(pygcmc_atom, 'get_occupancy'):
        if not math.isclose(direct_atom.occupancy, pygcmc_atom.get_occupancy(), abs_tol=tolerance):
            differences.append(f"Occupancy: direct={direct_atom.occupancy:.2f}, pygcmc={pygcmc_atom.get_occupancy():.2f}")
    
    if hasattr(pygcmc_atom, 'get_tempfactor'):
        if not math.isclose(direct_atom.temp_factor, pygcmc_atom.get_tempfactor(), abs_tol=tolerance):
            differences.append(f"Temperature factor: direct={direct_atom.temp_factor:.2f}, pygcmc={pygcmc_atom.get_tempfactor():.2f}")
    
    return differences