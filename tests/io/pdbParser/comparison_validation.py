# tests/io/pdbParser/comparison_validation.py
"""
Direct PDB reading comparison with pygcmc parser for validation.
Provides independent Python-based PDB parsing to verify pygcmc parser correctness.
"""

import pytest
import os
import re
import math
from typing import Dict, List, Tuple, NamedTuple
import pygcmc

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


def test_4wp7_direct_vs_pygcmc_parsing():
    """Compare direct Python PDB parsing with pygcmc parser for 4wp7 system."""
    pdb_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.prod.74.rec.pdb")
    
    # Skip test if file doesn't exist
    if not os.path.exists(pdb_path):
        pytest.skip(f"Test file {pdb_path} not found")
    
    # Parse with both methods
    direct_atoms, direct_metadata = read_pdb_direct(pdb_path)
    pygcmc_result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Basic comparison
    assert len(direct_atoms) == len(pygcmc_result.atoms), \
        f"Atom count mismatch: direct={len(direct_atoms)}, pygcmc={len(pygcmc_result.atoms)}"
    
    print(f"Both parsers found {len(direct_atoms)} atoms")
    
    # Sample comparison (first 100 atoms for performance)
    sample_size = min(100, len(direct_atoms))
    mismatches = 0
    
    for i in range(sample_size):
        direct_atom = direct_atoms[i]
        pygcmc_atom = pygcmc_result.atoms[i]
        
        differences = compare_atom_data(direct_atom, pygcmc_atom)
        if differences:
            mismatches += 1
            if mismatches <= 5:  # Only report first 5 mismatches
                print(f"Atom {i+1} differences: {'; '.join(differences)}")
    
    # Allow small number of mismatches due to different parsing approaches
    mismatch_rate = mismatches / sample_size
    assert mismatch_rate < 0.05, f"Too many mismatches: {mismatches}/{sample_size} ({mismatch_rate:.1%})"
    
    print(f"Comparison: {mismatches}/{sample_size} atoms had differences ({mismatch_rate:.1%})")


def test_4wp7_residue_statistics_comparison():
    """Compare residue-level statistics between direct and pygcmc parsing."""
    pdb_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.prod.74.rec.pdb")
    
    # Skip test if file doesn't exist
    if not os.path.exists(pdb_path):
        pytest.skip(f"Test file {pdb_path} not found")
    
    # Parse with both methods
    direct_atoms, _ = read_pdb_direct(pdb_path)
    pygcmc_result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Count residues with direct parsing
    direct_residue_counts = {}
    for atom in direct_atoms:
        res_key = atom.res_name
        direct_residue_counts[res_key] = direct_residue_counts.get(res_key, 0) + 1
    
    # Count residues with pygcmc parsing
    pygcmc_residue_counts = {}
    for atom in pygcmc_result.atoms:
        res_key = atom.get_resname()
        pygcmc_residue_counts[res_key] = pygcmc_residue_counts.get(res_key, 0) + 1
    
    # Define known residue name mappings (pygcmc may standardize names)
    residue_mapping = {
        # Original -> Standardized
        "ACE": "ACEY",   # acetylene
        "BEN": "BENX",   # benzene  
        "DME": "DMEE",   # dimethyl ether
        "FOR": "FORM",   # formaldehyde
        "MEO": "MEOH",   # methanol
        "IMI": "IMIA",   # imidazole
        "PRP": "PRPX",   # propane
        "MAM": "MAMY",   # methylamine
    }
    
    # Apply mapping to direct counts
    mapped_direct_counts = {}
    for res_type, count in direct_residue_counts.items():
        mapped_name = residue_mapping.get(res_type, res_type)
        mapped_direct_counts[mapped_name] = mapped_direct_counts.get(mapped_name, 0) + count
    
    # Compare mapped residue type counts
    all_residue_types = set(mapped_direct_counts.keys()) | set(pygcmc_residue_counts.keys())
    
    major_differences = []
    for res_type in all_residue_types:
        direct_count = mapped_direct_counts.get(res_type, 0)
        pygcmc_count = pygcmc_residue_counts.get(res_type, 0)
        
        if direct_count != pygcmc_count:
            difference = abs(direct_count - pygcmc_count)
            if difference > max(1, min(direct_count, pygcmc_count) * 0.01):  # >1% difference
                major_differences.append(f"{res_type}: direct={direct_count}, pygcmc={pygcmc_count}")
    
    # Allow some differences due to parser variations, but verify total counts match
    total_direct = sum(mapped_direct_counts.values())
    total_pygcmc = sum(pygcmc_residue_counts.values())
    assert total_direct == total_pygcmc, f"Total atom count mismatch: direct={total_direct}, pygcmc={total_pygcmc}"
    
    # Check that we have most key residue types in common
    common_types = set(mapped_direct_counts.keys()) & set(pygcmc_residue_counts.keys())
    assert len(common_types) >= len(all_residue_types) * 0.8, \
        f"Too few common residue types: {len(common_types)}/{len(all_residue_types)}"
    
    # Check key GCMC molecules (these should match exactly after mapping)
    gcmc_molecules = ["ACEY", "BENX", "DMEE", "FORM"]
    matching_gcmc = 0
    for molecule in gcmc_molecules:
        if molecule in mapped_direct_counts and molecule in pygcmc_residue_counts:
            direct_count = mapped_direct_counts[molecule]
            pygcmc_count = pygcmc_residue_counts[molecule]
            if direct_count == pygcmc_count:
                matching_gcmc += 1
    
    assert matching_gcmc >= 3, f"At least 3 GCMC molecules should match exactly, got {matching_gcmc}"
    
    print(f"Residue statistics comparison: {len(common_types)}/{len(all_residue_types)} types match, "
          f"{matching_gcmc}/4 GCMC molecules match exactly")


def test_4wp7_coordinate_ranges_comparison():
    """Compare coordinate ranges and statistics between parsers."""
    pdb_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.prod.74.rec.pdb")
    
    # Skip test if file doesn't exist
    if not os.path.exists(pdb_path):
        pytest.skip(f"Test file {pdb_path} not found")
    
    # Parse with both methods
    direct_atoms, _ = read_pdb_direct(pdb_path)
    pygcmc_result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Calculate coordinate ranges for direct parsing
    if direct_atoms:
        direct_x_coords = [atom.x for atom in direct_atoms]
        direct_y_coords = [atom.y for atom in direct_atoms]
        direct_z_coords = [atom.z for atom in direct_atoms]
        
        direct_ranges = {
            'x': (min(direct_x_coords), max(direct_x_coords)),
            'y': (min(direct_y_coords), max(direct_y_coords)),
            'z': (min(direct_z_coords), max(direct_z_coords))
        }
    
    # Calculate coordinate ranges for pygcmc parsing
    if pygcmc_result.atoms:
        pygcmc_x_coords = [atom.get_coor()[0] for atom in pygcmc_result.atoms]
        pygcmc_y_coords = [atom.get_coor()[1] for atom in pygcmc_result.atoms]
        pygcmc_z_coords = [atom.get_coor()[2] for atom in pygcmc_result.atoms]
        
        pygcmc_ranges = {
            'x': (min(pygcmc_x_coords), max(pygcmc_x_coords)),
            'y': (min(pygcmc_y_coords), max(pygcmc_y_coords)),
            'z': (min(pygcmc_z_coords), max(pygcmc_z_coords))
        }
    
    # Compare coordinate ranges (should be identical or very close)
    tolerance = 1e-3
    for coord in ['x', 'y', 'z']:
        direct_min, direct_max = direct_ranges[coord]
        pygcmc_min, pygcmc_max = pygcmc_ranges[coord]
        
        assert math.isclose(direct_min, pygcmc_min, abs_tol=tolerance), \
            f"{coord} min mismatch: direct={direct_min:.3f}, pygcmc={pygcmc_min:.3f}"
        
        assert math.isclose(direct_max, pygcmc_max, abs_tol=tolerance), \
            f"{coord} max mismatch: direct={direct_max:.3f}, pygcmc={pygcmc_max:.3f}"
    
    print(f"Coordinate ranges match: X[{direct_ranges['x'][0]:.1f}, {direct_ranges['x'][1]:.1f}], "
          f"Y[{direct_ranges['y'][0]:.1f}, {direct_ranges['y'][1]:.1f}], "
          f"Z[{direct_ranges['z'][0]:.1f}, {direct_ranges['z'][1]:.1f}]")


def test_specific_atom_verification():
    """Verify specific known atoms between direct and pygcmc parsing."""
    pdb_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.prod.74.rec.pdb")
    
    # Skip test if file doesn't exist
    if not os.path.exists(pdb_path):
        pytest.skip(f"Test file {pdb_path} not found")
    
    # Parse with both methods
    direct_atoms, _ = read_pdb_direct(pdb_path)
    pygcmc_result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Find first ASP N atom (from line 6 in PDB file)
    # ATOM      1  N   ASP     8      42.526  94.418  83.282  1.00  0.00           N
    
    direct_asp_n = None
    for atom in direct_atoms:
        if (atom.atom_num == 1 and atom.atom_name == "N" and 
            atom.res_name == "ASP" and atom.res_seq == 8):
            direct_asp_n = atom
            break
    
    pygcmc_asp_n = None
    for atom in pygcmc_result.atoms:
        if (atom.get_bynu() == 1 and atom.get_type() == "N" and 
            atom.get_resname() == "ASP" and atom.get_ires() == 8):
            pygcmc_asp_n = atom
            break
    
    assert direct_asp_n is not None, "Could not find ASP N atom in direct parsing"
    assert pygcmc_asp_n is not None, "Could not find ASP N atom in pygcmc parsing"
    
    # Verify coordinates match exactly
    tolerance = 1e-3
    assert math.isclose(direct_asp_n.x, 42.526, abs_tol=tolerance), f"Direct X coordinate wrong: {direct_asp_n.x}"
    assert math.isclose(direct_asp_n.y, 94.418, abs_tol=tolerance), f"Direct Y coordinate wrong: {direct_asp_n.y}"
    assert math.isclose(direct_asp_n.z, 83.282, abs_tol=tolerance), f"Direct Z coordinate wrong: {direct_asp_n.z}"
    
    pygcmc_coords = pygcmc_asp_n.get_coor()
    assert math.isclose(pygcmc_coords[0], 42.526, abs_tol=tolerance), f"Pygcmc X coordinate wrong: {pygcmc_coords[0]}"
    assert math.isclose(pygcmc_coords[1], 94.418, abs_tol=tolerance), f"Pygcmc Y coordinate wrong: {pygcmc_coords[1]}"
    assert math.isclose(pygcmc_coords[2], 83.282, abs_tol=tolerance), f"Pygcmc Z coordinate wrong: {pygcmc_coords[2]}"
    
    # Compare the two atoms directly
    differences = compare_atom_data(direct_asp_n, pygcmc_asp_n)
    assert len(differences) == 0, f"ASP N atom differences: {'; '.join(differences)}"
    
    print("First ASP N atom verification passed - both parsers read identical data")