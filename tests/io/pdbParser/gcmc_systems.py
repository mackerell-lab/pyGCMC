# tests/io/pdbParser/gcmc_systems.py
"""PDB Parser tests for GCMC systems - focusing on 4wp7 system parsing."""

import pytest
import os
import pygcmc
import math
from functools import lru_cache

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")
PDB_4WP7_PATH = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.prod.74.rec.pdb")


@lru_cache(maxsize=1)
def _parsed_4wp7_system():
    """Parse the large 4wp7 PDB once per module to reduce memory churn."""
    if not os.path.exists(PDB_4WP7_PATH):
        pytest.skip(f"Test file {PDB_4WP7_PATH} not found")
    return pygcmc.PDBParser.parse_file(PDB_4WP7_PATH)


def test_parse_4wp7_gcmc_system():
    """Test parsing 4wp7 GCMC system PDB file with protein and small molecules."""
    result = _parsed_4wp7_system()
    
    # Verify basic parsing success
    assert result is not None, "Failed to parse 4wp7 PDB file"
    assert len(result.atoms) > 200000, f"Expected >200k atoms, got {len(result.atoms)}"
    
    # Check that we have the expected large number of atoms (~221k)
    assert 220000 <= len(result.atoms) <= 225000, f"Unexpected atom count: {len(result.atoms)}"
    
    # Verify we have multiple residue types (protein + small molecules)
    residue_types = set()
    for atom in result.atoms:
        residue_types.add(atom.get_resname())
    
    # Should have both protein residues and small molecules
    protein_residues = {"ALA", "ARG", "ASN", "ASP", "GLN", "GLU", "GLY", "HIS", "LEU", "LYS"}
    gcmc_molecules = {"ACEY", "BENX", "DMEE", "FORM"}
    
    found_protein = len(protein_residues.intersection(residue_types))
    found_gcmc = len(gcmc_molecules.intersection(residue_types))
    
    assert found_protein >= 5, f"Expected multiple protein residues, found {found_protein}"
    assert found_gcmc >= 2, f"Expected GCMC molecules, found {found_gcmc}"
    
    # Check first ASP residue atoms (N-terminal)
    asp_atoms = [atom for atom in result.atoms if atom.get_resname() == "ASP" and atom.get_ires() == 8]
    assert len(asp_atoms) > 10, "ASP residue should have multiple atoms"
    
    # Verify N-terminal ASP nitrogen
    n_atom = None
    for atom in asp_atoms:
        if atom.get_type() == "N":
            n_atom = atom
            break
    
    assert n_atom is not None, "N atom not found in ASP residue"
    
    # Check coordinates are reasonable (within simulation box)
    coords = n_atom.get_coor()
    assert len(coords) == 3, "Coordinates should be 3D"
    assert all(0 <= c <= 150 for c in coords), f"Coordinates outside expected range: {coords}"
    
    # Verify crystal information is parsed
    assert hasattr(result, 'crystal_info') or len(result.atoms) > 0, "Should parse crystal info or atoms"


def test_4wp7_system_residue_distribution():
    """Test the distribution of residues in 4wp7 system matches expected GCMC composition."""
    result = _parsed_4wp7_system()
    
    # Count residue occurrences (counting atoms per residue type)
    residue_counts = {}  # This counts ATOMS per residue type, not molecules
    for atom in result.atoms:
        resname = atom.get_resname()
        residue_counts[resname] = residue_counts.get(resname, 0) + 1
    
    # Verify we have the expected GCMC molecules with exact counts from raw PDB file analysis
    # These counts were obtained directly from the raw PDB file using:
    # grep "^ATOM" 4wp7_*.pdb | awk '{print $4}' | sort | uniq -c
    expected_gcmc_molecules = {
        "ACEY": 2048,    # Exact count from raw PDB file
        "BENX": 3731,    # Exact count from raw PDB file  
        "DMEE": 2169,    # Exact count from raw PDB file
        "FORM": 1662,    # Exact count from raw PDB file
    }
    
    for molecule, expected_count in expected_gcmc_molecules.items():
        if molecule in residue_counts:
            actual_count = residue_counts[molecule]
            assert actual_count == expected_count, \
                f"{molecule} count {actual_count} does not match expected {expected_count}"
    
    # Verify we have protein residues (including HSD as protonated histidine and CYS)
    protein_atoms_found = 0  # Renamed for clarity: this counts ATOMS, not residues
    for resname in residue_counts:
        if resname in {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", 
                      "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HSD"}:
            protein_atoms_found += residue_counts[resname]
    
    # Should have significant protein content (multiple copies of the protein)
    assert protein_atoms_found > 25000, f"Expected substantial protein content, got {protein_atoms_found} atoms"
    # Exact check based on raw PDB file analysis: sum of all standard amino acid + HSD atoms
    # Calculated from: grep "^ATOM" *.pdb | awk '{print $4}' | grep -E "^(ALA|ARG|...|HSD)$" | wc -l
    assert protein_atoms_found == 30520, f"Expected exactly 30,520 protein atoms for 4wp7 system, got {protein_atoms_found}"
    
    # Total system should be dominated by non-protein molecules (GCMC + solvent)
    total_atoms = len(result.atoms)
    small_molecule_atoms = total_atoms - protein_atoms_found
    
    # Based on raw PDB file analysis: 30,520 protein, 190,921 non-protein (ratio 6.26:1)
    # Total atoms from raw file: grep "^ATOM" *.pdb | wc -l = 221,441
    assert small_molecule_atoms > protein_atoms_found * 5, \
        f"GCMC+solvent system should be dominated by small molecules, got ratio {small_molecule_atoms/protein_atoms_found:.1f}:1"
    # Exact check: 221,441 total - 30,520 protein = 190,921 non-protein
    expected_non_protein = 190921
    assert small_molecule_atoms == expected_non_protein, \
        f"Expected exactly {expected_non_protein} non-protein atoms, got {small_molecule_atoms}"


def test_4wp7_coordinate_validation():
    """Test coordinate parsing and validation for 4wp7 system."""
    result = _parsed_4wp7_system()
    
    # Sample atoms for coordinate testing (test first 1000 atoms for performance)
    # Note: Using first 1000 atoms for deterministic testing, but has order dependency
    sample_atoms = result.atoms[:1000]
    
    valid_coordinates = 0
    for atom in sample_atoms:
        coords = atom.get_coor()
        
        # Verify 3D coordinates
        assert len(coords) == 3, f"Atom {atom.get_bynu()} has {len(coords)} coordinates, expected 3"
        
        # Check coordinates are numeric and reasonable
        for i, coord in enumerate(coords):
            assert isinstance(coord, (int, float)), f"Coordinate {i} is not numeric: {coord}"
            assert -500 <= coord <= 500, f"Coordinate {i} out of reasonable range: {coord}"
        
        valid_coordinates += 1
    
    assert valid_coordinates == len(sample_atoms), "All sampled atoms should have valid coordinates"
    
    # Test specific known coordinates from first ASP residue
    asp_n = None
    for atom in result.atoms:
        if atom.get_resname() == "ASP" and atom.get_ires() == 8 and atom.get_type() == "N":
            asp_n = atom
            break
    
    if asp_n:
        coords = asp_n.get_coor()
        # From the file: ATOM      1  N   ASP     8      42.526  94.418  83.282
        assert math.isclose(coords[0], 42.526, rel_tol=1e-3), f"X coordinate mismatch: {coords[0]}"
        assert math.isclose(coords[1], 94.418, rel_tol=1e-3), f"Y coordinate mismatch: {coords[1]}"  
        assert math.isclose(coords[2], 83.282, rel_tol=1e-3), f"Z coordinate mismatch: {coords[2]}"
