# tests/io/pdbParser/enhanced_4wp7_tests.py
"""
Enhanced PDB Parser tests for 4wp7 system - Deep analysis of GCMC PDB parsing.

This file contains comprehensive tests for the 4wp7 PDB file that go beyond
basic parsing verification to validate molecular integrity, crystal structure,
and chemical composition.

Tests include:
- Molecular integrity and structure validation
- Crystal cell and spatial distribution analysis
- Protein structure deep verification
- Chemical composition and atom type validation
"""

import pytest
import os
import time
import math
from pygcmc.io import PDBParser

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_4wp7_molecule_integrity():
    """Test molecular integrity and structure validation for GCMC molecules and water."""
    pdb_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.prod.74.rec.pdb")
    
    # Skip test if file doesn't exist
    if not os.path.exists(pdb_path):
        pytest.skip(f"Test file {pdb_path} not found")
    
    result = PDBParser.parse_file(pdb_path)
    assert result is not None, "Failed to parse 4wp7 PDB file"
    
    # Expected molecular structures based on actual PDB parser output
    expected_gcmc_molecules = {
        "BENX": {"atoms_per_molecule": 13, "description": "Benzene with hydrogens"},
        "PRPX": {"atoms_per_molecule": 12, "description": "Propane with explicit hydrogens"},
        "IMIA": {"atoms_per_molecule": 9, "description": "Imidazole"},
        "MAMY": {"atoms_per_molecule": 9, "description": "Methylammonium"},
        "DMEE": {"atoms_per_molecule": 9, "description": "Dimethyl ether"},
        "ACEY": {"atoms_per_molecule": 8, "description": "Acetate"},
        "MEOH": {"atoms_per_molecule": 6, "description": "Methanol"},
        "FORM": {"atoms_per_molecule": 6, "description": "Formamide"}
    }
    
    # Count atoms per residue - now using actual residue objects from parser
    residue_atom_counts = {}
    residue_instances = {}
    
    # Use the residue objects from the parser instead of counting atoms manually
    for residue in result.residues:
        res_name = residue.get_resname()
        if res_name in expected_gcmc_molecules or res_name == "SOL":
            atom_count = residue.atom_count()
            res_id = f"{res_name}_{residue.get_ires()}_{id(residue)}"  # Include object ID to make unique
            
            residue_atom_counts[res_id] = atom_count
            if res_name not in residue_instances:
                residue_instances[res_name] = []
            residue_instances[res_name].append(res_id)
    
    # Verify GCMC molecule structures
    for mol_name, expected_info in expected_gcmc_molecules.items():
        if mol_name in residue_instances:
            # Check a sample of molecules of this type
            sample_molecules = residue_instances[mol_name][:5]  # Check first 5 instances
            
            for res_id in sample_molecules:
                actual_atoms = residue_atom_counts[res_id]
                expected_atoms = expected_info["atoms_per_molecule"]
                
                assert actual_atoms == expected_atoms, \
                    f"{mol_name} molecule {res_id} ({expected_info['description']}): " \
                    f"expected {expected_atoms} atoms, got {actual_atoms}"
    
    # Verify TIP3P water structure (SOL residues should have exactly 3 atoms)
    if "SOL" in residue_instances:
        water_sample = residue_instances["SOL"][:10]  # Check first 10 water molecules
        
        for res_id in water_sample:
            actual_atoms = residue_atom_counts[res_id]
            assert actual_atoms == 3, \
                f"Water molecule {res_id} should have 3 atoms (OW + 2×HW), got {actual_atoms}"
    
    print(f"✓ Molecular integrity validated: GCMC molecules and water structures correct")


def test_4wp7_crystal_cell_parsing():
    """Test crystal cell information parsing and coordinate validation."""
    pdb_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.prod.74.rec.pdb")
    
    # Skip test if file doesn't exist
    if not os.path.exists(pdb_path):
        pytest.skip(f"Test file {pdb_path} not found")
    
    result = PDBParser.parse_file(pdb_path)
    assert result is not None, "Failed to parse 4wp7 PDB file"
    
    # Read raw file to extract CRYST1 record
    cryst1_info = None
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("CRYST1"):
                # Parse CRYST1 record: a, b, c, alpha, beta, gamma
                parts = line.split()
                cryst1_info = {
                    'a': float(parts[1]),
                    'b': float(parts[2]), 
                    'c': float(parts[3]),
                    'alpha': float(parts[4]),
                    'beta': float(parts[5]),
                    'gamma': float(parts[6])
                }
                break
    
    assert cryst1_info is not None, "CRYST1 record not found in PDB file"
    
    # Verify expected crystal cell dimensions (from file header)
    expected_dims = {'a': 132.126, 'b': 132.213, 'c': 122.070}
    for dim, expected_val in expected_dims.items():
        actual_val = cryst1_info[dim]
        assert abs(actual_val - expected_val) < 0.01, \
            f"Crystal cell {dim}: expected {expected_val}, got {actual_val}"
    
    # Verify orthogonal cell (90 degree angles)
    for angle in ['alpha', 'beta', 'gamma']:
        assert abs(cryst1_info[angle] - 90.0) < 0.01, \
            f"Expected orthogonal cell, {angle} = {cryst1_info[angle]}"
    
    # Verify atom coordinates are within reasonable bounds
    x_coords = [atom.get_x() for atom in result.atoms]
    y_coords = [atom.get_y() for atom in result.atoms]
    z_coords = [atom.get_z() for atom in result.atoms]
    
    # Coordinates should generally be within cell dimensions (allowing for some flexibility)
    x_range = max(x_coords) - min(x_coords)
    y_range = max(y_coords) - min(y_coords)
    z_range = max(z_coords) - min(z_coords)
    
    assert x_range <= cryst1_info['a'] + 10, \
        f"X coordinate range {x_range:.1f} exceeds cell dimension a={cryst1_info['a']}"
    assert y_range <= cryst1_info['b'] + 10, \
        f"Y coordinate range {y_range:.1f} exceeds cell dimension b={cryst1_info['b']}"
    assert z_range <= cryst1_info['c'] + 10, \
        f"Z coordinate range {z_range:.1f} exceeds cell dimension c={cryst1_info['c']}"
    
    # Verify reasonable coordinate distribution (no clustering at origin)
    center_x = sum(x_coords) / len(x_coords)
    center_y = sum(y_coords) / len(y_coords)
    center_z = sum(z_coords) / len(z_coords)
    
    # Center should be roughly in middle of cell
    assert 20 <= center_x <= cryst1_info['a'] - 20, \
        f"X center {center_x:.1f} suggests poor distribution in cell"
    assert 20 <= center_y <= cryst1_info['b'] - 20, \
        f"Y center {center_y:.1f} suggests poor distribution in cell"
    assert 20 <= center_z <= cryst1_info['c'] - 20, \
        f"Z center {center_z:.1f} suggests poor distribution in cell"
    
    print(f"✓ Crystal cell validated: {cryst1_info['a']:.1f}×{cryst1_info['b']:.1f}×{cryst1_info['c']:.1f} Å³")


def test_4wp7_protein_structure_analysis():
    """Test detailed protein structure analysis and residue connectivity."""
    pdb_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.prod.74.rec.pdb")
    
    # Skip test if file doesn't exist
    if not os.path.exists(pdb_path):
        pytest.skip(f"Test file {pdb_path} not found")
    
    result = PDBParser.parse_file(pdb_path)
    assert result is not None, "Failed to parse 4wp7 PDB file"
    
    # Standard amino acid types
    standard_amino_acids = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
    }
    
    # Special residue types in this system
    special_residues = {"HSD"}  # Protonated histidine
    
    protein_residues = {}
    protein_atoms = []
    
    # Collect protein atoms and residues
    for atom in result.atoms:
        if atom.residue_name() in standard_amino_acids or atom.residue_name() in special_residues:
            protein_atoms.append(atom)
            res_key = f"{atom.residue_name()}_{atom.residue_number()}"
            if res_key not in protein_residues:
                protein_residues[res_key] = []
            protein_residues[res_key].append(atom)
    
    assert len(protein_atoms) > 28000, f"Expected >28k protein atoms, got {len(protein_atoms)}"
    
    # Verify protein residue sequence continuity
    residue_numbers = sorted(set(atom.residue_number() for atom in protein_atoms))
    
    # Check for reasonable sequence continuity (allowing some gaps)
    min_res = min(residue_numbers)
    max_res = max(residue_numbers)
    expected_range = max_res - min_res + 1
    actual_count = len(residue_numbers)
    
    # Should have most residues in the range (allowing ~10% gaps)
    assert actual_count >= expected_range * 0.9, \
        f"Protein sequence has too many gaps: {actual_count}/{expected_range} residues"
    
    # Verify standard backbone atoms for amino acids
    backbone_atoms = {"N", "CA", "C", "O"}
    
    # Check a sample of protein residues for backbone completeness
    sample_residues = list(protein_residues.keys())[:20]  # Check first 20 residues
    
    for res_key in sample_residues:
        atoms_in_residue = protein_residues[res_key]
        atom_names = {atom.atom_name() for atom in atoms_in_residue}
        res_name = atoms_in_residue[0].residue_name()
        
        # All non-proline residues should have N, CA, C, O
        if res_name != "PRO":  # Proline has different N
            missing_backbone = backbone_atoms - atom_names
            assert len(missing_backbone) == 0, \
                f"Residue {res_key} missing backbone atoms: {missing_backbone}"
        
        # All residues should have CA and C
        essential_atoms = {"CA", "C"}
        missing_essential = essential_atoms - atom_names
        assert len(missing_essential) == 0, \
            f"Residue {res_key} missing essential atoms: {missing_essential}"
    
    # Verify special residue HSD (protonated histidine)
    hsd_residues = [res for res in protein_residues.keys() if res.startswith("HSD")]
    if hsd_residues:
        # HSD should have characteristic atoms
        hsd_sample = protein_residues[hsd_residues[0]]
        hsd_atom_names = {atom.get_name() for atom in hsd_sample}
        
        # HSD should have histidine ring atoms
        expected_his_atoms = {"CG", "ND1", "CD2", "CE1", "NE2"}
        found_his_atoms = expected_his_atoms.intersection(hsd_atom_names)
        assert len(found_his_atoms) >= 3, \
            f"HSD residue should have histidine ring atoms, found: {found_his_atoms}"
    
    print(f"✓ Protein structure validated: {len(protein_residues)} residues, {len(protein_atoms)} atoms")


def test_4wp7_chemical_composition():
    """Test chemical composition and atom type validation."""
    pdb_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.prod.74.rec.pdb")
    
    # Skip test if file doesn't exist
    if not os.path.exists(pdb_path):
        pytest.skip(f"Test file {pdb_path} not found")
    
    result = PDBParser.parse_file(pdb_path)
    assert result is not None, "Failed to parse 4wp7 PDB file"
    
    # Count atom types by element
    element_counts = {}
    
    for atom in result.atoms:
        # Extract element from atom name (last character usually indicates element)
        atom_name = atom.atom_name().strip()
        
        # Determine element from atom name
        if atom_name.startswith(('C', 'CA', 'CB', 'CG', 'CD', 'CE', 'CZ')):
            element = 'C'
        elif atom_name.startswith(('N', 'ND', 'NE', 'NH', 'NZ')):
            element = 'N'
        elif atom_name.startswith(('O', 'OD', 'OE', 'OG', 'OH', 'OW')):
            element = 'O'
        elif atom_name.startswith(('H', 'HA', 'HB', 'HG', 'HD', 'HE', 'HZ', 'HW')):
            element = 'H'
        elif atom_name.startswith('S'):
            element = 'S'
        else:
            element = 'Other'
        
        element_counts[element] = element_counts.get(element, 0) + 1
    
    total_atoms = len(result.atoms)
    
    # Verify reasonable elemental composition for biomolecular system
    # Expected roughly: H > O > C > N (due to water dominance)
    
    assert 'H' in element_counts, "No hydrogen atoms found"
    assert 'O' in element_counts, "No oxygen atoms found"
    assert 'C' in element_counts, "No carbon atoms found"
    assert 'N' in element_counts, "No nitrogen atoms found"
    
    # Hydrogen should be most abundant (from water and organic molecules)
    h_fraction = element_counts['H'] / total_atoms
    assert 0.4 <= h_fraction <= 0.7, \
        f"Hydrogen fraction {h_fraction:.3f} unexpected for biomolecular system"
    
    # Oxygen should be second most abundant (from water)
    o_fraction = element_counts['O'] / total_atoms
    assert 0.2 <= o_fraction <= 0.4, \
        f"Oxygen fraction {o_fraction:.3f} unexpected for water-dominated system"
    
    # Carbon should be substantial (from protein and GCMC molecules)
    c_fraction = element_counts['C'] / total_atoms
    assert 0.05 <= c_fraction <= 0.2, \
        f"Carbon fraction {c_fraction:.3f} unexpected for organic system"
    
    # Nitrogen should be present but less abundant
    n_fraction = element_counts['N'] / total_atoms
    assert 0.01 <= n_fraction <= 0.1, \
        f"Nitrogen fraction {n_fraction:.3f} unexpected for biomolecular system"
    
    # Verify occupancy and B-factor fields are reasonable
    occupancies = [atom.occupancy for atom in result.atoms if hasattr(atom, 'occupancy')]
    b_factors = [atom.b_factor for atom in result.atoms if hasattr(atom, 'b_factor')]
    
    if occupancies:
        # Most occupancies should be 1.0 for crystal/simulation structures
        avg_occupancy = sum(occupancies) / len(occupancies)
        assert 0.8 <= avg_occupancy <= 1.0, \
            f"Average occupancy {avg_occupancy:.3f} suggests incomplete structure"
    
    if b_factors:
        # B-factors should be reasonable (typically 0-100 for simulation data)
        avg_b_factor = sum(b_factors) / len(b_factors)
        assert 0 <= avg_b_factor <= 100, \
            f"Average B-factor {avg_b_factor:.1f} outside expected range"
    
    print(f"✓ Chemical composition validated: H={h_fraction:.1%}, O={o_fraction:.1%}, C={c_fraction:.1%}, N={n_fraction:.1%}")


