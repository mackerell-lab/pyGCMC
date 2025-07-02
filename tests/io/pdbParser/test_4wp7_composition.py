# tests/io/pdbParser/test_4wp7_composition.py
"""
4wp7 Composition Tests - Protein structure and chemical analysis.

This file tests protein structure analysis and chemical composition of 4wp7 PDB file
including:
- Protein structure analysis and residue connectivity
- Chemical composition and elemental distribution
"""

import pytest
import os
import time
import math
from pygcmc.io import PDBParser

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


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


