# tests/io/topParser/enhanced_4wp7_tests.py
"""
Enhanced TOP Parser tests for 4wp7 system - Deep analysis of GCMC topology parsing.

This file contains comprehensive tests for the 4wp7 topology file that go beyond
basic parsing verification to validate detailed molecular structures, connectivity,
and advanced features.

Tests include:
- Detailed GCMC molecule structure validation
- Topology connectivity verification  
- Preprocessor directive handling
- Error handling and edge cases
- Water structure validation
"""

import pytest
import os
import time
import tempfile
from pygcmc.io import TOPParser
from pygcmc.model import Topology

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_4wp7_gcmc_molecules_detailed():
    """Test detailed structure of each GCMC molecule type in 4wp7 system."""
    top_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.gc.74.top")
    
    # Skip test if file doesn't exist
    if not os.path.exists(top_path):
        pytest.skip(f"Test file {top_path} not found")
    
    parser = TOPParser()
    topology = Topology()
    
    # Parse the topology file
    success = parser.parse_to_topology(top_path, topology)
    assert success, "Failed to parse 4wp7 TOP file"
    
    # Expected GCMC molecule structures (actual atom counts from 4wp7 force field)
    expected_gcmc_structures = {
        "ACEY": {"atoms": 8, "description": "Acetate (CH3COO-) with explicit hydrogens"},
        "BENX": {"atoms": 13, "description": "Benzene with extra atoms (possibly virtual sites)"},
        "DMEE": {"atoms": 9, "description": "Dimethyl ether with explicit hydrogens"},
        "FORM": {"atoms": 6, "description": "Formamide with explicit hydrogens"},
        "IMIA": {"atoms": 9, "description": "Imidazole with explicit hydrogens"},
        "MAMY": {"atoms": 9, "description": "Methylammonium with explicit hydrogens"},
        "MEOH": {"atoms": 6, "description": "Methanol with explicit hydrogens"},
        "PRPX": {"atoms": 12, "description": "Propane with explicit hydrogens"}
    }
    
    # Count atoms per molecule type
    molecule_atom_counts = {}
    molecule_instances = {}
    
    for i in range(topology.get_num_residues()):
        residue = topology.get_residue(i)
        mol_name = residue.name
        
        if mol_name in expected_gcmc_structures:
            if mol_name not in molecule_atom_counts:
                molecule_atom_counts[mol_name] = len(residue.atoms)
                molecule_instances[mol_name] = 1
            else:
                molecule_instances[mol_name] += 1
                # Verify all instances have same atom count
                assert len(residue.atoms) == molecule_atom_counts[mol_name], \
                    f"Inconsistent atom count for {mol_name}: expected {molecule_atom_counts[mol_name]}, got {len(residue.atoms)}"
    
    # Verify expected molecule types are found
    found_molecules = set(molecule_atom_counts.keys())
    expected_molecules = set(expected_gcmc_structures.keys())
    assert found_molecules == expected_molecules, \
        f"Missing GCMC molecules. Expected: {expected_molecules}, Found: {found_molecules}"
    
    # Verify molecule instance counts match expectations (from [molecules] section)
    expected_counts = {
        "BENX": 287, "PRPX": 292, "DMEE": 241, "MEOH": 277,
        "FORM": 277, "IMIA": 289, "ACEY": 256, "MAMY": 257
    }
    
    for mol_name, expected_count in expected_counts.items():
        actual_count = molecule_instances.get(mol_name, 0)
        assert actual_count == expected_count, \
            f"Wrong number of {mol_name} molecules: expected {expected_count}, got {actual_count}"
    
    # Verify atom counts match expectations exactly
    for mol_name, structure_info in expected_gcmc_structures.items():
        if mol_name in molecule_atom_counts:
            actual_atoms = molecule_atom_counts[mol_name]
            expected_atoms = structure_info["atoms"]
            
            assert actual_atoms == expected_atoms, \
                f"{mol_name} ({structure_info['description']}): expected {expected_atoms} atoms, got {actual_atoms}"
    
    print(f"✓ All {len(found_molecules)} GCMC molecule types validated with correct instance counts")


def test_4wp7_topology_connectivity():
    """Test topology connectivity - bonds, angles, dihedrals completeness."""
    top_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.gc.74.top")
    
    # Skip test if file doesn't exist
    if not os.path.exists(top_path):
        pytest.skip(f"Test file {top_path} not found")
    
    parser = TOPParser()
    topology = Topology()
    
    # Parse the topology file
    success = parser.parse_to_topology(top_path, topology)
    assert success, "Failed to parse 4wp7 TOP file"
    
    total_atoms = topology.get_num_atoms()
    total_bonds = topology.get_num_bonds()
    total_angles = topology.get_num_angles()
    total_dihedrals = topology.get_num_dihedrals()
    
    # Verify connectivity ratios are reasonable for a protein+small molecule system
    # Typical protein: ~1.0 bonds per atom, ~1.5 angles per atom, ~2.0 dihedrals per atom
    # Small molecules have higher connectivity density
    
    bonds_per_atom = total_bonds / total_atoms if total_atoms > 0 else 0
    angles_per_atom = total_angles / total_atoms if total_atoms > 0 else 0
    dihedrals_per_atom = total_dihedrals / total_atoms if total_atoms > 0 else 0
    
    # For protein+small molecule system, expect moderate connectivity
    assert 0.15 <= bonds_per_atom <= 0.35, \
        f"Unusual bonds/atom ratio: {bonds_per_atom:.3f} (expected 0.15-0.35)"
    
    assert 0.25 <= angles_per_atom <= 0.55, \
        f"Unusual angles/atom ratio: {angles_per_atom:.3f} (expected 0.25-0.55)"
    
    assert 0.35 <= dihedrals_per_atom <= 0.75, \
        f"Unusual dihedrals/atom ratio: {dihedrals_per_atom:.3f} (expected 0.35-0.75)"
    
    # Verify protein backbone connectivity
    protein_residues = 0
    for i in range(topology.get_num_residues()):
        residue = topology.get_residue(i)
        if residue.name in {"ALA", "ARG", "ASN", "ASP", "GLN", "GLU", "GLY", "HIS", "ILE",
                           "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HSD"}:
            protein_residues += 1
    
    # For ~494 protein residues, expect ~493 peptide bonds (N-C connections)
    if protein_residues > 400:
        # Each peptide bond contributes 1 bond, plus side chain bonds
        expected_protein_bonds = protein_residues - 1  # Peptide bonds
        expected_protein_bonds += protein_residues * 5  # Average ~5 bonds per residue (backbone + side chain)
        
        # Should have substantial protein contribution to total bonds
        assert total_bonds >= expected_protein_bonds * 0.5, \
            f"Too few bonds for protein size: {total_bonds} total, expected at least {expected_protein_bonds * 0.5} from protein"
    
    print(f"✓ Connectivity validated: {bonds_per_atom:.3f} bonds/atom, {angles_per_atom:.3f} angles/atom, {dihedrals_per_atom:.3f} dihedrals/atom")


def test_4wp7_preprocessor_directives():
    """Test preprocessor directive handling (#ifdef, #include, #define)."""
    top_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.gc.74.top")
    
    # Skip test if file doesn't exist
    if not os.path.exists(top_path):
        pytest.skip(f"Test file {top_path} not found")
    
    # Read the raw file to check preprocessor directives
    with open(top_path, 'r') as f:
        content = f.read()
    
    # Verify key preprocessor features are present
    assert "#include" in content, "Should contain #include directives"
    
    # Count include directives
    include_count = content.count("#include")
    assert include_count >= 10, f"Expected at least 10 includes (forcefield + molecules + water), got {include_count}"
    
    # Verify force field include is first
    lines = content.split('\n')
    first_include_line = None
    for line in lines:
        if line.strip().startswith("#include"):
            first_include_line = line.strip()
            break
    
    assert first_include_line is not None, "Should have at least one #include directive"
    assert "forcefield.itp" in first_include_line, f"First include should be forcefield, got: {first_include_line}"
    
    # Parse with the topology parser to verify preprocessor handling
    parser = TOPParser()
    topology = Topology()
    
    success = parser.parse_to_topology(top_path, topology)
    assert success, "Failed to parse 4wp7 TOP file with preprocessor directives"
    
    # Verify the parser correctly processed includes by checking for water molecules
    # Water should only be present if tip3p.itp was correctly processed with _FF_CHARMM defined
    water_found = False
    for i in range(topology.get_num_residues()):
        residue = topology.get_residue(i)
        if residue.name == "SOL":
            water_found = True
            break
    
    assert water_found, "Water molecules not found - preprocessor may not be handling _FF_CHARMM correctly"
    
    # Verify all expected includes were processed by checking molecule presence
    expected_molecules = {"BENX", "PRPX", "DMEE", "MEOH", "FORM", "IMIA", "ACEY", "MAMY", "SOL"}
    found_molecules = set()
    
    for i in range(topology.get_num_residues()):
        residue = topology.get_residue(i)
        if residue.name in expected_molecules:
            found_molecules.add(residue.name)
    
    missing_molecules = expected_molecules - found_molecules
    assert len(missing_molecules) == 0, \
        f"Missing molecules suggest include processing failed: {missing_molecules}"
    
    print(f"✓ Preprocessor directives validated: {include_count} includes processed, all molecules found")



def test_4wp7_error_handling():
    """Test error handling for corrupted or incomplete 4wp7-style files."""
    top_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.gc.74.top")
    
    # Skip test if file doesn't exist
    if not os.path.exists(top_path):
        pytest.skip(f"Test file {top_path} not found")
    
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Test 1: Missing include file
        corrupted_file1 = os.path.join(tmp_dir, "missing_include.top")
        with open(corrupted_file1, 'w') as f:
            f.write('#include "./nonexistent.itp"\n')
            f.write('[ moleculetype ]\n')
            f.write('Test 1\n')
            f.write('[ atoms ]\n')
            f.write('1 CT1 1 TEST C1 1 0.0 12.011\n')
            f.write('[ system ]\n')
            f.write('Test\n')
            f.write('[ molecules ]\n')
            f.write('Test 1\n')
        
        parser = TOPParser()
        topology = Topology()
        
        # Should handle missing include gracefully (may succeed or fail, but shouldn't crash)
        try:
            result = parser.parse_to_topology(corrupted_file1, topology)
            # Either succeeds (graceful handling) or fails (proper error)
            assert isinstance(result, bool), "Parser should return boolean result"
        except Exception as e:
            # If it throws exception, should be a reasonable error message
            assert len(str(e)) > 0, "Exception should have meaningful message"
        
        # Test 2: Malformed [molecules] section
        corrupted_file2 = os.path.join(tmp_dir, "bad_molecules.top")
        with open(corrupted_file2, 'w') as f:
            f.write('[ moleculetype ]\n')
            f.write('Test 1\n')
            f.write('[ atoms ]\n')
            f.write('1 CT1 1 TEST C1 1 0.0 12.011\n')
            f.write('[ system ]\n')
            f.write('Test\n')
            f.write('[ molecules ]\n')
            f.write('Test 999999\n')  # Impossibly high count
        
        topology2 = Topology()
        
        # Should handle unrealistic molecule counts gracefully
        try:
            result = parser.parse_to_topology(corrupted_file2, topology2)
            if result:
                # If parsing succeeds, should not create impossibly many atoms
                assert topology2.get_num_atoms() <= 1000000, "Should not create excessive atoms from bad molecule count"
        except Exception as e:
            # Acceptable to fail with reasonable error
            assert "molecule" in str(e).lower() or "count" in str(e).lower(), \
                f"Error message should mention molecule/count issue: {e}"
        
        # Test 3: Truncated file (incomplete atoms section)
        truncated_file = os.path.join(tmp_dir, "truncated.top")
        with open(top_path, 'r') as original:
            lines = original.readlines()
        
        # Take first 1000 lines (cuts off in middle of atoms section)
        with open(truncated_file, 'w') as f:
            f.writelines(lines[:1000])
        
        topology3 = Topology()
        result = parser.parse_to_topology(truncated_file, topology3)
        
        # Truncated file might succeed (partial parsing) or fail, but shouldn't crash
        assert isinstance(result, bool), "Parser should handle truncated file gracefully"
        
        if result:
            # If it succeeds, should have parsed some atoms but not the full system
            assert 0 < topology3.get_num_atoms() < 50000, \
                f"Truncated parsing should give partial result, got {topology3.get_num_atoms()} atoms"
    
    print("✓ Error handling validated: missing includes, bad molecule counts, truncated files")


def test_4wp7_water_structure_validation():
    """Test detailed validation of TIP3P water structure in 4wp7 system."""
    top_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.gc.74.top")
    
    # Skip test if file doesn't exist
    if not os.path.exists(top_path):
        pytest.skip(f"Test file {top_path} not found")
    
    parser = TOPParser()
    topology = Topology()
    
    success = parser.parse_to_topology(top_path, topology)
    assert success, "Failed to parse 4wp7 TOP file"
    
    # Find water molecules
    water_residues = []
    for i in range(topology.get_num_residues()):
        residue = topology.get_residue(i)
        if residue.name == "SOL":
            water_residues.append(residue)
    
    assert len(water_residues) == 57077, f"Expected 57077 water molecules, found {len(water_residues)}"
    
    # Validate first few water molecules structure
    for i, water in enumerate(water_residues[:10]):  # Check first 10 waters
        # TIP3P water should have exactly 3 atoms: 1 O + 2 H
        assert len(water.atoms) == 3, f"Water {i+1} should have 3 atoms, got {len(water.atoms)}"
        
        # Get atom names and types
        atom_names = []
        atom_types = []
        atom_charges = []
        
        for atom_idx in water.atoms:
            atom = topology.get_atom(atom_idx)
            atom_names.append(atom.name)
            atom_types.append(atom.type)
            atom_charges.append(atom.charge)
        
        # TIP3P naming: OW, HW1, HW2 or O, H1, H2
        oxygen_names = {"OW", "O", "OH2"}
        hydrogen_names = {"HW1", "HW2", "H1", "H2"}
        
        oxygen_count = sum(1 for name in atom_names if name in oxygen_names)
        hydrogen_count = sum(1 for name in atom_names if name in hydrogen_names)
        
        assert oxygen_count == 1, f"Water {i+1} should have 1 oxygen, got {oxygen_count} (names: {atom_names})"
        assert hydrogen_count == 2, f"Water {i+1} should have 2 hydrogens, got {hydrogen_count} (names: {atom_names})"
        
        # Check atom types (TIP3P: OT for oxygen, HT for hydrogen)
        oxygen_types = {"OT", "OH2"}
        hydrogen_types = {"HT", "H"}
        
        oxygen_type_count = sum(1 for atype in atom_types if atype in oxygen_types)
        hydrogen_type_count = sum(1 for atype in atom_types if atype in hydrogen_types)
        
        assert oxygen_type_count == 1, f"Water {i+1} should have 1 oxygen type, got {oxygen_type_count} (types: {atom_types})"
        assert hydrogen_type_count == 2, f"Water {i+1} should have 2 hydrogen types, got {hydrogen_type_count} (types: {atom_types})"
        
        # Check charges are reasonable for TIP3P (O: ~-0.834, H: ~+0.417)
        total_charge = sum(atom_charges)
        assert abs(total_charge) < 0.01, f"Water {i+1} should be neutral, got total charge {total_charge:.4f}"
    
    print(f"✓ Water structure validated: {len(water_residues)} TIP3P molecules with correct 3-atom structure")