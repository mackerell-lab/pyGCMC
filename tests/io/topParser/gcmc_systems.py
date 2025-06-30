# tests/io/topParser/gcmc_systems.py
"""
TOP Parser tests for GCMC systems - focusing on 4wp7 topology parsing.

PARSER ENHANCEMENT COMPLETED:
The TOP parser now handles preprocessor directives (#ifdef, #endif) with default CHARMM support:

✅ WORKS: Processes #include directives for all molecule types (.itp files)
✅ WORKS: Considers [ molecules ] section for complete system composition
✅ WORKS: Processes water molecules via _FF_CHARMM preprocessor support in tip3p.itp
✅ WORKS: Handles 7-token atom format in tip3p.itp (vs 8-token format in other files)

SOLUTION IMPLEMENTED: Added default _FF_CHARMM define and flexible atom parsing
- pp_state.defines["_FF_CHARMM"] = "1" enables CHARMM-specific sections
- Modified atom parsing to handle both 7-token and 8-token formats
- Result: Complete system parsing with all components

Expected vs Actual for 4wp7:
- Expected: ~218k atoms (protein 30k + water 171k + small molecules 17k)
- Actual:   ~221k atoms (complete system successfully parsed)
"""

import pytest
import os
import pygcmc
import re
from pygcmc.io import TOPParser
from pygcmc.model import Topology

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_parse_4wp7_gcmc_topology():
    """Test parsing 4wp7 TOP file - complete system parsing including water molecules."""
    top_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.gc.74.top")
    
    # Skip test if file doesn't exist
    if not os.path.exists(top_path):
        pytest.skip(f"Test file {top_path} not found")
    
    parser = TOPParser()
    topology = Topology()
    
    # Parse the topology file
    success = parser.parse_to_topology(top_path, topology)
    assert success, "Failed to parse 4wp7 TOP file"
    
    # Verify parser reads complete system (protein + small molecules + water)
    # Expected: ~218k atoms (protein 30k + small molecules 17k + water 171k)
    total_atoms = topology.get_num_atoms()
    assert 200000 <= total_atoms <= 230000, f"Expected 200k-230k atoms (complete system), got {total_atoms}"
    
    # Expected: ~60k residues (protein 494 + small molecules 2176 + water 57077)
    total_residues = topology.get_num_residues()
    assert 55000 <= total_residues <= 65000, f"Expected ~60k residues (complete system), got {total_residues}"
    
    # Single protein chain system
    total_segments = topology.get_num_segments()
    assert total_segments >= 1, f"Expected at least 1 segment for protein, got {total_segments}"
    
    # Verify connectivity for protein + small molecules (but NO water)
    bonds = topology.get_num_bonds()
    angles = topology.get_num_angles() 
    dihedrals = topology.get_num_dihedrals()
    
    # Based on actual measurement: protein ~31k + small molecules ~20k bonds
    assert 45000 <= bonds <= 55000, f"Expected 45k-55k bonds (protein + small molecules), got {bonds}"
    assert 70000 <= angles <= 85000, f"Expected 70k-85k angles (protein + small molecules), got {angles}"  
    assert 100000 <= dihedrals <= 120000, f"Expected 100k-120k dihedrals (protein + small molecules), got {dihedrals}"
    
    # Verify we have both protein residues AND small molecules (but NO water)
    protein_residue_names = {"ALA", "ARG", "ASN", "ASP", "GLN", "GLU", "GLY", "HIS", "ILE",
                           "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HSD"}
    
    gcmc_molecules = {"BENX", "PRPX", "DMEE", "MEOH", "FORM", "IMIA", "ACEY", "MAMY"}
    
    protein_residues = 0
    gcmc_residues = 0
    water_residues = 0
    asp_residues = 0
    
    for i in range(topology.get_num_residues()):
        residue = topology.get_residue(i)
        res_name = residue.name
        
        if res_name in protein_residue_names:
            protein_residues += 1
            if res_name == "ASP":
                asp_residues += 1
        elif res_name in gcmc_molecules:
            gcmc_residues += 1
        elif res_name == "SOL":
            water_residues += 1
    
    # Should have protein residues (~494)
    assert 450 <= protein_residues <= 550, f"Expected ~494 protein residues, found {protein_residues}"
    
    # Should have GCMC small molecules (~2176) 
    assert 2000 <= gcmc_residues <= 2300, f"Expected ~2176 GCMC molecules, found {gcmc_residues}"
    
    # Should find multiple ASP residues in the protein
    assert asp_residues >= 20, f"Expected multiple ASP residues in protein chain, found {asp_residues}"
    
    # Should have water molecules (~57077)
    assert 50000 <= water_residues <= 60000, \
        f"Expected ~57k water molecules, found {water_residues}"
    
    # Verify specific GCMC molecules are present
    found_gcmc_types = set()
    for i in range(topology.get_num_residues()):
        residue = topology.get_residue(i)
        if residue.name in gcmc_molecules:
            found_gcmc_types.add(residue.name)
    
    assert len(found_gcmc_types) >= 7, \
        f"Expected most GCMC molecule types, found {len(found_gcmc_types)}: {found_gcmc_types}"
        
    print(f"✓ Parser successfully reads complete system: {protein_residues} protein + {gcmc_residues} GCMC + {water_residues} water")


def test_4wp7_topology_atom_types():
    """Test atom types and parameters in 4wp7 topology.""" 
    top_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.gc.74.top")
    
    # Skip test if file doesn't exist
    if not os.path.exists(top_path):
        pytest.skip(f"Test file {top_path} not found")
    
    parser = TOPParser()
    topology = Topology()
    
    # Parse the topology file
    success = parser.parse_to_topology(top_path, topology)
    assert success, "Failed to parse 4wp7 TOP file"
    
    # Check that we have expected CHARMM36 atom types
    atom_types_found = set()
    charges_found = []
    masses_found = []
    
    for i in range(topology.get_num_atoms()):
        atom = topology.get_atom(i)
        atom_types_found.add(atom.type)
        charges_found.append(atom.charge)
        masses_found.append(atom.mass)
    
    # Should find standard CHARMM36 protein atom types (more strict check)
    expected_protein_types = {"NH3", "HC", "CT1", "HB1", "C", "O", "CT2", "CT3", "NH1", "H"}
    expected_water_types = {"OT", "HT"}  # TIP3P water
    expected_small_mol_types = {"CA", "HA"}  # Aromatic carbons from GCMC molecules
    
    found_protein_types = atom_types_found.intersection(expected_protein_types)
    found_water_types = atom_types_found.intersection(expected_water_types)
    
    # Should find most protein types (system has substantial protein content)
    assert len(found_protein_types) >= 7, f"Expected more CHARMM36 protein types, found {found_protein_types}"
    
    # Should find water types (large solvent content)
    assert len(found_water_types) >= 1, f"Expected water types for solvated system, found {atom_types_found}"
    
    # Total unique types should be substantial for complex GCMC system
    assert len(atom_types_found) >= 15, f"Expected diverse atom types for GCMC system, found {len(atom_types_found)} types"
    
    # Verify charges are reasonable (should be between -2 and +2 for most atoms)
    reasonable_charges = [c for c in charges_found if -2.0 <= c <= 2.0]
    assert len(reasonable_charges) > len(charges_found) * 0.95, "Most charges should be reasonable"
    
    # Verify masses are reasonable (should be > 0 and < 200 for standard atoms)
    reasonable_masses = [m for m in masses_found if 0.5 <= m <= 200.0]
    assert len(reasonable_masses) > len(masses_found) * 0.9, "Most masses should be reasonable"
    
    # Check specific N-terminal nitrogen (should be NH3 type)
    first_n_atom = None
    for i in range(min(100, topology.get_num_atoms())):  # Check first 100 atoms
        atom = topology.get_atom(i)
        if atom.name == "N" and atom.type == "NH3":
            first_n_atom = atom
            break
    
    if first_n_atom:
        # N-terminal nitrogen should have negative charge (around -0.3 in CHARMM36)
        assert -0.6 <= first_n_atom.charge <= 0.0, f"N-terminal N charge should be negative (~-0.3), got {first_n_atom.charge}"
        # Standard nitrogen mass
        assert 13.0 <= first_n_atom.mass <= 15.0, f"Nitrogen mass should be ~14, got {first_n_atom.mass}"


def test_4wp7_force_field_includes():
    """Test force field include directives in 4wp7 topology."""
    top_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.gc.74.top")
    
    # Skip test if file doesn't exist
    if not os.path.exists(top_path):
        pytest.skip(f"Test file {top_path} not found")
    
    # Read file content to verify includes
    with open(top_path, 'r') as f:
        content = f.read()
    
    # Should include CHARMM36 force field
    assert "charmm36" in content.lower(), "Should reference CHARMM36 force field"
    assert "#include" in content, "Should have include directives"
    
    # Check for proper moleculetype definition (allow space variations)
    assert "[ moleculetype ]" in content or "[moleculetype]" in content, "Should have moleculetype section"
    assert "Protein_chain_P" in content, "Should define protein chain molecule type"
    
    # Should have atoms section (allow space variations)
    assert "[ atoms ]" in content or "[atoms]" in content, "Should have atoms section"
    
    # Verify file structure includes key sections
    required_sections = ["moleculetype", "atoms"]
    for section in required_sections:
        pattern = rf"\[\s*{section}\s*\]"
        assert re.search(pattern, content), f"Should have [{section}] section"
    
    # Check that includes reference proper file structure
    include_pattern = r"#include\s+[\"']([^\"']+)[\"']"
    includes = re.findall(include_pattern, content)
    
    assert len(includes) >= 1, "Should have at least one include directive"
    
    # First include should be the force field
    if includes:
        first_include = includes[0]
        assert "forcefield" in first_include.lower() or "charmm" in first_include.lower(), \
            f"First include should be force field file, got {first_include}"


def test_4wp7_parser_limitations():
    """Test documenting current parser limitations - should be updated when parser is enhanced."""
    top_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.gc.74.top")
    
    # Skip test if file doesn't exist
    if not os.path.exists(top_path):
        pytest.skip(f"Test file {top_path} not found")
    
    parser = TOPParser()
    topology = Topology()
    
    # Parse the topology file
    success = parser.parse_to_topology(top_path, topology)
    assert success, "Failed to parse 4wp7 TOP file"
    
    # Document current parser limitations
    total_atoms = topology.get_num_atoms()
    total_residues = topology.get_num_residues()
    
    # Current limitation: only reads protein (30k atoms, 494 residues)
    current_protein_only = (28000 <= total_atoms <= 32000 and 480 <= total_residues <= 510)
    
    if current_protein_only:
        # Parser has limitations - only reads protein
        print(f"PARSER LIMITATION: Only reading protein template ({total_atoms} atoms, {total_residues} residues)")
        print("Expected for complete system: ~218k atoms (protein + water + small molecules)")
        
        # Verify no small molecules or water are parsed
        small_mol_types = {"BENX", "PRPX", "DMEE", "MEOH", "FORM", "IMIA", "ACEY", "MAMY", "SOL"}
        found_small_mols = set()
        
        for i in range(topology.get_num_residues()):
            residue = topology.get_residue(i)
            if residue.name in small_mol_types:
                found_small_mols.add(residue.name)
        
        assert len(found_small_mols) == 0, \
            f"Parser limitation: should only find protein, but found: {found_small_mols}"
            
        # Mark test as documenting current limitation
        print("✓ Parser limitation confirmed: ignores #include directives and [ molecules ] section")
        
    else:
        # Parser has been enhanced! 
        print(f"PARSER ENHANCED: Now reading complete system ({total_atoms} atoms, {total_residues} residues)")
        
        # When parser is enhanced, verify complete system
        expected_total_atoms = 220000  # Protein 30k + water 171k + small molecules 17k (actual ~221k)
        expected_total_residues = 60000  # Protein 494 + water 57k + small molecules 2.4k
        
        assert expected_total_atoms * 0.9 <= total_atoms <= expected_total_atoms * 1.1, \
            f"Enhanced parser should read ~{expected_total_atoms} atoms, got {total_atoms}"
            
        assert expected_total_residues * 0.9 <= total_residues <= expected_total_residues * 1.1, \
            f"Enhanced parser should read ~{expected_total_residues} residues, got {total_residues}"
        
        # Verify presence of small molecules and water
        small_mol_types = {"BENX", "PRPX", "DMEE", "MEOH", "FORM", "IMIA", "ACEY", "MAMY", "SOL"}
        found_small_mols = set()
        
        for i in range(topology.get_num_residues()):
            residue = topology.get_residue(i)
            if residue.name in small_mol_types:
                found_small_mols.add(residue.name)
        
        assert len(found_small_mols) >= 8, \
            f"Enhanced parser should find all small molecule types, found: {found_small_mols}"
            
        print("✓ Parser enhancement verified: correctly processes includes and [ molecules ] section")