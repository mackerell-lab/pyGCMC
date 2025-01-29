# tests/io/test_pdbParser.py

import pytest
import os
import pygcmc
import math

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

def test_parse_simple_pdb():
    """Test parsing a simple PDB file with basic ATOM records."""
    pdb_path = os.path.join(TEST_DATA_DIR, "simple.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Check if atoms were parsed
    assert len(result.atoms) == 8  # MET has 8 atoms
    assert len(result.residues) == 1  # One MET residue
    
    # Check first atom's properties (N)
    atom = result.atoms[0]
    assert atom.getBynu() == 1
    assert atom.getType() == "N"
    assert atom.getResname() == "MET"
    assert atom.getChain() == "A"
    assert atom.getIres() == 1
    
    # Check coordinates
    coords = atom.getCoor()
    assert len(coords) == 3
    assert math.isclose(coords[0], 27.340, rel_tol=1e-5)
    assert math.isclose(coords[1], 24.430, rel_tol=1e-5)
    assert math.isclose(coords[2], 2.614, rel_tol=1e-5)

def test_parse_hetatm():
    """Test parsing HETATM records."""
    pdb_path = os.path.join(TEST_DATA_DIR, "hetatm.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Check counts
    protein_atoms = [atom for atom in result.atoms if not atom.isHetatm()]
    water_atoms = [atom for atom in result.atoms if atom.isHetatm()]
    assert len(protein_atoms) == 5  # ALA has 5 atoms
    assert len(water_atoms) == 3    # 3 water molecules
    
    # Check water properties
    water = water_atoms[0]
    assert water.isHetatm()
    assert water.getResname() == "HOH"
    assert water.getType() == "O"
    
    # Check water coordinates
    coords = water.getCoor()
    assert math.isclose(coords[0], 15.168, rel_tol=1e-5)
    assert math.isclose(coords[1], 20.391, rel_tol=1e-5)
    assert math.isclose(coords[2], 18.649, rel_tol=1e-5)

def test_parse_ter():
    """Test parsing TER records and chain termination."""
    pdb_path = os.path.join(TEST_DATA_DIR, "multichain.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Check chain counts
    chains = set(atom.getChain() for atom in result.atoms)
    assert chains == {"A", "B", "C"}
    
    # Check residues in different chains
    chain_residues = {}
    for res in result.residues:
        chain = res.getChain()
        if chain not in chain_residues:
            chain_residues[chain] = []
        chain_residues[chain].append(res)
    
    # Verify chain contents
    assert len(chain_residues["A"]) == 1  # GLY
    assert len(chain_residues["B"]) == 1  # ALA
    assert len(chain_residues["C"]) == 1  # VAL
    
    # Check residue types
    assert chain_residues["A"][0].getResname() == "GLY"
    assert chain_residues["B"][0].getResname() == "ALA"
    assert chain_residues["C"][0].getResname() == "VAL"

def test_parse_secondary_structure():
    """Test parsing HELIX and SHEET records."""
    pdb_path = os.path.join(TEST_DATA_DIR, "secondary.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Check HELIX records
    assert len(result.helices) == 1  # Chain A has helices
    assert len(result.helices["A"]) == 2  # Two helices in chain A
    
    # Check helix classes
    helices = result.helices["A"]
    assert helices[0].helixClass == 1  # Right-handed alpha
    assert helices[1].helixClass == 1  # Right-handed alpha
    
    # Check SHEET records
    assert len(result.sheets) == 1  # Chain B has sheets
    assert len(result.sheets["B"]) == 2  # Two strands in chain B
    
    # Check sheet info format
    sheet_info = result.sheets["B"][1]  # Second strand
    parts = sheet_info.split(":")
    assert len(parts) == 3
    assert parts[0] == "S1"  # Sheet ID
    assert parts[1] == "2"   # Strand number
    assert parts[2] == "-1"  # Anti-parallel sense

def test_parse_ssbond():
    """Test parsing SSBOND records."""
    pdb_path = os.path.join(TEST_DATA_DIR, "ssbond.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)

    # Check SSBOND records
    assert len(result.ssbonds) == 2

    # Check first bond (intra-chain)
    bond1 = result.ssbonds[0]
    parts = bond1.split("-")
    assert parts[0] == "A:3 "   # First CYS: chain A, residue 3
    assert parts[1] == "A:20 "  # Second CYS: chain A, residue 20

    # Check second bond (inter-chain)
    bond2 = result.ssbonds[1]
    parts = bond2.split("-")
    assert parts[0] == "A:15 "  # First CYS: chain A, residue 15
    assert parts[1] == "B:5 "   # Second CYS: chain B, residue 5

    # Check CYS residue structure
    cys_residues = [res for res in result.residues if res.getResname() == "CYS"]
    assert len(cys_residues) > 0
    assert cys_residues[0].getResname() == "CYS"

def test_parse_invalid_pdb():
    """Test handling of invalid PDB files."""
    # Test non-existent file
    with pytest.raises(RuntimeError):
        pygcmc.PDBParser.parse_file("nonexistent.pdb")
    
    # Test invalid atom record
    invalid_pdb = """
ATOM   INVALID  N   MET A   1      27.340  24.430   2.614  1.00  0.00
"""
    with pytest.raises(RuntimeError):
        pygcmc.PDBParser.parse_string(invalid_pdb)
    
    # Test invalid coordinates
    invalid_coords = """
ATOM      1  N   MET A   1      XXXXX  24.430   2.614  1.00  0.00           N  
"""
    with pytest.raises(RuntimeError):
        pygcmc.PDBParser.parse_string(invalid_coords)

def test_residue_atom_association():
    """Test if atoms are correctly associated with residues."""
    pdb_path = os.path.join(TEST_DATA_DIR, "simple.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Check MET residue
    met = result.residues[0]
    assert met.getResname() == "MET"
    assert len(met.getAtoms()) == 8
    
    # Check atom types in residue
    atom_types = [atom.getType().strip() for atom in met.getAtoms()]
    expected_types = ["N", "CA", "C", "O", "CB", "CG", "SD", "CE"]
    assert sorted(atom_types) == sorted(expected_types)
    
    # Check atom-residue consistency
    for atom in met.getAtoms():
        assert atom.getResname() == met.getResname()
        assert atom.getIres() == met.getIres()
        assert atom.getChain() == met.getChain()
        assert atom.getInscode() == met.getInscode()

def test_coordinate_parsing():
    """Test parsing of atomic coordinates."""
    pdb_path = os.path.join(TEST_DATA_DIR, "simple.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    for atom in result.atoms:
        coords = atom.getCoor()
        # Check coordinate format
        assert len(coords) == 3
        assert all(isinstance(c, float) for c in coords)
        assert all(not isinstance(c, str) for c in coords)
        # Check coordinate ranges
        assert all(-1000 < c < 1000 for c in coords)
        # Check precision (PDB format: 8.3f)
        for c in coords:
            assert abs(c - round(c, 3)) < 1e-3

def test_occupancy_and_tempfactor():
    """Test parsing of occupancy and temperature factor."""
    pdb_path = os.path.join(TEST_DATA_DIR, "simple.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    for atom in result.atoms:
        # Check occupancy (default 1.00)
        assert math.isclose(atom.getOccupancy(), 1.00, rel_tol=1e-5)
        # Check temperature factor (default 0.00)
        assert math.isclose(atom.getTempfactor(), 0.00, rel_tol=1e-5)

def test_parse_protein_fragment():
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    # Check total number of protein atoms
    protein_atoms = [atom for atom in result.atoms if not atom.isHetatm()]
    assert len(protein_atoms) == 196, f"Expected 196 protein atoms, got {len(protein_atoms)}"
    
    # Check number of residues (only protein residues)
    protein_residues = set((atom.getResname(), atom.getIres()) 
                          for atom in protein_atoms 
                          if atom.getResname() in ["ALA", "VAL", "PRO", "ASN", "GLN"])
    assert len(protein_residues) == 9, f"Expected 9 protein residues, got {len(protein_residues)}"
    
    # Check first residue (ALA 7)
    ala_atoms = [atom for atom in protein_atoms if atom.getResname() == "ALA" and atom.getIres() == 7]
    assert len(ala_atoms) == 12, f"Expected 12 atoms in ALA 7, got {len(ala_atoms)}"
    
    # Check first atom properties
    first_atom = protein_atoms[0]
    assert first_atom.getType() == "N"
    assert first_atom.getResname() == "ALA"
    assert first_atom.getIres() == 7
    assert first_atom.getChain() == " "
    assert not first_atom.isHetatm()

def test_parse_crystal_info():
    """Test parsing crystallographic information."""
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Check box dimensions
    assert hasattr(result, "boxDimensions")
    assert len(result.boxDimensions) == 6
    assert math.isclose(result.boxDimensions[0], 127.022, rel_tol=1e-5)  # a
    assert math.isclose(result.boxDimensions[1], 133.419, rel_tol=1e-5)  # b
    assert math.isclose(result.boxDimensions[2], 132.854, rel_tol=1e-5)  # c
    assert math.isclose(result.boxDimensions[3], 90.0, rel_tol=1e-5)     # alpha
    assert math.isclose(result.boxDimensions[4], 90.0, rel_tol=1e-5)     # beta
    assert math.isclose(result.boxDimensions[5], 90.0, rel_tol=1e-5)     # gamma

def test_parse_solvent_and_ligands():
    result = pygcmc.PDBParser.parse_file(os.path.join(TEST_DATA_DIR, "test.pdb"))
    
    # Check water molecules (SOL)
    water_atoms = [atom for atom in result.atoms if atom.getResname() == "SOL"]
    assert len(water_atoms) == 30, f"Expected 30 water atoms (10 molecules), got {len(water_atoms)}"
    
    # Check BENX ligand (marked as ATOM in test.pdb)
    benx_atoms = [atom for atom in result.atoms if atom.getResname() == "BENX"]
    assert len(benx_atoms) == 13, f"Expected 13 BENX atoms, got {len(benx_atoms)}"
    
    # Check PRPX ligand (marked as ATOM in test.pdb)
    prpx_atoms = [atom for atom in result.atoms if atom.getResname() == "PRPX"]
    assert len(prpx_atoms) == 24, f"Expected 24 PRPX atoms (2 molecules), got {len(prpx_atoms)}"
    
    # Check first water molecule
    first_water = water_atoms[0]
    assert first_water.getType() == "OW"
    assert first_water.getResname() == "SOL"
    assert not first_water.isHetatm()  # In test.pdb, these are ATOM records
    
    # Check first BENX molecule
    first_benx = benx_atoms[0]
    assert first_benx.getType() == "CG"
    assert first_benx.getResname() == "BENX"
    assert not first_benx.isHetatm()  # In test.pdb, these are ATOM records
    
    # Check first PRPX molecule
    first_prpx = prpx_atoms[0]
    assert first_prpx.getType() == "H11"
    assert first_prpx.getResname() == "PRPX"
    assert not first_prpx.isHetatm()  # In test.pdb, these are ATOM records

def test_hydrogen_atoms():
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Get all atoms from first ALA residue (ALA 7)
    ala_atoms = [atom for atom in result.atoms if atom.getResname() == "ALA" and atom.getIres() == 7]
    h_atoms = [atom for atom in ala_atoms if atom.getType().strip()[0] == "H"]
    
    assert len(h_atoms) == 7, f"Expected 7 hydrogen atoms in ALA 7, got {len(h_atoms)}"
    
    # Check types of hydrogen atoms
    h_types = set(atom.getType().strip() for atom in h_atoms)
    expected_types = {"H1", "H2", "H3", "HA", "HB1", "HB2", "HB3"}
    assert h_types == expected_types, f"Expected H types {expected_types}, got {h_types}"
    
    # Check coordinates of first hydrogen
    h1_atom = next(atom for atom in h_atoms if atom.getType().strip() == "H1")
    coords = h1_atom.getCoor()
    assert len(coords) == 3
    assert abs(coords[0] - 76.044) < 0.001
    assert abs(coords[1] - 92.324) < 0.001
    assert abs(coords[2] - 93.379) < 0.001

def test_center_of_mass():
    """Test center of mass calculation for residues."""
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Test ALA 7 center of mass
    ala7 = next(res for res in result.residues if res.getResname() == "ALA" and res.getIres() == 7)
    
    # Debug information
    print("\nALA 7 atoms:")
    total_mass = 0.0
    weighted_pos = [0.0, 0.0, 0.0]
    for atom in ala7.getAtoms():
        coords = atom.getCoor()
        mass = atom.getMass()
        print(f"Atom {atom.getType()}: mass={mass:.3f}, coords=({coords[0]:.3f}, {coords[1]:.3f}, {coords[2]:.3f})")
        total_mass += mass
        for i in range(3):
            weighted_pos[i] += mass * coords[i]
    
    if total_mass > 0:
        for i in range(3):
            weighted_pos[i] /= total_mass
    print(f"\nCalculated COM: ({weighted_pos[0]:.3f}, {weighted_pos[1]:.3f}, {weighted_pos[2]:.3f})")
    
    # Get the COM from the residue
    com = ala7.getCenterOfMass()
    print(f"Residue COM: ({com[0]:.3f}, {com[1]:.3f}, {com[2]:.3f})")
    
    # Verify COM coordinates
    assert len(com) == 3
    # The COM should be somewhere in the middle of the residue
    assert 77.0 < com[0] < 79.0, f"COM x-coordinate {com[0]} not in range (77.0, 79.0)"
    assert 91.0 < com[1] < 93.0, f"COM y-coordinate {com[1]} not in range (91.0, 93.0)"
    assert 92.0 < com[2] < 94.0, f"COM z-coordinate {com[2]} not in range (92.0, 94.0)"
    
    # Test VAL 8 center of mass
    val8 = next(res for res in result.residues if res.getResname() == "VAL" and res.getIres() == 8)
    com = val8.getCenterOfMass()
    
    # Verify COM coordinates
    assert len(com) == 3
    assert 78.0 < com[0] < 80.0, f"VAL8 COM x-coordinate {com[0]} not in range (78.0, 80.0)"
    assert 94.0 < com[1] < 96.0, f"VAL8 COM y-coordinate {com[1]} not in range (94.0, 96.0)"
    assert 89.0 < com[2] < 91.0, f"VAL8 COM z-coordinate {com[2]} not in range (89.0, 91.0)"
    
    # Test water molecule (SOL)
    sol = next(res for res in result.residues if res.getResname() == "SOL")
    com = sol.getCenterOfMass()
    
    # Verify COM coordinates for water
    assert len(com) == 3
    # The COM should be close to the oxygen atom position for water
    assert all(isinstance(x, float) for x in com)
    assert all(not math.isnan(x) for x in com)
    assert all(not math.isinf(x) for x in com)

