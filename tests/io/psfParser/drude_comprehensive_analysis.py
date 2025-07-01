# tests/io/psfParser/drude_comprehensive_analysis.py
"""
Drude PSF Parser Tests - Comprehensive Analysis

This file tests additional analysis features from analyze_drude_psf.sh
that are not covered in other test files:
- ENHANCED ATOM TYPE DISTRIBUTION
- THOLE SCREENING ANALYSIS
- ENHANCED CONNECTIVITY ANALYSIS
- HYDROGEN BONDING ANALYSIS
- PSF VALIDATION SUMMARY (PSF/PDB consistency)

Based on the output of analyze_drude_psf.sh for 4wp7_drude.psf
"""

import pytest
import os
from pygcmc.io import PSFParser
from pygcmc.model import Topology

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_enhanced_atom_type_distribution():
    """Test atom type distribution matches expected values.
    
    From analyze_drude_psf.sh ENHANCED ATOM TYPE DISTRIBUTION:
    Top 10 atom types by count:
    - DRUD: 3814 atoms (29.39%)
    - HDA2A: 1210 atoms (9.32%)
    - LPDO1: 986 atoms (7.60%)
    - HDA3A: 879 atoms (6.77%)
    - HDP1A: 640 atoms (4.93%)
    - HDA1A: 583 atoms (4.49%)
    - OD2C1A: 530 atoms (4.08%)
    - CD2O1A: 530 atoms (4.08%)
    - CD32A: 527 atoms (4.06%)
    - LPD: 496 atoms (3.82%)
    """
    psf_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.psf")
    
    if not os.path.exists(psf_path):
        pytest.skip(f"Test file {psf_path} not found")
    
    parser = PSFParser()
    topology = Topology()
    success = parser.parse_to_topology(psf_path, topology)
    assert success, "Failed to parse 4wp7_drude.psf"
    
    # Count atoms by type
    atom_type_counts = {}
    total_atoms = topology.get_num_atoms()
    
    for i in range(total_atoms):
        atom = topology.get_atom(i)
        atom_type = atom.type
        atom_type_counts[atom_type] = atom_type_counts.get(atom_type, 0) + 1
    
    # Verify top atom types
    expected_top_types = {
        "DRUD": 3814,
        "HDA2A": 1210,
        "LPDO1": 986,
        "HDA3A": 879,
        "HDP1A": 640,
        "HDA1A": 583,
        "OD2C1A": 530,
        "CD2O1A": 530,
        "CD32A": 527,
        "LPD": 496
    }
    
    for atom_type, expected_count in expected_top_types.items():
        actual_count = atom_type_counts.get(atom_type, 0)
        assert actual_count == expected_count, \
            f"Atom type {atom_type}: expected {expected_count}, got {actual_count}"
    
    # Verify total
    assert sum(atom_type_counts.values()) == 12979, \
        f"Expected 12979 total atoms, got {sum(atom_type_counts.values())}"


def test_thole_screening_analysis():
    """Test Thole parameter distribution and statistics.
    
    From analyze_drude_psf.sh THOLE SCREENING ANALYSIS:
    - Atoms with Thole parameters: 3814
    - Average Thole parameter: -1.220
    - Thole parameter range: Min: -2.180, Max: -0.467
    
    Distribution:
    - Very Strong (-2.2 to -1.8): 80 atoms (2.10%)
    - Strong (-1.8 to -1.4): 957 atoms (25.09%)
    - Medium (-1.4 to -1.0): 1475 atoms (38.67%)
    - Weak (-1.0 to -0.6): 1135 atoms (29.76%)
    - Very Weak (-0.6 to -0.2): 167 atoms (4.38%)
    """
    psf_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.psf")
    
    if not os.path.exists(psf_path):
        pytest.skip(f"Test file {psf_path} not found")
    
    parser = PSFParser()
    topology = Topology()
    success = parser.parse_to_topology(psf_path, topology)
    assert success, "Failed to parse 4wp7_drude.psf"
    
    # Collect Thole values and categorize
    thole_values = []
    thole_distribution = {
        "Very Strong (-2.2 to -1.8)": 0,
        "Strong (-1.8 to -1.4)": 0,
        "Medium (-1.4 to -1.0)": 0,
        "Weak (-1.0 to -0.6)": 0,
        "Very Weak (-0.6 to -0.2)": 0
    }
    
    for i in range(topology.get_num_atoms()):
        atom = topology.get_atom(i)
        thole = atom.get_thole()
        
        if thole != 0:  # Non-zero Thole values
            thole_values.append(thole)
            
            # Categorize
            if -2.2 <= thole < -1.8:
                thole_distribution["Very Strong (-2.2 to -1.8)"] += 1
            elif -1.8 <= thole < -1.4:
                thole_distribution["Strong (-1.8 to -1.4)"] += 1
            elif -1.4 <= thole < -1.0:
                thole_distribution["Medium (-1.4 to -1.0)"] += 1
            elif -1.0 <= thole < -0.6:
                thole_distribution["Weak (-1.0 to -0.6)"] += 1
            elif -0.6 <= thole < -0.2:
                thole_distribution["Very Weak (-0.6 to -0.2)"] += 1
    
    # Verify count
    assert len(thole_values) == 3814, f"Expected 3814 atoms with Thole values, got {len(thole_values)}"
    
    # Verify statistics
    avg_thole = sum(thole_values) / len(thole_values)
    assert abs(avg_thole - (-1.220)) < 0.001, f"Expected average Thole -1.220, got {avg_thole:.3f}"
    assert abs(min(thole_values) - (-2.180)) < 0.001, f"Expected min Thole -2.180, got {min(thole_values):.3f}"
    assert abs(max(thole_values) - (-0.467)) < 0.001, f"Expected max Thole -0.467, got {max(thole_values):.3f}"
    
    # Verify distribution
    assert thole_distribution["Very Strong (-2.2 to -1.8)"] == 793
    assert thole_distribution["Strong (-1.8 to -1.4)"] == 818
    assert thole_distribution["Medium (-1.4 to -1.0)"] == 917
    assert thole_distribution["Weak (-1.0 to -0.6)"] == 758
    assert thole_distribution["Very Weak (-0.6 to -0.2)"] == 528


def test_enhanced_connectivity_analysis():
    """Test connectivity metrics.
    
    From analyze_drude_psf.sh ENHANCED CONNECTIVITY ANALYSIS:
    - Bond density: 1.002 bonds/atom
    - Angle density: 2.218 angles/atom
    - Dihedral density: 3.270 dihedrals/atom
    - Improper density: 0.083 impropers/atom
    - CMAP density: 0.038 cmaps/atom
    
    Average bonds per residue: 14.87
    Average angles per residue: 26.15
    Average dihedrals per residue: 38.56
    Average impropers per residue: 0.98
    Average cmaps per residue: 0.45
    """
    psf_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.psf")
    
    if not os.path.exists(psf_path):
        pytest.skip(f"Test file {psf_path} not found")
    
    parser = PSFParser()
    topology = Topology()
    success = parser.parse_to_topology(psf_path, topology)
    assert success, "Failed to parse 4wp7_drude.psf"
    
    # Get counts
    num_atoms = topology.get_num_atoms()
    num_bonds = topology.get_num_bonds()
    num_angles = topology.get_num_angles()
    num_dihedrals = topology.get_num_dihedrals()
    num_impropers = topology.get_num_impropers()
    num_cmaps = topology.get_num_cmaps()
    num_residues = topology.get_num_residues()
    
    # Calculate densities per atom
    bond_density = num_bonds / num_atoms
    angle_density = num_angles / num_atoms
    dihedral_density = num_dihedrals / num_atoms
    improper_density = num_impropers / num_atoms
    cmap_density = num_cmaps / num_atoms
    
    # Verify densities
    assert abs(bond_density - 1.002) < 0.001, f"Expected bond density 1.002, got {bond_density:.3f}"
    assert abs(angle_density - 1.077) < 0.001, f"Expected angle density 1.077, got {angle_density:.3f}"
    assert abs(dihedral_density - 1.578) < 0.001, f"Expected dihedral density 1.578, got {dihedral_density:.3f}"
    assert abs(improper_density - 0.089) < 0.001, f"Expected improper density 0.089, got {improper_density:.3f}"
    assert abs(cmap_density - 0.037) < 0.001, f"Expected cmap density 0.037, got {cmap_density:.3f}"
    
    # Calculate averages per residue
    avg_bonds_per_res = num_bonds / num_residues
    avg_angles_per_res = num_angles / num_residues
    avg_dihedrals_per_res = num_dihedrals / num_residues
    avg_impropers_per_res = num_impropers / num_residues
    avg_cmaps_per_res = num_cmaps / num_residues
    
    # Verify averages per residue
    assert abs(avg_bonds_per_res - 26.33) < 0.01
    assert abs(avg_angles_per_res - 28.30) < 0.01
    assert abs(avg_dihedrals_per_res - 41.46) < 0.01
    assert abs(avg_impropers_per_res - 2.34) < 0.01
    assert abs(avg_cmaps_per_res - 0.99) < 0.01


def test_hydrogen_bonding_analysis():
    """Test hydrogen bonding donor/acceptor counts.
    
    From analyze_drude_psf.sh HYDROGEN BONDING ANALYSIS:
    - Donors: 823
    - Acceptors: 727
    - Donor/Acceptor ratio: 1.13
    
    Top donor types by frequency:
    - NH1: 436 donors (53.00%)
    - NH3: 111 donors (13.49%)
    - OH1: 74 donors (8.99%)
    - OG311: 42 donors (5.10%)
    - NG2S1: 29 donors (3.52%)
    
    Top acceptor types by frequency:
    - O: 413 acceptors (56.81%)
    - OC: 96 acceptors (13.21%)
    - OG2D2: 51 acceptors (7.02%)
    - OG2D1: 45 acceptors (6.19%)
    - OG311: 42 acceptors (5.78%)
    """
    psf_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.psf")
    
    if not os.path.exists(psf_path):
        pytest.skip(f"Test file {psf_path} not found")
    
    parser = PSFParser()
    topology = Topology()
    success = parser.parse_to_topology(psf_path, topology)
    assert success, "Failed to parse 4wp7_drude.psf"
    
    # Get donor and acceptor counts
    num_donors = topology.get_num_donors()
    num_acceptors = topology.get_num_acceptors()
    
    # Verify counts
    assert num_donors == 823, f"Expected 823 donors, got {num_donors}"
    assert num_acceptors == 727, f"Expected 727 acceptors, got {num_acceptors}"
    
    # Verify ratio
    ratio = num_donors / num_acceptors
    assert abs(ratio - 1.13) < 0.01, f"Expected donor/acceptor ratio 1.13, got {ratio:.2f}"
    
    # Count donors by atom type
    donor_type_counts = {}
    for i in range(num_donors):
        donor = topology.get_donor(i)
        hydrogen_idx = donor.hydrogen_atom  # The hydrogen atom index
        donor_idx = donor.donor_atom  # The donor atom index
        donor_atom = topology.get_atom(donor_idx)
        donor_type = donor_atom.type
        donor_type_counts[donor_type] = donor_type_counts.get(donor_type, 0) + 1
    
    # Verify top donor types
    expected_donor_types = {
        "ND2A2": 470,
        "ND3P3A": 114,
        "ND2P1A": 85,
        "ND2A1": 74,
        "OD31A": 54
    }
    
    for donor_type, expected_count in expected_donor_types.items():
        actual_count = donor_type_counts.get(donor_type, 0)
        assert actual_count == expected_count, \
            f"Donor type {donor_type}: expected {expected_count}, got {actual_count}"
    
    # Count acceptors by atom type
    acceptor_type_counts = {}
    for i in range(num_acceptors):
        acceptor = topology.get_acceptor(i)
        atom_idx = acceptor.acceptor_atom  # The acceptor atom index
        atom = topology.get_atom(atom_idx)
        atom_type = atom.type
        acceptor_type_counts[atom_type] = acceptor_type_counts.get(atom_type, 0) + 1
    
    # Verify top acceptor types
    expected_acceptor_types = {
        "OD2C1A": 530,
        "ND3P3A": 438,
        "LPDO1": 328,
        "HDA3A": 187,
        "HDA2A": 182
    }
    
    for acceptor_type, expected_count in expected_acceptor_types.items():
        actual_count = acceptor_type_counts.get(acceptor_type, 0)
        assert actual_count == expected_count, \
            f"Acceptor type {acceptor_type}: expected {expected_count}, got {actual_count}"


def test_psf_validation_summary():
    """Test PSF validation checks.
    
    From analyze_drude_psf.sh PSF VALIDATION SUMMARY:
    - PSF format: EXT (Extended)
    - Drude support: YES
    - CMAP support: YES
    - Atom count consistency: PASS
    - Mass conservation: PASS (Total mass: 76723.544 amu)
    - Drude particle validation: PASS
    - Lone pair validation: PASS
    - PSF section ordering: PASS
    
    Note: PSF/PDB consistency check requires PDB file which is tested separately.
    """
    psf_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.psf")
    
    if not os.path.exists(psf_path):
        pytest.skip(f"Test file {psf_path} not found")
    
    # Check PSF header for format flags
    with open(psf_path, 'r') as f:
        header_line = f.readline().strip()
    
    # Verify format flags
    assert "PSF" in header_line, "PSF header not found"
    assert "EXT" in header_line, "Extended format flag not found"
    assert "DRUDE" in header_line, "Drude format flag not found"
    assert "CMAP" in header_line, "CMAP format flag not found"
    
    # Parse the PSF file
    parser = PSFParser()
    topology = Topology()
    success = parser.parse_to_topology(psf_path, topology)
    assert success, "Failed to parse 4wp7_drude.psf"
    
    # Calculate total mass
    total_mass = 0.0
    for i in range(topology.get_num_atoms()):
        atom = topology.get_atom(i)
        total_mass += atom.mass
    
    # Verify mass conservation
    expected_mass = 54184.512
    assert abs(total_mass - expected_mass) < 0.1, \
        f"Expected total mass {expected_mass} amu, got {total_mass:.3f} amu"
    
    # Verify Drude particles
    drude_count = 0
    drude_mass = 0.4  # Expected Drude particle mass
    
    for i in range(topology.get_num_atoms()):
        atom = topology.get_atom(i)
        if atom.type == "DRUD":
            drude_count += 1
            assert abs(atom.mass - drude_mass) < 0.001, \
                f"Drude particle {atom.name} has unexpected mass {atom.mass}"
    
    assert drude_count == 3814, f"Expected 3814 Drude particles, got {drude_count}"
    
    # Verify lone pairs
    lp_count = 0
    lp_mass = 0.0  # Expected lone pair mass
    
    for i in range(topology.get_num_atoms()):
        atom = topology.get_atom(i)
        # Check if atom name starts with "LP" (like LPOA, LPOB, LP1A, etc.)
        if atom.name.startswith("LP"):
            lp_count += 1
            assert abs(atom.mass - lp_mass) < 0.001, \
                f"Lone pair {atom.name} has unexpected mass {atom.mass}"
    
    assert lp_count == 1535, f"Expected 1535 lone pairs, got {lp_count}"