# tests/io/psfParser/drude_comprehensive_analysis.py
"""
Drude PSF Parser Tests - Hydrogen Bonding and Validation

This file tests hydrogen bonding and validation features from analyze_drude_psf.sh
that are not covered in other test files:
- Hydrogen bonding analysis
- PSF validation summary (PSF/PDB consistency)

Based on the output of analyze_drude_psf.sh for 4wp7_drude.psf
"""

import pytest
import os
from pygcmc.io import PSFParser
from pygcmc.model import Topology

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_hydrogen_bonding_analysis():
    """Test hydrogen bonding donor/acceptor counts.
    
    From analyze_drude_psf.sh HYDROGEN BONDING ANALYSIS:
    - Donors: 823
    - Acceptors: 727
    - Donor/Acceptor ratio: 1.13
    
    Top donor types by frequency (from script):
    - ND2A2: 470 donors
    - ND3P3A: 114 donors
    - ND2P1A: 85 donors
    - ND2A1: 74 donors
    - OD31A: 54 donors
    
    Top acceptor types by frequency (corrected - script has a bug):
    - OD2C1A: 530 acceptors (all OD2C1A atoms)
    - OD2C2A: 116 acceptors
    - OD31A: 54 acceptors
    - SD31B: 11 acceptors
    - ND2R5B: 8 acceptors
    
    Note: The script incorrectly counts both columns in acceptor pairs,
    leading to wrong statistics (e.g., showing ND3P3A: 438 when there
    are only 39 ND3P3A atoms total in the PSF).
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
    # Note: The analyze_drude_psf.sh script has a bug - it counts both columns in acceptor pairs
    # The correct counts based on actual PSF parsing are:
    expected_acceptor_types = {
        "OD2C1A": 530,  # All 530 OD2C1A atoms are acceptors
        "OD2C2A": 116,  # Second most common acceptor type
        "OD31A": 54,    # Third most common
        "SD31B": 11,    # Fourth most common
        "ND2R5B": 8     # Fifth most common
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
    - Mass conservation: PASS (Total mass: 54184.512 amu)
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