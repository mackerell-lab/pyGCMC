# tests/io/psfParser/drude_comprehensive_analysis.py
"""
Drude PSF Parser Tests - Atom Type and Connectivity Analysis

This file tests atom type distribution and connectivity features from analyze_drude_psf.sh
that are not covered in other test files:
- Enhanced atom type distribution
- Thole screening analysis
- Enhanced connectivity analysis

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
    - Very Strong (-2.2 to -1.8): 793 atoms (20.79%)
    - Strong (-1.8 to -1.4): 818 atoms (21.45%)
    - Medium (-1.4 to -1.0): 917 atoms (24.04%)
    - Weak (-1.0 to -0.6): 758 atoms (19.87%)
    - Very Weak (-0.6 to -0.2): 528 atoms (13.84%)
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
    - Angle density: 1.077 angles/atom
    - Dihedral density: 1.578 dihedrals/atom
    - Improper density: 0.089 impropers/atom
    - CMAP density: 0.037 cmaps/atom
    
    Average bonds per residue: 26.33
    Average angles per residue: 28.30
    Average dihedrals per residue: 41.46
    Average impropers per residue: 2.34
    Average cmaps per residue: 0.99
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


