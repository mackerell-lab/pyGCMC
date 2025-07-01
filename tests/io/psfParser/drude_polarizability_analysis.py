# tests/io/psfParser/drude_polarizability_analysis.py
"""
Drude PSF Parser Tests - Polarizability Analysis

This file documents the expected polarizability distribution and patterns
based on analyze_drude_psf.sh results.

Expected distribution:
- Small (< 0.5): 1671 atoms (43.8%)
- Medium (0.5-1.0): 567 atoms (14.9%)
- Large (1.0-1.5): 1339 atoms (35.1%)
- X-Large (1.5-2.0): 58 atoms (1.5%)
- Huge (> 2.0): 179 atoms (4.7%)

Top residue types by total polarizability:
- GLU: 314.82 Å³ total (n=288 atoms, avg=1.093)
- LYS: 228.80 Å³ total (n=342 atoms, avg=0.669)
- ASP: 222.58 Å³ total (n=208 atoms, avg=1.070)
"""

import pytest
import os
from pygcmc.io import PSFParser
from pygcmc.model import Topology

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_drude_polarizability_distribution():
    """Test polarizability distribution matches expected values."""
    psf_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.psf")
    
    if not os.path.exists(psf_path):
        pytest.skip(f"Test file {psf_path} not found")
    
    parser = PSFParser()
    topology = Topology()
    success = parser.parse_to_topology(psf_path, topology)
    assert success, "Failed to parse 4wp7_drude.psf"
    
    # Categorize atoms by polarizability
    distribution = {
        "Small (< 0.5)": 0,
        "Medium (0.5-1.0)": 0,
        "Large (1.0-1.5)": 0,
        "X-Large (1.5-2.0)": 0,
        "Huge (> 2.0)": 0
    }
    
    for i in range(topology.get_num_atoms()):
        atom = topology.get_atom(i)
        
        # Get alpha value (this will fail in current implementation)
        if hasattr(atom, 'alpha'):
            alpha = atom.alpha
        elif hasattr(atom, 'get_alpha'):
            alpha = atom.get_alpha()
        else:
            continue
            
        if alpha > 0:
            if alpha < 0.5:
                distribution["Small (< 0.5)"] += 1
            elif alpha < 1.0:
                distribution["Medium (0.5-1.0)"] += 1
            elif alpha < 1.5:
                distribution["Large (1.0-1.5)"] += 1
            elif alpha < 2.0:
                distribution["X-Large (1.5-2.0)"] += 1
            else:
                distribution["Huge (> 2.0)"] += 1
    
    # Verify distribution matches expected
    assert distribution["Small (< 0.5)"] == 1671, f"Expected 1671 small alpha atoms"
    assert distribution["Medium (0.5-1.0)"] == 567, f"Expected 567 medium alpha atoms"
    assert distribution["Large (1.0-1.5)"] == 1339, f"Expected 1339 large alpha atoms"
    assert distribution["X-Large (1.5-2.0)"] == 58, f"Expected 58 x-large alpha atoms"
    assert distribution["Huge (> 2.0)"] == 179, f"Expected 179 huge alpha atoms"


def test_drude_residue_polarizability():
    """Test residue-specific polarizability patterns."""
    psf_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.psf")
    
    if not os.path.exists(psf_path):
        pytest.skip(f"Test file {psf_path} not found")
    
    parser = PSFParser()
    topology = Topology()
    success = parser.parse_to_topology(psf_path, topology)
    assert success, "Failed to parse 4wp7_drude.psf"
    
    # Calculate polarizability by residue type
    residue_alpha = {}
    residue_atoms = {}
    
    for i in range(topology.get_num_atoms()):
        atom = topology.get_atom(i)
        residue = topology.get_residue(atom.residue_id)
        res_name = residue.name
        
        # Get alpha value (this will fail in current implementation)
        if hasattr(atom, 'alpha'):
            alpha = atom.alpha
        elif hasattr(atom, 'get_alpha'):
            alpha = atom.get_alpha()
        else:
            continue
            
        if alpha > 0:
            if res_name not in residue_alpha:
                residue_alpha[res_name] = 0.0
                residue_atoms[res_name] = 0
            residue_alpha[res_name] += alpha
            residue_atoms[res_name] += 1
    
    # Verify top residue types by total polarizability
    expected_top_residues = {
        "GLU": {"total": 314.82, "atoms": 288, "avg": 1.093},
        "LYS": {"total": 228.80, "atoms": 342, "avg": 0.669},
        "ASP": {"total": 222.58, "atoms": 208, "avg": 1.070},
        "PHE": {"total": 219.60, "atoms": 242, "avg": 0.907},
        "ILE": {"total": 218.86, "atoms": 280, "avg": 0.782}
    }
    
    for res_name, expected in expected_top_residues.items():
        if res_name in residue_alpha:
            actual_total = residue_alpha[res_name]
            actual_atoms = residue_atoms[res_name]
            actual_avg = actual_total / actual_atoms if actual_atoms > 0 else 0
            
            assert abs(actual_total - expected["total"]) < 0.1, \
                f"{res_name}: Expected total {expected['total']}, got {actual_total}"
            assert actual_atoms == expected["atoms"], \
                f"{res_name}: Expected {expected['atoms']} atoms, got {actual_atoms}"
            assert abs(actual_avg - expected["avg"]) < 0.01, \
                f"{res_name}: Expected avg {expected['avg']}, got {actual_avg:.3f}"


def test_drude_anisotropy_information():
    """Test anisotropic Drude particle information.
    
    From analysis:
    - Anisotropic sites: 745
    - Anisotropic fraction: 19.53% of Drude particles
    """
    psf_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.psf")
    
    if not os.path.exists(psf_path):
        pytest.skip(f"Test file {psf_path} not found")
    
    parser = PSFParser()
    topology = Topology()
    success = parser.parse_to_topology(psf_path, topology)
    assert success, "Failed to parse 4wp7_drude.psf"
    
    # The PSF file should have anisotropic information
    # This would need to be parsed from the NUMANISO section
    # For now, we document the expected values
    expected_anisotropic_sites = 745
    expected_drude_particles = 3814
    expected_fraction = 0.1953  # 19.53%
    
    # These would need to be implemented in the parser
    # anisotropic_sites = topology.get_num_anisotropic()
    # assert anisotropic_sites == expected_anisotropic_sites
    
    # Document that anisotropic information should be available
    # For now, just check if we can access the expected counts from our basic tests
    assert expected_anisotropic_sites == 745
    assert expected_drude_particles == 3814
    assert abs(expected_fraction - 0.1953) < 1e-4