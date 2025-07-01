"""Shared fixtures for io tests - imports fixtures from subdirectories"""

import pytest
import os
from typing import Dict

# Directly import the analysis function from pdbParser conftest
def analyze_drude_structure_from_raw_file() -> Dict:
    """
    Analyze Drude structure directly from raw PDB file to get ground truth data.
    This avoids circular validation by not using pygcmc to get test expectations.
    """
    # Get the test data directory - we're in tests/io/
    test_data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
    pdb_path = os.path.join(test_data_dir, "4wp7", "4wp7_drude.pdb")
    
    if not os.path.exists(pdb_path):
        return {}
    
    stats = {
        'total_atoms': 0,
        'parent_atoms': 0, 
        'drude_particles': 0,
        'lone_pairs': 0,
        'residue_types': set(),
        'drude_atom_types': set(),
        'lone_pair_types': set(),
        'polarizable_residues': set()
    }
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                stats['total_atoms'] += 1
                atom_name = line[12:16].strip()
                res_name = line[17:21].strip()
                
                stats['residue_types'].add(res_name)
                
                if atom_name.startswith('D') and len(atom_name) > 1:
                    # Drude particle
                    stats['drude_particles'] += 1
                    stats['drude_atom_types'].add(atom_name)
                    stats['polarizable_residues'].add(res_name)
                elif atom_name.startswith('LP'):
                    # Lone pair
                    stats['lone_pairs'] += 1
                    stats['lone_pair_types'].add(atom_name)
                else:
                    # Parent atom
                    stats['parent_atoms'] += 1
    
    return stats


@pytest.fixture
def drude_ground_truth():
    """Fixture providing ground truth data from raw PDB file analysis"""
    return analyze_drude_structure_from_raw_file()