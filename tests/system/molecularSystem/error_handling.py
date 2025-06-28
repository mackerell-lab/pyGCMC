# tests/system/molecularSystem/error_handling.py

import pytest
from .helpers import *

def test_combine_null_inputs():
    """Test handling of null inputs."""
    mol_system = pygcmc.MolecularSystem()
    
    with pytest.raises(ValueError):
        mol_system.combine(None, None)
    
    with pytest.raises(ValueError):
        mol_system.combine(pygcmc.Structure(), None)
    
    with pytest.raises(ValueError):
        mol_system.combine(None, TOPParser.parse_string(""))

def test_combine_mismatched_data(test_structure):
    """Test handling of mismatched Structure and Topology data."""
    mol_system = pygcmc.MolecularSystem()
    
    # Create a topology with different number of atoms
    topology = TOPParser.parse_file(os.path.join(TEST_DATA_DIR, "mols", "sol.itp"))
    topology.add_atom("CA", "CT", 0.0, 12.01, "ALA", 1, "PROT")
    topology.add_atom("CB", "CT", 0.0, 12.01, "ALA", 1, "PROT")
    
    with pytest.raises(RuntimeError):
        mol_system.combine(test_structure, topology)

def test_combine_empty_data():
    """Test combining empty Structure and Topology objects."""
    mol_system = pygcmc.MolecularSystem()
    structure = pygcmc.Structure()
    topology = TOPParser.parse_string("")
    
    molecular = mol_system.combine(structure, topology)
    
    assert molecular.get_num_atoms() == 0
    assert molecular.get_num_residues() == 0
    assert molecular.get_num_segments() == 0
    assert molecular.get_num_bonds() == 0
    assert molecular.get_num_angles() == 0
    assert len(molecular.segment_map) == 0
    assert len(molecular.residue_map) == 0
    assert len(molecular.atom_map) == 0

