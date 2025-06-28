# tests/system/molecularSystem/basic_combination.py

import pytest
from .helpers import *

def test_combine_structure_data(test_structure, test_topology):
    """Test if Structure data is correctly combined into Molecular object."""
    mol_system = pygcmc.MolecularSystem()
    molecular = mol_system.combine(test_structure, test_topology)
    
    # Check structure data
    assert len(molecular.atoms) == len(test_structure.atoms)
    assert len(molecular.residues) == len(test_structure.residues)
    
    # Check specific atom properties
    first_atom = molecular.atoms[0]
    first_struct_atom = test_structure.atoms[0]
    assert first_atom.get_type() == first_struct_atom.get_type()
    assert first_atom.get_resname() == first_struct_atom.get_resname()
    assert first_atom.get_chain() == first_struct_atom.get_chain()
    
    # Check box dimensions
    assert molecular.boxDimensions == pytest.approx(test_structure.box_dimensions)

def test_combine_topology_data(test_structure, test_topology):
    """Test if Topology data is correctly combined into Molecular object."""
    mol_system = pygcmc.MolecularSystem()
    molecular = mol_system.combine(test_structure, test_topology)
    
    # Check topology data
    assert len(molecular.topology_atoms) == test_topology.get_num_atoms()
    assert len(molecular.topology_residues) == test_topology.get_num_residues()
    assert len(molecular.segments) == test_topology.get_num_segments()
    
    # Check specific topology properties
    first_top_atom = molecular.topology_atoms[0]
    orig_top_atom = test_topology.get_atom(0)
    assert first_top_atom.name == orig_top_atom.name
    assert first_top_atom.type == orig_top_atom.type
    assert first_top_atom.charge == pytest.approx(orig_top_atom.charge)
    assert first_top_atom.mass == pytest.approx(orig_top_atom.mass)

def test_combine_bonds_and_angles(test_structure, test_topology):
    """Test if bonds and angles are correctly combined into Molecular object."""
    mol_system = pygcmc.MolecularSystem()
    molecular = mol_system.combine(test_structure, test_topology)
    
    # Check bonds
    assert molecular.get_num_bonds() == test_topology.get_num_bonds()
    if molecular.get_num_bonds() > 0:
        first_bond = molecular.bonds[0]
        assert len(first_bond) == 2  # Bond should connect 2 atoms
        assert all(isinstance(idx, int) for idx in first_bond)
    
    # Check angles
    assert molecular.get_num_angles() == test_topology.get_num_angles()
    if molecular.get_num_angles() > 0:
        first_angle = molecular.angles[0]
        assert len(first_angle) == 3  # Angle should involve 3 atoms
        assert all(isinstance(idx, int) for idx in first_angle)

def test_combine_mapping_data(test_structure, test_topology):
    """Test if mapping data is correctly combined into Molecular object."""
    mol_system = pygcmc.MolecularSystem()
    molecular = mol_system.combine(test_structure, test_topology)
    
    # Check segment mapping
    assert len(molecular.segment_map) == test_topology.get_num_segments()
    
    # Check residue mapping
    assert len(molecular.residue_map) == test_topology.get_num_residues()
    
    # Check atom mapping
    assert len(molecular.atom_map) == test_topology.get_num_atoms()
    
    # Check specific mapping
    if test_topology.get_num_atoms() > 0:
        first_atom = test_topology.get_atom(0)
        first_res = test_topology.get_residue(first_atom.residue_id)
        key = (first_res.name, first_res.number, first_atom.name)
        assert key in molecular.atom_map
        assert molecular.atom_map[key] == 0

