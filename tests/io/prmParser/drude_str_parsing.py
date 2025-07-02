#!/usr/bin/env python3
"""
Test PRM Parser's ability to correctly read force field parameters from Drude STR files
"""

import pytest
import os
import sys

# Add pygcmc to path
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../..'))

import pygcmc


@pytest.fixture
def str_file_path():
    """Return Drude STR file path"""
    return os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 
        '../../../data/forcefields/charmm/c36_jul24/drude/drude_toppar_2023/toppar_drude_main_protein_2023a.str'
    )

def test_parse_drude_str_file(str_file_path):
    """Test parsing of Drude STR file"""
    # Use PRM Parser to parse STR file
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # Check basic information if successfully loaded
    assert ff is not None

def test_bond_parameters_from_str(str_file_path):
    """Test reading of bond parameters"""
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # Test specific bond parameters
    bond_params = ff.get_bond_params('ND2A2', 'CD2O1A')
    assert bond_params.kb == 376.20
    assert bond_params.b0 == 1.285
    
    bond_params = ff.get_bond_params('CD31C', 'CD32A')
    assert bond_params.kb == 222.50
    assert bond_params.b0 == 1.528  # Actual value from STR file
    
    bond_params = ff.get_bond_params('ODW', 'HDW')
    assert bond_params.kb == 450.00
    assert bond_params.b0 == 0.9572

def test_angle_parameters_from_str(str_file_path):
    """Test reading of angle parameters"""
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # Test specific angle parameters
    # Using parameters that appear first in ANGLES section
    angle_params = ff.get_angle_params('CD2O1A', 'ND2A2', 'CD31C')
    assert angle_params.ktheta == 40.90
    assert angle_params.theta0 == 116.10
    
    angle_params = ff.get_angle_params('HDP1A', 'ND2A1', 'HDP1A')
    assert angle_params.ktheta == 24.00
    assert angle_params.theta0 == 113.00

def test_dihedral_parameters_from_str(str_file_path):
    """Test reading of dihedral parameters"""
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # Test specific dihedral parameters
    dihedral_params = ff.get_dihedral_params('HDP1A', 'ND2A1', 'CD2O1A', 'OD2C1A')
    assert len(dihedral_params) > 0
    # Check first parameter
    assert dihedral_params[0].kchi == 2.000
    assert dihedral_params[0].n == 2
    assert dihedral_params[0].delta == 180.0

def test_improper_parameters_from_str(str_file_path):
    """Test reading of improper dihedral parameters"""
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # Test specific improper dihedral parameters
    # Using actual existing improper parameters
    improper_params = ff.get_improper_params('CD2O1A', 'CD32A', 'ND2A2', 'OD2C1A')
    assert improper_params.kpsi == 100.00
    assert improper_params.psi0 == 0.0

def test_nonbonded_parameters_from_str(str_file_path):
    """Test reading of nonbonded parameters"""
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # Test LJ parameters for specific atom types
    lj_params = ff.get_lj_params('HDP1A')
    assert lj_params.epsilon == -0.0100
    assert lj_params.rmin_half == 0.4000
    
    lj_params = ff.get_lj_params('CD31C')
    assert lj_params.epsilon == -0.0320
    assert lj_params.rmin_half == 1.8000
    
    lj_params = ff.get_lj_params('ND2A2')
    assert lj_params.epsilon == -0.2000
    assert lj_params.rmin_half == 1.8300
    
    lj_params = ff.get_lj_params('OD2C1A')
    assert lj_params.epsilon == -0.2000
    assert lj_params.rmin_half == 1.7800

def test_nbfix_parameters_from_str(str_file_path):
    """Test reading of NBFIX parameters"""
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # Test specific NBFIX parameters
    epsilon, rmin, found = ff.get_nbfix('ODW', 'CD2O3A')
    assert found == True
    assert epsilon == -0.11528  # Actual output value from analysis script
    assert rmin == 3.4869
    
    epsilon, rmin, found = ff.get_nbfix('ODW', 'ND2A2')
    assert found == True
    assert epsilon == -0.2054
    assert rmin == 3.6369

def test_drude_alpha_thole_parameters(str_file_path):
    """Test reading of Drude ALPHA/THOLE parameters"""
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # Get alpha and thole parameters
    alpha_params = ff.get_alpha_params('ODW')
    assert alpha_params.alpha == -0.97825258
    assert alpha_params.thole == 1.3
    
    alpha_params = ff.get_alpha_params('ND2A2')
    assert alpha_params.alpha == -1.858
    assert alpha_params.thole == 0.126

def test_lonepair_definitions(str_file_path):
    """Test reading of LONEPAIR definitions"""
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # Get lonepair definitions
    lonepairs = ff.get_lonepairs()
    assert len(lonepairs) == 143

def test_anisotropy_definitions(str_file_path):
    """Test reading of ANISOTROPY definitions"""
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # Get anisotropy definitions
    anisotropies = ff.get_anisotropies()
    assert len(anisotropies) == 68

def test_atom_type_count(str_file_path):
    """Test atom type count"""
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # Check if all atom types are read
    # According to analysis results, there should be 175 atom types
    # Use get_num_lj_params to count atom types
    num_atom_types = ff.get_num_lj_params()
    assert num_atom_types >= 175

def test_specific_drude_atom_types(str_file_path):
    """Test specific Drude atom types"""
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # Test Drude particle types
    lj_params = ff.get_lj_params('DRUD')
    assert lj_params.epsilon == 0.0
    assert lj_params.rmin_half == 0.01
    
    # Test water Drude particle
    lj_params = ff.get_lj_params('DOH2')
    assert lj_params.epsilon == 0.0
    assert lj_params.rmin_half == 0.01
    
    # Test lone pair electron
    lj_params = ff.get_lj_params('LPD')
    assert lj_params.epsilon == 0.0
    assert lj_params.rmin_half == 0.01
    
    lj_params = ff.get_lj_params('LPDW')
    assert lj_params.epsilon == 0.0
    assert lj_params.rmin_half == 0.01

if __name__ == "__main__":
    pytest.main([__file__, "-v"])