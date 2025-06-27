# tests/io/prmParser/basic_parsing.py
"""PRM Parser basic parsing tests."""

import os
import pytest
import pygcmc

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_parse_file_direct():
    """Test the new parse_file method that returns a ForceField object directly."""
    param_file = os.path.join(TEST_DATA_DIR, "toppar_water_ions.str")

    # Test direct parsing
    ff = pygcmc.PRMParser.parse_file(param_file)

    # Test some known values from the file
    params = ff.get_nonbonded_params()
    assert params.nbxmod == 5
    assert params.cutnb == pytest.approx(14.0)

    # Test some LJ parameters
    sod_params = ff.get_lj_params("SOD")
    assert sod_params.epsilon == pytest.approx(-0.0469)
    assert sod_params.rmin_half == pytest.approx(1.41075)

    # Test some NBFIX parameters
    epsilon, rmin, found = ff.get_nbfix("SOD", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.0839)
    assert rmin == pytest.approx(3.731)


def test_parse_nonbonded_from_string():
    content = """
NONBONDED nbxmod 5 cdiel fshift vatom vdistance vfswitch -
cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5

!TIP3P LJ parameters
HT       0.0       -0.046     0.2245
OT       0.0       -0.1521    1.7682
"""
    ff = pygcmc.ForceField()
    pygcmc.PRMParser.parse_string(content, ff)

    # Test nonbonded parameters
    params = ff.get_nonbonded_params()
    assert params.nbxmod == 5
    assert params.cdiel == True
    assert params.fshift == True
    assert params.vatom == True
    assert params.vdistance == True
    assert params.vfswitch == True
    assert params.cutnb == pytest.approx(14.0)
    assert params.ctofnb == pytest.approx(12.0)
    assert params.ctonnb == pytest.approx(10.0)
    assert params.eps == pytest.approx(1.0)
    assert params.e14fac == pytest.approx(1.0)
    assert params.wmin == pytest.approx(1.5)

    # Test LJ parameters
    ht_params = ff.get_lj_params("HT")
    assert ht_params.epsilon == pytest.approx(-0.046)
    assert ht_params.rmin_half == pytest.approx(0.2245)

    ot_params = ff.get_lj_params("OT")
    assert ot_params.epsilon == pytest.approx(-0.1521)
    assert ot_params.rmin_half == pytest.approx(1.7682)


def test_parse_nbfix_from_string():
    content = """
NBFIX
SOD    CLA      -0.083875   3.731
POT    CLA      -0.114236   4.081
"""
    ff = pygcmc.ForceField()
    pygcmc.PRMParser.parse_string(content, ff)

    # Test NBFIX parameters
    epsilon, rmin, found = ff.get_nbfix("SOD", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.083875)
    assert rmin == pytest.approx(3.731)

    # Test symmetry
    epsilon, rmin, found = ff.get_nbfix("CLA", "SOD")
    assert found == True
    assert epsilon == pytest.approx(-0.083875)
    assert rmin == pytest.approx(3.731)

    epsilon, rmin, found = ff.get_nbfix("POT", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.114236)
    assert rmin == pytest.approx(4.081)


def test_parse_from_file():
    """Test the original parse_file method that takes a ForceField object."""
    param_file = os.path.join(TEST_DATA_DIR, "toppar_water_ions.str")

    ff = pygcmc.ForceField()
    pygcmc.PRMParser.parse_file_to_forcefield(param_file, ff)  # Changed from parse_file to parse_file_to_forcefield

    # Test some known values from the file
    params = ff.get_nonbonded_params()
    assert params.nbxmod == 5
    assert params.cutnb == pytest.approx(14.0)

    # Test some LJ parameters
    sod_params = ff.get_lj_params("SOD")
    assert sod_params.epsilon == pytest.approx(-0.0469)
    assert sod_params.rmin_half == pytest.approx(1.41075)

    # Test some NBFIX parameters
    epsilon, rmin, found = ff.get_nbfix("SOD", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.0839)
    assert rmin == pytest.approx(3.731)


def test_parse_comments_and_empty_lines():
    content = """
! This is a comment
   ! This is an indented comment

NONBONDED nbxmod  5 atom cdiel fshift vatom vdistance vfswitch -
cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5

! Comment between sections
SOD      0.0       -0.0469    1.41075  ! inline comment
CLA      0.0       -0.150      2.27    ! another comment
   
NBFIX
! Another comment
SOD    CLA      -0.083875   3.731 ! inline comment
"""
    ff = pygcmc.ForceField()
    pygcmc.PRMParser.parse_string(content, ff)
    
    # Test that comments didn't affect parsing
    epsilon, rmin, found = ff.get_nbfix("SOD", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.083875)
    assert rmin == pytest.approx(3.731)

    # Test that atom parameters were parsed correctly
    sod_params = ff.get_lj_params("SOD")
    assert sod_params.epsilon == pytest.approx(-0.0469)
    assert sod_params.rmin_half == pytest.approx(1.41075)


def test_special_formatting():
    content = """
NONBONDED nbxmod  5 atom cdiel fshift vatom vdistance vfswitch -
cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5

!Values with scientific notation and different spacing
HT          0.0    -4.6e-2     0.2245
OT     0.0         -0.1521        1.7682
"""
    ff = pygcmc.ForceField()
    pygcmc.PRMParser.parse_string(content, ff)

    # Test scientific notation parsing
    ht_params = ff.get_lj_params("HT")
    assert ht_params.epsilon == pytest.approx(-0.046)
    assert ht_params.rmin_half == pytest.approx(0.2245)

    # Test irregular spacing parsing
    ot_params = ff.get_lj_params("OT")
    assert ot_params.epsilon == pytest.approx(-0.1521)
    assert ot_params.rmin_half == pytest.approx(1.7682)


# def test_atom_priority_ordering():
#     """Test if parameter sorting rules conform to CHARMM specifications"""
#     ff = pygcmc.ForceField()
#     
#     # Test bond parameter sorting - using correct CHARMM format
#     content = """
# BONDS
# !V(bond) = Kb(b - b0)**2
# !
# !Kb: kcal/mole/A**2
# !b0: A
# !
# !atom type Kb          b0
# CT1  CS    222.500   1.5380 ! 
# """
#     pygcmc.PRMParser.parse_string(content, ff)
#     assert ff.get_bond_params("CT1", "CS") is not None
#     assert ff.get_bond_params("CS", "CT1") is None  # Should automatically sort when storing
#     
#     # Test dihedral angle inversion rules
#     content = """
# DIHEDRALS
# !V(dihedral) = Kchi(1 + cos(n(chi) - delta))
# !
# !Kchi: kcal/mole
# !n: multiplicity
# !delta: degrees
# !
# !atom types             Kchi    n   delta
# CS   CT1  CT2  HA2    0.200   3     0.00 ! 
# """
#     pygcmc.PRMParser.parse_string(content, ff)
#     assert ff.get_dihedral_params("HA2", "CT2", "CT1", "CS") is not None