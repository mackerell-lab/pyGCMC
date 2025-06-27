# tests/io/prmParser/multiple_files_basic.py
"""PRM Parser basic multiple files handling tests."""

import os
import pytest
import pygcmc

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_multiple_nbfix_combinations():
    content = """
NONBONDED nbxmod  5 atom cdiel fshift vatom vdistance vfswitch -
cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5

SOD      0.0       -0.0469    1.41075
POT      0.0       -0.0870    1.76375
CAL      0.0       -0.120     1.367
CLA      0.0       -0.150     2.27
O2L      0.0       -0.120     1.700

NBFIX
SOD    CLA      -0.083875   3.731
POT    CLA      -0.114236   4.081
CAL    CLA      -0.134164   3.727
CAL    O2L      -0.12       3.256
END
"""
    ff = pygcmc.ForceField()
    pygcmc.PRMParser.parse_string(content, ff)

    # First verify LJ parameters are correctly parsed
    sod_params = ff.get_lj_params("SOD")
    assert sod_params.epsilon == pytest.approx(-0.0469)
    assert sod_params.rmin_half == pytest.approx(1.41075)

    cla_params = ff.get_lj_params("CLA")
    assert cla_params.epsilon == pytest.approx(-0.150)
    assert cla_params.rmin_half == pytest.approx(2.27)

    # Then test NBFIX combinations
    epsilon, rmin, found = ff.get_nbfix("SOD", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.083875)
    assert rmin == pytest.approx(3.731)

    # Test reverse order - should still work
    epsilon, rmin, found = ff.get_nbfix("CLA", "SOD")
    assert found == True
    assert epsilon == pytest.approx(-0.083875)
    assert rmin == pytest.approx(3.731)

    epsilon, rmin, found = ff.get_nbfix("POT", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.114236)
    assert rmin == pytest.approx(4.081)

    epsilon, rmin, found = ff.get_nbfix("CAL", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.134164)
    assert rmin == pytest.approx(3.727)

    epsilon, rmin, found = ff.get_nbfix("CAL", "O2L")
    assert found == True
    assert epsilon == pytest.approx(-0.12)
    assert rmin == pytest.approx(3.256)

    # Test non-existent combinations
    epsilon, rmin, found = ff.get_nbfix("SOD", "POT")
    assert found == False


def test_multiple_file_parsing():
    water_ions_file = os.path.join(TEST_DATA_DIR, "toppar_water_ions.str")
    silcs_file = os.path.join(TEST_DATA_DIR, "silcs.str")

    ff = pygcmc.ForceField()
    
    # Parse water_ions file first
    pygcmc.PRMParser.parse_file_to_forcefield(water_ions_file, ff)
    
    # Test parameters from water_ions file
    sod_params = ff.get_lj_params("SOD")
    assert sod_params.epsilon == pytest.approx(-0.0469)
    assert sod_params.rmin_half == pytest.approx(1.41075)

    cla_params = ff.get_lj_params("CLA")
    assert cla_params.epsilon == pytest.approx(-0.150)
    assert cla_params.rmin_half == pytest.approx(2.27)

    epsilon, rmin, found = ff.get_nbfix("SOD", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.0839)
    assert rmin == pytest.approx(3.731)
    
    # Parse silcs file
    pygcmc.PRMParser.parse_file_to_forcefield(silcs_file, ff)

    # Test that original parameters are preserved
    epsilon, rmin, found = ff.get_nbfix("SOD", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.0839)
    assert rmin == pytest.approx(3.731)

    # Test new parameters from silcs file
    lp_params = ff.get_lj_params("LP")
    assert lp_params.epsilon == pytest.approx(0.0)
    assert lp_params.rmin_half == pytest.approx(0.0)

    lq_params = ff.get_lj_params("LQ")
    assert lq_params.epsilon == pytest.approx(0.0)
    assert lq_params.rmin_half == pytest.approx(0.0)

    # Test NBFIX parameters from silcs file
    epsilon, rmin, found = ff.get_nbfix("LP", "LP")
    assert found == True
    assert epsilon == pytest.approx(-0.01)
    assert rmin == pytest.approx(12.0)

    epsilon, rmin, found = ff.get_nbfix("LQ", "LQ")
    assert found == True
    assert epsilon == pytest.approx(-0.01)
    assert rmin == pytest.approx(12.0)

    # Test reverse order access
    epsilon, rmin, found = ff.get_nbfix("LP", "LP")
    assert found == True
    assert epsilon == pytest.approx(-0.01)
    assert rmin == pytest.approx(12.0)


def test_random_parameter_combinations():
    """Test random parameter combinations from multiple parameter files."""
    water_ions_file = os.path.join(TEST_DATA_DIR, "toppar_water_ions.str")
    silcs_file = os.path.join(TEST_DATA_DIR, "silcs.str")

    ff = pygcmc.ForceField()
    
    # Parse both files
    pygcmc.PRMParser.parse_file_to_forcefield(water_ions_file, ff)
    pygcmc.PRMParser.parse_file_to_forcefield(silcs_file, ff)

    # Test nonbonded parameters from water_ions file
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

    # Test random LJ parameters from water_ions file
    # Test water parameters
    ht_params = ff.get_lj_params("HT")
    assert ht_params.epsilon == pytest.approx(-0.046)
    assert ht_params.rmin_half == pytest.approx(0.2245)

    ot_params = ff.get_lj_params("OT")
    assert ot_params.epsilon == pytest.approx(-0.1521)
    assert ot_params.rmin_half == pytest.approx(1.7682)

    # Test ion parameters
    sod_params = ff.get_lj_params("SOD")
    assert sod_params.epsilon == pytest.approx(-0.0469)
    assert sod_params.rmin_half == pytest.approx(1.41075)

    cal_params = ff.get_lj_params("CAL")
    assert cal_params.epsilon == pytest.approx(-0.120)
    assert cal_params.rmin_half == pytest.approx(1.367)

    cla_params = ff.get_lj_params("CLA")
    assert cla_params.epsilon == pytest.approx(-0.150)
    assert cla_params.rmin_half == pytest.approx(2.27)

    # Test random NBFIX parameters from water_ions file
    # Test ion-ion interactions
    epsilon, rmin, found = ff.get_nbfix("SOD", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.0839)
    assert rmin == pytest.approx(3.731)

    epsilon, rmin, found = ff.get_nbfix("CAL", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.134164)
    assert rmin == pytest.approx(3.727)

    epsilon, rmin, found = ff.get_nbfix("CAL", "O2L")
    assert found == True
    assert epsilon == pytest.approx(-0.12)
    assert rmin == pytest.approx(3.256)

    # Test parameters from silcs file
    # Test LJ parameters
    lp_params = ff.get_lj_params("LP")
    assert lp_params.epsilon == pytest.approx(0.0)
    assert lp_params.rmin_half == pytest.approx(0.0)

    # Test NBFIX parameters
    epsilon, rmin, found = ff.get_nbfix("LP", "LP")
    assert found == True
    assert epsilon == pytest.approx(-0.01)
    assert rmin == pytest.approx(12.0)

    epsilon, rmin, found = ff.get_nbfix("LQ", "LQ")
    assert found == True
    assert epsilon == pytest.approx(-0.01)
    assert rmin == pytest.approx(12.0)

    # Test non-existent combinations
    epsilon, rmin, found = ff.get_nbfix("LP", "SOD")
    assert found == False

    epsilon, rmin, found = ff.get_nbfix("LQ", "CLA")
    assert found == False

    # Test parameter overriding
    # Parse water_ions file again to ensure parameters are not duplicated or corrupted
    pygcmc.PRMParser.parse_file_to_forcefield(water_ions_file, ff)
    
    # Verify parameters remain consistent
    sod_params = ff.get_lj_params("SOD")
    assert sod_params.epsilon == pytest.approx(-0.0469)
    assert sod_params.rmin_half == pytest.approx(1.41075)

    epsilon, rmin, found = ff.get_nbfix("LP", "LP")
    assert found == True
    assert epsilon == pytest.approx(-0.01)
    assert rmin == pytest.approx(12.0)


def test_parse_multiple_files():
    """Test parsing multiple parameter files at once."""
    # Define the files to parse
    water_ions_file = os.path.join(TEST_DATA_DIR, "toppar_water_ions.str")
    silcs_file = os.path.join(TEST_DATA_DIR, "silcs.str")
    cgenff_file = os.path.join(TEST_DATA_DIR, "par_all36_cgenff.prm")
    
    # Parse all files at once
    ff = pygcmc.PRMParser.parse_files([water_ions_file, silcs_file, cgenff_file])
    
    # Test parameters from water_ions.str
    # Test water parameters
    ht_params = ff.get_lj_params("HT")
    assert ht_params.epsilon == pytest.approx(-0.046)
    assert ht_params.rmin_half == pytest.approx(0.2245)
    
    # Test ion parameters
    sod_params = ff.get_lj_params("SOD")
    assert sod_params.epsilon == pytest.approx(-0.0469)
    assert sod_params.rmin_half == pytest.approx(1.41075)
    
    # Test NBFIX parameters
    epsilon, rmin, found = ff.get_nbfix("SOD", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.0839)
    assert rmin == pytest.approx(3.731)
    
    # Test parameters from silcs.str
    lp_params = ff.get_lj_params("LP")
    assert lp_params.epsilon == pytest.approx(0.0)
    assert lp_params.rmin_half == pytest.approx(0.0)
    
    # Test parameters from par_all36_cgenff.prm
    # Test some CGenFF specific parameters
    assert ff.get_atom_mass("HGA1") == pytest.approx(1.00800)  # alphatic proton, CH
    assert ff.get_atom_mass("CG2R61") == pytest.approx(12.01100)  # 6-mem aromatic C
    assert ff.get_atom_mass("NG2S1") == pytest.approx(14.00700)  # peptide nitrogen


def test_parse_multiple_files_with_invalid():
    """Test parsing multiple files with an invalid file."""
    water_ions_file = os.path.join(TEST_DATA_DIR, "toppar_water_ions.str")
    
    # Try to parse with a non-existent file
    with pytest.raises(RuntimeError):
        _ = pygcmc.PRMParser.parse_files([water_ions_file, "nonexistent.str"])