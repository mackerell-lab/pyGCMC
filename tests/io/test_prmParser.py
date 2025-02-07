# tests/io/test_prmParser.py

import os
import pytest
import pygcmc

def test_parse_file_direct():
    """Test the new parse_file method that returns a ForceField object directly."""
    test_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(os.path.dirname(test_dir), "data")
    param_file = os.path.join(data_dir, "toppar_water_ions.str")

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
    test_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(os.path.dirname(test_dir), "data")
    param_file = os.path.join(data_dir, "toppar_water_ions.str")

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

def test_invalid_file():
    """Test invalid file handling with new API."""
    with pytest.raises(RuntimeError):
        _ = pygcmc.PRMParser.parse_file("nonexistent.str")

def test_invalid_file_old_api():
    """Test invalid file handling with old API."""
    ff = pygcmc.ForceField()
    with pytest.raises(RuntimeError):
        pygcmc.PRMParser.parse_file_to_forcefield("nonexistent.str", ff)

def test_invalid_atom_type():
    ff = pygcmc.ForceField()
    with pytest.raises(KeyError):
        _ = ff.lj_params["INVALID"]

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

def test_malformed_parameters():
    # Missing value
    with pytest.raises(RuntimeError):
        content = """
NONBONDED nbxmod  5 atom cdiel fshift vatom vdistance vfswitch -
cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5

HT       0.0       -0.046
END
"""
        ff = pygcmc.ForceField()
        pygcmc.PRMParser.parse_string(content, ff)

    # Invalid number format
    with pytest.raises(RuntimeError):
        content = """
NONBONDED nbxmod  5 atom cdiel fshift vatom vdistance vfswitch -
cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5

HT       0.0       -0.046     abc
END
"""
        ff = pygcmc.ForceField()
        pygcmc.PRMParser.parse_string(content, ff)

    # Test that invalid NBFIX lines are skipped with warning
    content = """
NONBONDED nbxmod  5 atom cdiel fshift vatom vdistance vfswitch -
cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5

SOD      0.0       -0.0469    1.41075
CLA      0.0       -0.150     2.27

NBFIX
SOD    CLA    invalid    3.731
END
"""
    ff = pygcmc.ForceField()
    pygcmc.PRMParser.parse_string(content, ff)
    # Verify that SOD and CLA parameters were still parsed correctly
    sod_params = ff.get_lj_params("SOD")
    assert sod_params.epsilon == pytest.approx(-0.0469)
    assert sod_params.rmin_half == pytest.approx(1.41075)

    cla_params = ff.get_lj_params("CLA")
    assert cla_params.epsilon == pytest.approx(-0.150)
    assert cla_params.rmin_half == pytest.approx(2.27)

    # Verify that no NBFIX parameters were added
    epsilon, rmin, found = ff.get_nbfix("SOD", "CLA")
    assert found == False

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

def test_multiple_file_parsing():
    test_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(os.path.dirname(test_dir), "data")
    water_ions_file = os.path.join(data_dir, "toppar_water_ions.str")
    silcs_file = os.path.join(data_dir, "silcs.str")

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
    test_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(os.path.dirname(test_dir), "data")
    water_ions_file = os.path.join(data_dir, "toppar_water_ions.str")
    silcs_file = os.path.join(data_dir, "silcs.str")

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

def test_multiple_parameter_files():
    """Test reading and combining multiple parameter files."""
    test_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(os.path.dirname(test_dir), "data")
    water_ions_file = os.path.join(data_dir, "toppar_water_ions.str")
    silcs_file = os.path.join(data_dir, "silcs.str")
    prot_file = os.path.join(data_dir, "par_all36m_prot.prm")

    ff = pygcmc.ForceField()
    
    # Parse all files
    pygcmc.PRMParser.parse_file_to_forcefield(water_ions_file, ff)
    pygcmc.PRMParser.parse_file_to_forcefield(silcs_file, ff)
    pygcmc.PRMParser.parse_file_to_forcefield(prot_file, ff)

    # 1. Verify nonbonded parameters (from water_ions.str)
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

    # 2. Verify water parameters (from water_ions.str)
    # TIP3P water
    ht_params = ff.get_lj_params("HT")
    assert ht_params.epsilon == pytest.approx(-0.046)
    assert ht_params.rmin_half == pytest.approx(0.2245)

    ot_params = ff.get_lj_params("OT")
    assert ot_params.epsilon == pytest.approx(-0.1521)
    assert ot_params.rmin_half == pytest.approx(1.7682)

    # 3. Verify ion parameters (from water_ions.str)
    # Sodium
    sod_params = ff.get_lj_params("SOD")
    assert sod_params.epsilon == pytest.approx(-0.0469)
    assert sod_params.rmin_half == pytest.approx(1.41075)

    # Calcium
    cal_params = ff.get_lj_params("CAL")
    assert cal_params.epsilon == pytest.approx(-0.120)
    assert cal_params.rmin_half == pytest.approx(1.367)

    # Chloride
    cla_params = ff.get_lj_params("CLA")
    assert cla_params.epsilon == pytest.approx(-0.150)
    assert cla_params.rmin_half == pytest.approx(2.27)

    # 4. Verify NBFIX parameters (from water_ions.str)
    # Ion-ion interactions
    epsilon, rmin, found = ff.get_nbfix("SOD", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.0839)
    assert rmin == pytest.approx(3.731)

    epsilon, rmin, found = ff.get_nbfix("CAL", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.134164)
    assert rmin == pytest.approx(3.727)

    # Ion-oxygen interactions
    epsilon, rmin, found = ff.get_nbfix("CAL", "O2L")
    assert found == True
    assert epsilon == pytest.approx(-0.12)
    assert rmin == pytest.approx(3.256)

    epsilon, rmin, found = ff.get_nbfix("SOD", "OC")
    assert found == True
    assert epsilon == pytest.approx(-0.07502)
    assert rmin == pytest.approx(3.23)

    # 5. Verify SILCS parameters (from silcs.str)
    # LP parameters
    lp_params = ff.get_lj_params("LP")
    assert lp_params.epsilon == pytest.approx(0.0)
    assert lp_params.rmin_half == pytest.approx(0.0)

    # LQ parameters
    lq_params = ff.get_lj_params("LQ")
    assert lq_params.epsilon == pytest.approx(0.0)
    assert lq_params.rmin_half == pytest.approx(0.0)

    # SILCS NBFIX parameters
    epsilon, rmin, found = ff.get_nbfix("LP", "LP")
    assert found == True
    assert epsilon == pytest.approx(-0.01)
    assert rmin == pytest.approx(12.0)

    epsilon, rmin, found = ff.get_nbfix("LQ", "LQ")
    assert found == True
    assert epsilon == pytest.approx(-0.01)
    assert rmin == pytest.approx(12.0)

    # 6. Verify additional NBFIX parameters
    epsilon, rmin, found = ff.get_nbfix("NC2", "OC")
    assert found == True
    assert epsilon == pytest.approx(-0.154919, abs=1e-5)  # From par_all36m_prot.prm
    assert rmin == pytest.approx(3.637)  # Fixed: actual value from par_all36m_prot.prm

    # 7. Verify non-existent combinations
    # Between SILCS and ions
    epsilon, rmin, found = ff.get_nbfix("LP", "SOD")
    assert found == False

    epsilon, rmin, found = ff.get_nbfix("LQ", "CLA")
    assert found == False

    # Between ions
    epsilon, rmin, found = ff.get_nbfix("SOD", "POT")
    assert found == False

    # 8. Verify parameter overriding behavior
    # Parse water_ions file again to ensure parameters are not duplicated or corrupted
    pygcmc.PRMParser.parse_file_to_forcefield(water_ions_file, ff)
    
    # Verify parameters remain consistent
    sod_params = ff.get_lj_params("SOD")
    assert sod_params.epsilon == pytest.approx(-0.0469)
    assert sod_params.rmin_half == pytest.approx(1.41075)

    # SILCS parameters should still be present
    epsilon, rmin, found = ff.get_nbfix("LP", "LP")
    assert found == True
    assert epsilon == pytest.approx(-0.01)
    assert rmin == pytest.approx(12.0)

def test_prm_and_str_files():
    """Test reading both .prm and .str files together."""
    test_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(os.path.dirname(test_dir), "data")
    
    # Load .str files
    water_ions_file = os.path.join(data_dir, "toppar_water_ions.str")
    silcs_file = os.path.join(data_dir, "silcs.str")

    ff = pygcmc.ForceField()
    
    # Parse files in a specific order
    # First load water and ion parameters
    pygcmc.PRMParser.parse_file_to_forcefield(water_ions_file, ff)
    # Then load SILCS parameters
    pygcmc.PRMParser.parse_file_to_forcefield(silcs_file, ff)

    # 1. Verify water parameters (from water_ions.str)
    # TIP3P water
    ht_params = ff.get_lj_params("HT")
    assert ht_params.epsilon == pytest.approx(-0.046)
    assert ht_params.rmin_half == pytest.approx(0.2245)

    ot_params = ff.get_lj_params("OT")
    assert ot_params.epsilon == pytest.approx(-0.1521)
    assert ot_params.rmin_half == pytest.approx(1.7682)

    # 2. Verify ion parameters (from water_ions.str)
    # Sodium
    sod_params = ff.get_lj_params("SOD")
    assert sod_params.epsilon == pytest.approx(-0.0469)
    assert sod_params.rmin_half == pytest.approx(1.41075)

    # Calcium
    cal_params = ff.get_lj_params("CAL")
    assert cal_params.epsilon == pytest.approx(-0.120)
    assert cal_params.rmin_half == pytest.approx(1.367)

    # Chloride
    cla_params = ff.get_lj_params("CLA")
    assert cla_params.epsilon == pytest.approx(-0.150)
    assert cla_params.rmin_half == pytest.approx(2.27)

    # 3. Verify NBFIX parameters from water_ions.str
    # Ion-ion interactions
    epsilon, rmin, found = ff.get_nbfix("SOD", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.0839)
    assert rmin == pytest.approx(3.731)

    epsilon, rmin, found = ff.get_nbfix("CAL", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.134164)
    assert rmin == pytest.approx(3.727)

    # Ion-oxygen interactions
    epsilon, rmin, found = ff.get_nbfix("CAL", "O2L")
    assert found == True
    assert epsilon == pytest.approx(-0.12)
    assert rmin == pytest.approx(3.256)

    epsilon, rmin, found = ff.get_nbfix("SOD", "OC")
    assert found == True
    assert epsilon == pytest.approx(-0.07502)
    assert rmin == pytest.approx(3.23)

    # 4. Verify SILCS parameters (from silcs.str)
    # LP parameters
    lp_params = ff.get_lj_params("LP")
    assert lp_params.epsilon == pytest.approx(0.0)
    assert lp_params.rmin_half == pytest.approx(0.0)

    # LQ parameters
    lq_params = ff.get_lj_params("LQ")
    assert lq_params.epsilon == pytest.approx(0.0)
    assert lq_params.rmin_half == pytest.approx(0.0)

    # SILCS NBFIX parameters
    epsilon, rmin, found = ff.get_nbfix("LP", "LP")
    assert found == True
    assert epsilon == pytest.approx(-0.01)
    assert rmin == pytest.approx(12.0)

    epsilon, rmin, found = ff.get_nbfix("LQ", "LQ")
    assert found == True
    assert epsilon == pytest.approx(-0.01)
    assert rmin == pytest.approx(12.0)

    # 5. Verify non-existent combinations
    # Between SILCS and ions
    epsilon, rmin, found = ff.get_nbfix("LP", "SOD")
    assert found == False

    epsilon, rmin, found = ff.get_nbfix("LQ", "CLA")
    assert found == False

    # Between ions
    epsilon, rmin, found = ff.get_nbfix("SOD", "POT")
    assert found == False

    # 6. Test parameter overriding and coexistence
    # Parse water_ions file again to ensure parameters are not duplicated or corrupted
    pygcmc.PRMParser.parse_file_to_forcefield(water_ions_file, ff)
    
    # Water parameters should remain unchanged
    ht_params = ff.get_lj_params("HT")
    assert ht_params.epsilon == pytest.approx(-0.046)
    assert ht_params.rmin_half == pytest.approx(0.2245)
    
    # SILCS parameters should still be present
    epsilon, rmin, found = ff.get_nbfix("LP", "LP")
    assert found == True
    assert epsilon == pytest.approx(-0.01)
    assert rmin == pytest.approx(12.0)

    # 7. Verify nonbonded parameters are properly maintained
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

def test_charmm_prm_files():
    """Test parsing of CHARMM force field files (par_all36m_prot.prm)."""
    test_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(os.path.dirname(test_dir), "data")
    prm_file = os.path.join(data_dir, "par_all36m_prot.prm")

    ff = pygcmc.ForceField()
    pygcmc.PRMParser.parse_file_to_forcefield(prm_file, ff)

    # Test header section (ATOMS/MASS)
    # Test some hydrogen atoms
    assert ff.get_atom_mass("H") == pytest.approx(1.00800)    # polar H
    assert ff.get_atom_mass("HC") == pytest.approx(1.00800)   # N-ter H
    assert ff.get_atom_mass("HA") == pytest.approx(1.00800)   # nonpolar H
    assert ff.get_atom_mass("HP") == pytest.approx(1.00800)   # aromatic H
    assert ff.get_atom_mass("HB1") == pytest.approx(1.00800)  # backbone H
    assert ff.get_atom_mass("HB2") == pytest.approx(1.00800)  # aliphatic backbone H

    # Test some carbon atoms
    assert ff.get_atom_mass("C") == pytest.approx(12.01100)    # carbonyl C, peptide backbone
    assert ff.get_atom_mass("CA") == pytest.approx(12.01100)   # aromatic C
    assert ff.get_atom_mass("CT1") == pytest.approx(12.01100)  # aliphatic sp3 C for CH
    assert ff.get_atom_mass("CT2") == pytest.approx(12.01100)  # aliphatic sp3 C for CH2
    assert ff.get_atom_mass("CT3") == pytest.approx(12.01100)  # aliphatic sp3 C for CH3
    assert ff.get_atom_mass("CPH1") == pytest.approx(12.01100) # his CG and CD2 carbons

    # Test some nitrogen atoms
    assert ff.get_atom_mass("N") == pytest.approx(14.00700)    # proline N
    assert ff.get_atom_mass("NH1") == pytest.approx(14.00700)  # peptide nitrogen
    assert ff.get_atom_mass("NH2") == pytest.approx(14.00700)  # amide nitrogen
    assert ff.get_atom_mass("NH3") == pytest.approx(14.00700)  # ammonium nitrogen
    assert ff.get_atom_mass("NR1") == pytest.approx(14.00700)  # neutral his protonated ring nitrogen

    # Test some oxygen atoms
    assert ff.get_atom_mass("O") == pytest.approx(15.99940)    # carbonyl oxygen
    assert ff.get_atom_mass("OB") == pytest.approx(15.99940)   # carbonyl oxygen in acetic acid
    assert ff.get_atom_mass("OC") == pytest.approx(15.99940)   # carboxylate oxygen
    assert ff.get_atom_mass("OH1") == pytest.approx(15.99940)  # hydroxyl oxygen

    # Test sulfur atoms
    assert ff.get_atom_mass("S") == pytest.approx(32.06000)    # sulphur
    assert ff.get_atom_mass("SM") == pytest.approx(32.06000)   # sulfur C-S-S-C type
    assert ff.get_atom_mass("SS") == pytest.approx(32.06000)   # thiolate sulfur

    # Test BONDS section
    # Test some peptide backbone bonds
    key = pygcmc.ForceField.makeTypePair("N", "C")
    print("\nTesting bond N-C")
    print(f"Generated key: {key}")
    print(f"Available keys in bond_params: {list(ff.bond_params.keys())}")
    bond_params = ff.get_bond_params("N", "C")
    assert bond_params.kb == pytest.approx(260.000)
    assert bond_params.b0 == pytest.approx(1.3000)

    bond_params = ff.get_bond_params("C", "O")
    assert bond_params.kb == pytest.approx(620.000)
    assert bond_params.b0 == pytest.approx(1.2300)

    # Test some side chain bonds
    bond_params = ff.get_bond_params("CA", "CA")
    assert bond_params.kb == pytest.approx(305.000)
    assert bond_params.b0 == pytest.approx(1.3750)

    bond_params = ff.get_bond_params("CT2", "OH1")
    assert bond_params.kb == pytest.approx(428.000)
    assert bond_params.b0 == pytest.approx(1.4200)

    # Test some hydrogen bonds
    bond_params = ff.get_bond_params("CT3", "HA3")
    assert bond_params.kb == pytest.approx(322.000)
    assert bond_params.b0 == pytest.approx(1.1110)

    # Test ANGLES section
    # Test some peptide backbone angles
    angle_params = ff.get_angle_params("N", "C", "CT1")
    assert angle_params.ktheta == pytest.approx(20.000)  # N-C-CT1 angle force constant
    assert angle_params.theta0 == pytest.approx(112.5000)  # N-C-CT1 equilibrium angle

    angle_params = ff.get_angle_params("N", "C", "O")
    assert angle_params.ktheta == pytest.approx(80.000)
    assert angle_params.theta0 == pytest.approx(122.5000)

    # Test some side chain angles
    angle_params = ff.get_angle_params("CA", "CA", "CA")
    assert angle_params.ktheta == pytest.approx(40.000)
    assert angle_params.theta0 == pytest.approx(120.0000)

    # Test DIHEDRALS section
    # Test some peptide backbone dihedrals
    # Test 1: Neutral N-terminus dihedral
    dihedral_params = ff.get_dihedral_params("NH2", "CT1", "C", "O")
    assert len(dihedral_params) == 1
    assert dihedral_params[0].kchi == pytest.approx(0.0000)
    assert dihedral_params[0].n == 1
    assert dihedral_params[0].delta == pytest.approx(0.00)

    # Test 2: Proline ring dihedral
    dihedral_params = ff.get_dihedral_params("CT1", "C", "N", "CP1")
    assert len(dihedral_params) == 2  # This dihedral has two terms
    assert dihedral_params[0].kchi == pytest.approx(2.7500)
    assert dihedral_params[0].n == 2  # Changed from 4 to 2 to match parameter file
    assert dihedral_params[1].kchi == pytest.approx(0.3000)
    assert dihedral_params[1].n == 4
    assert dihedral_params[1].delta == pytest.approx(0.00)

    # Test 3: Histidine ring dihedral
    dihedral_params = ff.get_dihedral_params("CPH2", "NR1", "CPH1", "CPH1")
    assert len(dihedral_params) == 1
    assert dihedral_params[0].kchi == pytest.approx(14.0000)
    assert dihedral_params[0].n == 2
    assert dihedral_params[0].delta == pytest.approx(180.00)

    # Test 4: Generic dihedral with wildcard (X)
    dihedral_params = ff.get_dihedral_params("X", "CP1", "C", "X")
    assert len(dihedral_params) == 1
    assert dihedral_params[0].kchi == pytest.approx(0.0000)
    assert dihedral_params[0].n == 6
    assert dihedral_params[0].delta == pytest.approx(180.00)

    # Test NONBONDED parameters
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

    # Test some LJ parameters
    # Test some hydrogen LJ parameters
    h_params = ff.get_lj_params("H")
    assert h_params.epsilon == pytest.approx(-0.0460)
    assert h_params.rmin_half == pytest.approx(0.2245)

    # Test some carbon LJ parameters
    c_params = ff.get_lj_params("C")
    assert c_params.epsilon == pytest.approx(-0.1100)
    assert c_params.rmin_half == pytest.approx(2.0000)

    # Test some nitrogen LJ parameters
    n_params = ff.get_lj_params("N")
    assert n_params.epsilon == pytest.approx(-0.2000)
    assert n_params.rmin_half == pytest.approx(1.8500)

    # Test some oxygen LJ parameters
    o_params = ff.get_lj_params("O")
    assert o_params.epsilon == pytest.approx(-0.1200)
    assert o_params.rmin_half == pytest.approx(1.7000)

    # Test NBFIX parameters if present
    epsilon, rmin, found = ff.get_nbfix("SOD", "CLA")
    assert found == False  # Should be false if no NBFIX in the file

def test_cgenff_prm_file():
    """Test parsing of CHARMM General Force Field file (par_all36_cgenff.prm)."""
    test_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(os.path.dirname(test_dir), "data")
    prm_file = os.path.join(data_dir, "par_all36_cgenff.prm")

    ff = pygcmc.ForceField()
    pygcmc.PRMParser.parse_file_to_forcefield(prm_file, ff)

    # Test header section (ATOMS/MASS)
    # Test some hydrogen atoms
    assert ff.get_atom_mass("HGA1") == pytest.approx(1.00800)  # alphatic proton, CH
    assert ff.get_atom_mass("HGA2") == pytest.approx(1.00800)  # alphatic proton, CH2
    assert ff.get_atom_mass("HGA3") == pytest.approx(1.00800)  # alphatic proton, CH3
    assert ff.get_atom_mass("HGR61") == pytest.approx(1.00800) # aromatic H
    assert ff.get_atom_mass("HGP1") == pytest.approx(1.00800)  # polar H

    # Test some carbon atoms
    assert ff.get_atom_mass("CG2R61") == pytest.approx(12.01100)  # 6-mem aromatic C
    assert ff.get_atom_mass("CG2O1") == pytest.approx(12.01100)   # carbonyl C: amides
    assert ff.get_atom_mass("CG321") == pytest.approx(12.01100)   # aliphatic C for CH2
    assert ff.get_atom_mass("CG331") == pytest.approx(12.01100)   # aliphatic C for methyl group

    # Test some nitrogen atoms
    assert ff.get_atom_mass("NG2S1") == pytest.approx(14.00700)  # peptide nitrogen
    assert ff.get_atom_mass("NG2S2") == pytest.approx(14.00700)  # terminal amide nitrogen
    assert ff.get_atom_mass("NG301") == pytest.approx(14.00700)  # neutral trimethylamine nitrogen

    # Test some oxygen atoms
    assert ff.get_atom_mass("OG2D1") == pytest.approx(15.99940)  # carbonyl O: amides
    assert ff.get_atom_mass("OG2D2") == pytest.approx(15.99940)  # carbonyl O: negative groups
    assert ff.get_atom_mass("OG2D3") == pytest.approx(15.99940)  # carbonyl O: ketones

    # Test BONDS section
    # Test some typical CGenFF bonds
    bond_params = ff.get_bond_params("CG2R61", "CG2R61")
    assert bond_params.kb == pytest.approx(305.000)
    assert bond_params.b0 == pytest.approx(1.3750)

    bond_params = ff.get_bond_params("CG2R61", "HGR61")
    assert bond_params.kb == pytest.approx(340.000)
    assert bond_params.b0 == pytest.approx(1.0800)

    bond_params = ff.get_bond_params("CG2O1", "OG2D1")
    assert bond_params.kb == pytest.approx(620.000)
    assert bond_params.b0 == pytest.approx(1.2300)

    # Test ANGLES section
    # Test some typical CGenFF angles
    angle_params = ff.get_angle_params("CG2R61", "CG2R61", "CG2R61")
    assert angle_params.ktheta == pytest.approx(40.000)
    assert angle_params.theta0 == pytest.approx(120.0000)

    angle_params = ff.get_angle_params("HGR61", "CG2R61", "CG2R61")
    assert angle_params.ktheta == pytest.approx(30.000)
    assert angle_params.theta0 == pytest.approx(120.0000)

    # Test DIHEDRALS section
    # Test some typical CGenFF dihedrals
    dihedral_params = ff.get_dihedral_params("CG2R61", "CG2R61", "CG2R61", "CG2R61")
    assert len(dihedral_params) > 0
    assert dihedral_params[0].kchi == pytest.approx(3.1000)
    assert dihedral_params[0].n == 2
    assert dihedral_params[0].delta == pytest.approx(180.00)

    # Test NONBONDED parameters
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

    # Test some LJ parameters
    # Test some typical CGenFF LJ parameters
    cg2r61_params = ff.get_lj_params("CG2R61")
    assert cg2r61_params.epsilon == pytest.approx(-0.0700)
    assert cg2r61_params.rmin_half == pytest.approx(1.9924)

    hgr61_params = ff.get_lj_params("HGR61")
    assert hgr61_params.epsilon == pytest.approx(-0.0300)
    assert hgr61_params.rmin_half == pytest.approx(1.3582)

    og2d1_params = ff.get_lj_params("OG2D1")
    assert og2d1_params.epsilon == pytest.approx(-0.1200)
    assert og2d1_params.rmin_half == pytest.approx(1.7000)

    # Test NBFIX parameters if present
    epsilon, rmin, found = ff.get_nbfix("CG2R61", "OG2D1")
    assert found == False  # Should be false if no NBFIX in the file

def test_ion_ligand_nbfix():
    """Test ion-ligand NBFIX parameters from toppar_water_ions.str."""
    test_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(os.path.dirname(test_dir), "data")
    param_file = os.path.join(data_dir, "toppar_water_ions.str")

    ff = pygcmc.ForceField()
    pygcmc.PRMParser.parse_file_to_forcefield(param_file, ff)

    # Test SOD-OC interaction (sodium-carboxylate)
    epsilon, rmin, found = ff.get_nbfix("SOD", "OC")
    assert found == True
    assert epsilon == pytest.approx(-0.07502)
    assert rmin == pytest.approx(3.23)

    # Test SOD-O2L interaction (sodium-phosphate)
    epsilon, rmin, found = ff.get_nbfix("SOD", "O2L")
    assert found == True
    assert epsilon == pytest.approx(-0.07502)
    assert rmin == pytest.approx(3.16)  # Fixed: actual value from toppar_water_ions.str

    # Test CAL-OC interaction (calcium-carboxylate)
    epsilon, rmin, found = ff.get_nbfix("CAL", "OC")
    assert found == True
    assert epsilon == pytest.approx(-0.12)
    assert rmin == pytest.approx(3.232)  # Fixed: actual value from toppar_water_ions.str

    # Test special water parameters
    hper_params = ff.get_lj_params("HPER")
    assert hper_params.epsilon == pytest.approx(-0.046)
    assert hper_params.rmin_half == pytest.approx(0.2245)

    oper_params = ff.get_lj_params("OPER")
    assert oper_params.epsilon == pytest.approx(-0.20384)
    assert oper_params.rmin_half == pytest.approx(1.67423)

def test_heterocyclic_parameters():
    """Test parameters for heterocyclic compounds from par_all36_cgenff.prm."""
    test_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(os.path.dirname(test_dir), "data")
    prm_file = os.path.join(data_dir, "par_all36_cgenff.prm")

    ff = pygcmc.ForceField()
    pygcmc.PRMParser.parse_file_to_forcefield(prm_file, ff)

    # Test pyridine-like parameters
    angle_params = ff.get_angle_params("NG2R60", "CG2R64", "NG2S1")
    assert angle_params.ktheta == pytest.approx(40.00)
    assert angle_params.theta0 == pytest.approx(120.00)

    # Test pyrazole-like parameters
    dihedral_params = ff.get_dihedral_params("CG2R51", "CG2R51", "CG2R52", "NG2R50")
    assert len(dihedral_params) > 0
    assert dihedral_params[0].kchi == pytest.approx(8.5000)
    assert dihedral_params[0].n == 2
    assert dihedral_params[0].delta == pytest.approx(180.00)

    # Test halogen interactions
    epsilon, rmin, found = ff.get_nbfix("CLGR1", "OG2D2")
    assert found == True
    assert epsilon == pytest.approx(-2.50)
    assert rmin == pytest.approx(2.8)

    epsilon, rmin, found = ff.get_nbfix("CLGR1", "NG2R51")
    assert found == True
    assert epsilon == pytest.approx(-0.48)
    assert rmin == pytest.approx(3.75)  # Fixed: actual value from par_all36_cgenff.prm

def test_nucleic_parameters():
    """Test nucleic acid related parameters from par_all36_cgenff.prm."""
    test_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(os.path.dirname(test_dir), "data")
    prm_file = os.path.join(data_dir, "par_all36_cgenff.prm")

    ff = pygcmc.ForceField()
    pygcmc.PRMParser.parse_file_to_forcefield(prm_file, ff)

    # Test nucleic acid - halogen interactions
    epsilon, rmin, found = ff.get_nbfix("NN2G", "BRGR1")
    assert found == True
    assert epsilon == pytest.approx(-0.72)
    assert rmin == pytest.approx(3.8)

    epsilon, rmin, found = ff.get_nbfix("ON1C", "CLGR1")
    assert found == True
    assert epsilon == pytest.approx(-0.20)
    assert rmin == pytest.approx(3.4)  # Fixed: actual value from par_all36_cgenff.prm

def test_cross_forcefield_compatibility():
    """Test parameter compatibility between protein and CGenFF force fields."""
    test_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(os.path.dirname(test_dir), "data")
    prot_file = os.path.join(data_dir, "par_all36m_prot.prm")
    cgenff_file = os.path.join(data_dir, "par_all36_cgenff.prm")

    ff = pygcmc.ForceField()
    
    # Load both force fields
    pygcmc.PRMParser.parse_file_to_forcefield(prot_file, ff)
    pygcmc.PRMParser.parse_file_to_forcefield(cgenff_file, ff)

    # Test CT2A parameters (special carbon type in GLU/HSP)
    ct2a_params = ff.get_lj_params("CT2A")
    assert ct2a_params.epsilon == pytest.approx(-0.0560)
    assert ct2a_params.rmin_half == pytest.approx(2.010)

    # Test compatibility of common atom types between force fields
    # Aromatic carbon parameters should be consistent
    ca_prot = ff.get_lj_params("CA")
    cg2r61 = ff.get_lj_params("CG2R61")
    assert abs(ca_prot.epsilon - cg2r61.epsilon) < 0.01  # Should be similar
    assert abs(ca_prot.rmin_half - cg2r61.rmin_half) < 0.1  # Should be similar

    # Test compatibility of peptide backbone parameters
    c_prot = ff.get_lj_params("C")    # Protein carbonyl carbon
    cg2o1 = ff.get_lj_params("CG2O1") # CGenFF carbonyl carbon
    assert abs(c_prot.epsilon - cg2o1.epsilon) < 0.01  # Should be similar
    assert abs(c_prot.rmin_half - cg2o1.rmin_half) < 0.1  # Should be similar

def test_parse_multiple_files():
    """Test parsing multiple parameter files at once."""
    test_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(os.path.dirname(test_dir), "data")
    
    # Define the files to parse
    water_ions_file = os.path.join(data_dir, "toppar_water_ions.str")
    silcs_file = os.path.join(data_dir, "silcs.str")
    cgenff_file = os.path.join(data_dir, "par_all36_cgenff.prm")
    
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
    test_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(os.path.dirname(test_dir), "data")
    water_ions_file = os.path.join(data_dir, "toppar_water_ions.str")
    
    # Try to parse with a non-existent file
    with pytest.raises(RuntimeError):
        _ = pygcmc.PRMParser.parse_files([water_ions_file, "nonexistent.str"])

# def test_atom_priority_ordering():
#     """测试参数排序规则是否符合CHARMM规范"""
#     ff = pygcmc.ForceField()
    
#     # 测试键参数排序 - 使用正确的 CHARMM 格式
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
#     assert ff.get_bond_params("CS", "CT1") is None  # 应自动排序存储
    
#     # 测试二面角反转规则
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
