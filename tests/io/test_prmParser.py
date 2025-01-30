# tests/io/test_prmParser.py

import os
import pytest
import pygcmc

def test_parse_nonbonded_from_string():
    content = """
NONBONDED nbxmod 5 cdiel fshift vatom vdistance vfswitch -
cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5

!TIP3P LJ parameters
HT       0.0       -0.046     0.2245
OT       0.0       -0.1521    1.7682
"""
    ff = pygcmc.ForceField()
    pygcmc.PrmParser.parse_string(content, ff)

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
    assert ht_params.rmin == pytest.approx(0.2245)

    ot_params = ff.get_lj_params("OT")
    assert ot_params.epsilon == pytest.approx(-0.1521)
    assert ot_params.rmin == pytest.approx(1.7682)

def test_parse_nbfix_from_string():
    content = """
NBFIX
SOD    CLA      -0.083875   3.731
POT    CLA      -0.114236   4.081
"""
    ff = pygcmc.ForceField()
    pygcmc.PrmParser.parse_string(content, ff)

    # Test NBFIX parameters
    epsilon, found = ff.get_nbfix("SOD", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.083875)

    # Test symmetry
    epsilon, found = ff.get_nbfix("CLA", "SOD")
    assert found == True
    assert epsilon == pytest.approx(-0.083875)

    epsilon, found = ff.get_nbfix("POT", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.114236)

def test_parse_from_file():
    test_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(os.path.dirname(test_dir), "data")
    param_file = os.path.join(data_dir, "toppar_water_ions.str")

    ff = pygcmc.ForceField()
    pygcmc.PrmParser.parse_file(param_file, ff)

    # Test some known values from the file
    params = ff.get_nonbonded_params()
    assert params.nbxmod == 5
    assert params.cutnb == pytest.approx(14.0)

    # Test some LJ parameters
    sod_params = ff.get_lj_params("SOD")
    assert sod_params.epsilon == pytest.approx(-0.0469)
    assert sod_params.rmin == pytest.approx(1.41075)

    # Test some NBFIX parameters
    epsilon, found = ff.get_nbfix("SOD", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.0839)

def test_invalid_file():
    ff = pygcmc.ForceField()
    with pytest.raises(RuntimeError):
        pygcmc.PrmParser.parse_file("nonexistent.str", ff)

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
    pygcmc.PrmParser.parse_string(content, ff)
    
    # Test that comments didn't affect parsing
    epsilon, found = ff.get_nbfix("SOD", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.083875)

    # Test that atom parameters were parsed correctly
    sod_params = ff.get_lj_params("SOD")
    assert sod_params.epsilon == pytest.approx(-0.0469)
    assert sod_params.rmin == pytest.approx(1.41075)

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
    pygcmc.PrmParser.parse_string(content, ff)

    # First verify LJ parameters are correctly parsed
    sod_params = ff.get_lj_params("SOD")
    assert sod_params.epsilon == pytest.approx(-0.0469)
    assert sod_params.rmin == pytest.approx(1.41075)

    cla_params = ff.get_lj_params("CLA")
    assert cla_params.epsilon == pytest.approx(-0.150)
    assert cla_params.rmin == pytest.approx(2.27)

    # Then test NBFIX combinations
    epsilon, found = ff.get_nbfix("SOD", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.083875)

    # Test reverse order - should still work
    epsilon, found = ff.get_nbfix("CLA", "SOD")
    assert found == True
    assert epsilon == pytest.approx(-0.083875)

    epsilon, found = ff.get_nbfix("POT", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.114236)

    epsilon, found = ff.get_nbfix("CAL", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.134164)

    epsilon, found = ff.get_nbfix("CAL", "O2L")
    assert found == True
    assert epsilon == pytest.approx(-0.12)

    # Test non-existent combinations
    epsilon, found = ff.get_nbfix("SOD", "POT")
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
        pygcmc.PrmParser.parse_string(content, ff)

    # Invalid number format
    with pytest.raises(RuntimeError):
        content = """
NONBONDED nbxmod  5 atom cdiel fshift vatom vdistance vfswitch -
cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5

HT       0.0       -0.046     abc
END
"""
        ff = pygcmc.ForceField()
        pygcmc.PrmParser.parse_string(content, ff)

    # Invalid NBFIX format
    with pytest.raises(RuntimeError):
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
        pygcmc.PrmParser.parse_string(content, ff)

def test_special_formatting():
    content = """
NONBONDED nbxmod  5 atom cdiel fshift vatom vdistance vfswitch -
cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5

!Values with scientific notation and different spacing
HT          0.0    -4.6e-2     0.2245
OT     0.0         -0.1521        1.7682
"""
    ff = pygcmc.ForceField()
    pygcmc.PrmParser.parse_string(content, ff)

    # Test scientific notation parsing
    ht_params = ff.get_lj_params("HT")
    assert ht_params.epsilon == pytest.approx(-0.046)
    assert ht_params.rmin == pytest.approx(0.2245)

    # Test irregular spacing parsing
    ot_params = ff.get_lj_params("OT")
    assert ot_params.epsilon == pytest.approx(-0.1521)
    assert ot_params.rmin == pytest.approx(1.7682)

def test_multiple_file_parsing():
    test_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(os.path.dirname(test_dir), "data")
    water_ions_file = os.path.join(data_dir, "toppar_water_ions.str")
    silcs_file = os.path.join(data_dir, "silcs.str")

    ff = pygcmc.ForceField()
    
    # Parse water_ions file first
    pygcmc.PrmParser.parse_file(water_ions_file, ff)
    
    # Test parameters from water_ions file
    sod_params = ff.get_lj_params("SOD")
    assert sod_params.epsilon == pytest.approx(-0.0469)
    assert sod_params.rmin == pytest.approx(1.41075)

    cla_params = ff.get_lj_params("CLA")
    assert cla_params.epsilon == pytest.approx(-0.150)
    assert cla_params.rmin == pytest.approx(2.27)

    epsilon, found = ff.get_nbfix("SOD", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.0839)
    
    # Parse silcs file
    pygcmc.PrmParser.parse_file(silcs_file, ff)

    # Test that original parameters are preserved
    epsilon, found = ff.get_nbfix("SOD", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.0839)

    # Test new parameters from silcs file
    lp_params = ff.get_lj_params("LP")
    assert lp_params.epsilon == pytest.approx(0.0)
    assert lp_params.rmin == pytest.approx(0.0)

    lq_params = ff.get_lj_params("LQ")
    assert lq_params.epsilon == pytest.approx(0.0)
    assert lq_params.rmin == pytest.approx(0.0)

    # Test NBFIX parameters from silcs file
    epsilon, found = ff.get_nbfix("LP", "LP")
    assert found == True
    assert epsilon == pytest.approx(-0.01)

    epsilon, found = ff.get_nbfix("LQ", "LQ")
    assert found == True
    assert epsilon == pytest.approx(-0.01)

    # Test reverse order access
    epsilon, found = ff.get_nbfix("LP", "LP")
    assert found == True
    assert epsilon == pytest.approx(-0.01)

def test_random_parameter_combinations():
    """Test random parameter combinations from multiple parameter files."""
    test_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(os.path.dirname(test_dir), "data")
    water_ions_file = os.path.join(data_dir, "toppar_water_ions.str")
    silcs_file = os.path.join(data_dir, "silcs.str")

    ff = pygcmc.ForceField()
    
    # Parse both files
    pygcmc.PrmParser.parse_file(water_ions_file, ff)
    pygcmc.PrmParser.parse_file(silcs_file, ff)

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
    assert ht_params.rmin == pytest.approx(0.2245)

    ot_params = ff.get_lj_params("OT")
    assert ot_params.epsilon == pytest.approx(-0.1521)
    assert ot_params.rmin == pytest.approx(1.7682)

    # Test ion parameters
    sod_params = ff.get_lj_params("SOD")
    assert sod_params.epsilon == pytest.approx(-0.0469)
    assert sod_params.rmin == pytest.approx(1.41075)

    cal_params = ff.get_lj_params("CAL")
    assert cal_params.epsilon == pytest.approx(-0.120)
    assert cal_params.rmin == pytest.approx(1.367)

    cla_params = ff.get_lj_params("CLA")
    assert cla_params.epsilon == pytest.approx(-0.150)
    assert cla_params.rmin == pytest.approx(2.27)

    # Test random NBFIX parameters from water_ions file
    # Test ion-ion interactions
    epsilon, found = ff.get_nbfix("SOD", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.0839)

    epsilon, found = ff.get_nbfix("CAL", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.134164)

    epsilon, found = ff.get_nbfix("CAL", "O2L")
    assert found == True
    assert epsilon == pytest.approx(-0.12)

    # Test parameters from silcs file
    # Test LJ parameters
    lp_params = ff.get_lj_params("LP")
    assert lp_params.epsilon == pytest.approx(0.0)
    assert lp_params.rmin == pytest.approx(0.0)

    # Test NBFIX parameters
    epsilon, found = ff.get_nbfix("LP", "LP")
    assert found == True
    assert epsilon == pytest.approx(-0.01)

    # Test non-existent combinations
    epsilon, found = ff.get_nbfix("LP", "SOD")
    assert found == False

    epsilon, found = ff.get_nbfix("LQ", "CLA")
    assert found == False

    # Test parameter overriding
    # Parse water_ions file again to ensure parameters are not duplicated or corrupted
    pygcmc.PrmParser.parse_file(water_ions_file, ff)
    
    # Verify parameters remain consistent
    sod_params = ff.get_lj_params("SOD")
    assert sod_params.epsilon == pytest.approx(-0.0469)
    assert sod_params.rmin == pytest.approx(1.41075)

    epsilon, found = ff.get_nbfix("LP", "LP")
    assert found == True
    assert epsilon == pytest.approx(-0.01)

def test_multiple_parameter_files():
    """Test reading and combining multiple parameter files."""
    test_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(os.path.dirname(test_dir), "data")
    water_ions_file = os.path.join(data_dir, "toppar_water_ions.str")
    silcs_file = os.path.join(data_dir, "silcs.str")

    ff = pygcmc.ForceField()
    
    # Parse all files
    pygcmc.PrmParser.parse_file(water_ions_file, ff)
    pygcmc.PrmParser.parse_file(silcs_file, ff)

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
    assert ht_params.rmin == pytest.approx(0.2245)

    ot_params = ff.get_lj_params("OT")
    assert ot_params.epsilon == pytest.approx(-0.1521)
    assert ot_params.rmin == pytest.approx(1.7682)

    # 3. Verify ion parameters (from water_ions.str)
    # Sodium
    sod_params = ff.get_lj_params("SOD")
    assert sod_params.epsilon == pytest.approx(-0.0469)
    assert sod_params.rmin == pytest.approx(1.41075)

    # Calcium
    cal_params = ff.get_lj_params("CAL")
    assert cal_params.epsilon == pytest.approx(-0.120)
    assert cal_params.rmin == pytest.approx(1.367)

    # Chloride
    cla_params = ff.get_lj_params("CLA")
    assert cla_params.epsilon == pytest.approx(-0.150)
    assert cla_params.rmin == pytest.approx(2.27)

    # 4. Verify NBFIX parameters (from water_ions.str)
    # Ion-ion interactions
    epsilon, found = ff.get_nbfix("SOD", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.0839)

    epsilon, found = ff.get_nbfix("CAL", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.134164)

    # Ion-oxygen interactions
    epsilon, found = ff.get_nbfix("CAL", "O2L")
    assert found == True
    assert epsilon == pytest.approx(-0.12)

    epsilon, found = ff.get_nbfix("SOD", "OC")
    assert found == True
    assert epsilon == pytest.approx(-0.07502)

    # 5. Verify SILCS parameters (from silcs.str)
    # LP parameters
    lp_params = ff.get_lj_params("LP")
    assert lp_params.epsilon == pytest.approx(0.0)
    assert lp_params.rmin == pytest.approx(0.0)

    # LQ parameters
    lq_params = ff.get_lj_params("LQ")
    assert lq_params.epsilon == pytest.approx(0.0)
    assert lq_params.rmin == pytest.approx(0.0)

    # SILCS NBFIX parameters
    epsilon, found = ff.get_nbfix("LP", "LP")
    assert found == True
    assert epsilon == pytest.approx(-0.01)

    epsilon, found = ff.get_nbfix("LQ", "LQ")
    assert found == True
    assert epsilon == pytest.approx(-0.01)

    # 6. Verify non-existent combinations
    # Between SILCS and ions
    epsilon, found = ff.get_nbfix("LP", "SOD")
    assert found == False

    epsilon, found = ff.get_nbfix("LQ", "CLA")
    assert found == False

    # Between ions
    epsilon, found = ff.get_nbfix("SOD", "POT")
    assert found == False

    # 7. Verify parameter overriding behavior
    # Parse water_ions file again to ensure parameters are not duplicated or corrupted
    pygcmc.PrmParser.parse_file(water_ions_file, ff)
    
    # Verify parameters remain consistent
    sod_params = ff.get_lj_params("SOD")
    assert sod_params.epsilon == pytest.approx(-0.0469)
    assert sod_params.rmin == pytest.approx(1.41075)

    # SILCS parameters should still be present
    epsilon, found = ff.get_nbfix("LP", "LP")
    assert found == True
    assert epsilon == pytest.approx(-0.01)

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
    pygcmc.PrmParser.parse_file(water_ions_file, ff)
    # Then load SILCS parameters
    pygcmc.PrmParser.parse_file(silcs_file, ff)

    # 1. Verify water parameters (from water_ions.str)
    # TIP3P water
    ht_params = ff.get_lj_params("HT")
    assert ht_params.epsilon == pytest.approx(-0.046)
    assert ht_params.rmin == pytest.approx(0.2245)

    ot_params = ff.get_lj_params("OT")
    assert ot_params.epsilon == pytest.approx(-0.1521)
    assert ot_params.rmin == pytest.approx(1.7682)

    # 2. Verify ion parameters (from water_ions.str)
    # Sodium
    sod_params = ff.get_lj_params("SOD")
    assert sod_params.epsilon == pytest.approx(-0.0469)
    assert sod_params.rmin == pytest.approx(1.41075)

    # Calcium
    cal_params = ff.get_lj_params("CAL")
    assert cal_params.epsilon == pytest.approx(-0.120)
    assert cal_params.rmin == pytest.approx(1.367)

    # Chloride
    cla_params = ff.get_lj_params("CLA")
    assert cla_params.epsilon == pytest.approx(-0.150)
    assert cla_params.rmin == pytest.approx(2.27)

    # 3. Verify NBFIX parameters from water_ions.str
    # Ion-ion interactions
    epsilon, found = ff.get_nbfix("SOD", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.0839)

    epsilon, found = ff.get_nbfix("CAL", "CLA")
    assert found == True
    assert epsilon == pytest.approx(-0.134164)

    # Ion-oxygen interactions
    epsilon, found = ff.get_nbfix("CAL", "O2L")
    assert found == True
    assert epsilon == pytest.approx(-0.12)

    epsilon, found = ff.get_nbfix("SOD", "OC")
    assert found == True
    assert epsilon == pytest.approx(-0.07502)

    # 4. Verify SILCS parameters (from silcs.str)
    # LP parameters
    lp_params = ff.get_lj_params("LP")
    assert lp_params.epsilon == pytest.approx(0.0)
    assert lp_params.rmin == pytest.approx(0.0)

    # LQ parameters
    lq_params = ff.get_lj_params("LQ")
    assert lq_params.epsilon == pytest.approx(0.0)
    assert lq_params.rmin == pytest.approx(0.0)

    # SILCS NBFIX parameters
    epsilon, found = ff.get_nbfix("LP", "LP")
    assert found == True
    assert epsilon == pytest.approx(-0.01)

    epsilon, found = ff.get_nbfix("LQ", "LQ")
    assert found == True
    assert epsilon == pytest.approx(-0.01)

    # 5. Verify non-existent combinations
    # Between SILCS and ions
    epsilon, found = ff.get_nbfix("LP", "SOD")
    assert found == False

    epsilon, found = ff.get_nbfix("LQ", "CLA")
    assert found == False

    # Between ions
    epsilon, found = ff.get_nbfix("SOD", "POT")
    assert found == False

    # 6. Test parameter overriding and coexistence
    # Parse water_ions file again to ensure parameters are not duplicated or corrupted
    pygcmc.PrmParser.parse_file(water_ions_file, ff)
    
    # Water parameters should remain unchanged
    ht_params = ff.get_lj_params("HT")
    assert ht_params.epsilon == pytest.approx(-0.046)
    assert ht_params.rmin == pytest.approx(0.2245)
    
    # SILCS parameters should still be present
    epsilon, found = ff.get_nbfix("LP", "LP")
    assert found == True
    assert epsilon == pytest.approx(-0.01)

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
