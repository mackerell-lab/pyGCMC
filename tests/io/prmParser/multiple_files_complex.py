# tests/io/prmParser/multiple_files_complex.py
"""PRM Parser complex multiple files handling tests."""

import os
import pytest
import pygcmc

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_multiple_parameter_files():
    """Test reading and combining multiple parameter files."""
    water_ions_file = os.path.join(TEST_DATA_DIR, "toppar_water_ions.str")
    silcs_file = os.path.join(TEST_DATA_DIR, "silcs.str")
    prot_file = os.path.join(TEST_DATA_DIR, "par_all36m_prot.prm")

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
    # Load .str files
    water_ions_file = os.path.join(TEST_DATA_DIR, "toppar_water_ions.str")
    silcs_file = os.path.join(TEST_DATA_DIR, "silcs.str")

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