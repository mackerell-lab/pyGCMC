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
