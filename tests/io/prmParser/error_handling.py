# tests/io/prmParser/error_handling.py
"""PRM Parser error handling tests."""

import pytest
import pygcmc


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