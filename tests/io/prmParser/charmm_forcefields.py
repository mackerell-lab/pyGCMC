# tests/io/prmParser/charmm_forcefields.py
"""PRM Parser CHARMM force field interaction tests."""

import os
import pytest
import pygcmc

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_ion_ligand_nbfix():
    """Test ion-ligand NBFIX parameters from toppar_water_ions.str."""
    param_file = os.path.join(TEST_DATA_DIR, "toppar_water_ions.str")

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
    prm_file = os.path.join(TEST_DATA_DIR, "par_all36_cgenff.prm")

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
    prm_file = os.path.join(TEST_DATA_DIR, "par_all36_cgenff.prm")

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
    prot_file = os.path.join(TEST_DATA_DIR, "par_all36m_prot.prm")
    cgenff_file = os.path.join(TEST_DATA_DIR, "par_all36_cgenff.prm")

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