# tests/io/prmParser/charmm_protein.py
"""PRM Parser CHARMM protein force field tests."""

import os
import pytest
import pygcmc

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_charmm_prm_files():
    """Test parsing of CHARMM force field files (par_all36m_prot.prm)."""
    prm_file = os.path.join(TEST_DATA_DIR, "par_all36m_prot.prm")

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
    prm_file = os.path.join(TEST_DATA_DIR, "par_all36_cgenff.prm")

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