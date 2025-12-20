# tests/simulation/energyNB/cutoff_pbc_contracts.py
"""Cutoff and PBC minimum-image contracts using energy/ API."""

import pytest
import pygcmc

COULOMB = 138.935458  # kJ·nm/mol/e^2


def _build_two_residue_state(r, q1, q2, sigma, eps):
    state = pygcmc.MCState()
    state.forcefield.numTotalTypes = 1
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [eps]
    state.forcefield.ljSigma = [sigma]
    state.movementAtomTypes = [0]
    state.numMovementAtomTypes = 1

    atom1 = pygcmc.MCAtom()
    atom1.x = 0.0
    atom1.y = 0.0
    atom1.z = 0.0
    atom1.charge = q1
    atom1.type = 0

    atom2 = pygcmc.MCAtom()
    atom2.x = r
    atom2.y = 0.0
    atom2.z = 0.0
    atom2.charge = q2
    atom2.type = 0

    state.atoms = [atom1, atom2]
    state.activeAtomCount = 2

    res1 = pygcmc.MCResidue()
    res1.active = True
    res1.type = 0
    res1.atomStart = 0
    res1.atomCount = 1

    res2 = pygcmc.MCResidue()
    res2.active = True
    res2.type = 0
    res2.atomStart = 1
    res2.atomCount = 1

    state.residues = [res1, res2]
    state.activeResidueCount = 2
    return state


def test_cutoff_excludes_far_pair():
    cutoff = 1.0
    q1 = 0.5
    q2 = -0.5

    near_state = _build_two_residue_state(0.8, q1, q2, sigma=1.0, eps=0.0)
    near_state.info.cutoff = cutoff
    near_energy = pygcmc.computeSystemEnergyCutoff(near_state)
    expected_near = COULOMB * q1 * q2 / 0.8
    assert near_energy == pytest.approx(expected_near, rel=1e-6, abs=1e-8)

    far_state = _build_two_residue_state(1.2, q1, q2, sigma=1.0, eps=0.0)
    far_state.info.cutoff = cutoff
    far_energy = pygcmc.computeSystemEnergyCutoff(far_state)
    assert far_energy == pytest.approx(0.0, abs=1e-12)


def test_pbc_minimum_image_distance():
    box = 10.0
    q1 = 0.1
    q2 = -0.1
    state = _build_two_residue_state(0.0, q1, q2, sigma=1.0, eps=0.0)
    state.atoms[0].x = 0.1
    state.atoms[1].x = 9.9
    state.info.box = [box, box, box]

    pygcmc.computeSystemEnergyPBC(state)
    elec, vdw = pygcmc.getTotalEnergyComponents(state)
    total = elec + vdw

    expected_r = 0.2
    expected = COULOMB * q1 * q2 / expected_r
    assert total == pytest.approx(expected, rel=1e-5, abs=1e-7)
