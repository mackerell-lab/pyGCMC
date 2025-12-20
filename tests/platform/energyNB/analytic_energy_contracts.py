# tests/simulation/energyNB/analytic_energy_contracts.py
"""Analytic LJ/Coulomb energy contracts using energy/ API."""

import math
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


def _lj_energy(sigma, eps, r):
    sigma_r = sigma / r
    term6 = sigma_r ** 6
    term12 = term6 * term6
    return 4.0 * eps * (term12 - term6)


def test_system_energy_lj_only_analytic():
    r = 1.4
    sigma = 1.0
    eps = 0.6
    state = _build_two_residue_state(r, 0.0, 0.0, sigma, eps)

    total = pygcmc.computeSystemEnergy(state)
    expected = _lj_energy(sigma, eps, r)
    assert total == pytest.approx(expected, rel=1e-6, abs=1e-8)


def test_system_energy_coulomb_only_analytic():
    r = 0.8
    q1 = 0.7
    q2 = -0.3
    state = _build_two_residue_state(r, q1, q2, sigma=1.0, eps=0.0)

    total = pygcmc.computeSystemEnergy(state)
    expected = COULOMB * q1 * q2 / r
    assert total == pytest.approx(expected, rel=1e-6, abs=1e-8)


def test_system_energy_components_mixed_analytic():
    r = 1.1
    sigma = 0.9
    eps = 0.8
    q1 = 0.4
    q2 = -0.2
    state = _build_two_residue_state(r, q1, q2, sigma, eps)

    total = pygcmc.computeSystemEnergy(state)
    elec, vdw = pygcmc.getTotalEnergyComponents(state)

    expected_vdw = _lj_energy(sigma, eps, r)
    expected_elec = COULOMB * q1 * q2 / r
    expected_total = expected_vdw + expected_elec

    assert vdw == pytest.approx(expected_vdw, rel=1e-6, abs=1e-8)
    assert elec == pytest.approx(expected_elec, rel=1e-6, abs=1e-8)
    assert total == pytest.approx(expected_total, rel=1e-6, abs=1e-8)
