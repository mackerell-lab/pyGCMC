"""
Regression tests for Thole screening application.

Goal: lock in CHARMM/OpenMM-style "point charge" Thole screening semantics:
- StandardS1: apply S1(u) to all 4 charge-charge interactions between two Drude dipoles
             (P1-P2, P1-D2, D1-P2, D1-D2) using per-interaction distance.
- OpenMMCompat: apply S1(u) only to interactions involving at least one Drude particle
                (P1-D2, D1-P2, D1-D2), leaving P1-P2 unscreened.

These tests are intentionally "non-cheating":
- Use an independent Python implementation of S1(u).
- Disable SCF updates (maxIterations=0) so the system geometry is fixed and
  the expected energy can be computed analytically.
"""

from __future__ import annotations

import math

import pygcmc


def _s1(u: float) -> float:
    if u <= 0.0:
        return 0.0
    if u > 50.0:
        return 1.0
    return 1.0 - (1.0 + 0.5 * u) * math.exp(-u)


def _s1_from_r(r: float, alpha_i: float, alpha_j: float, thole: float) -> float:
    if thole == 0.0:
        return 1.0
    alpha_eff = (alpha_i * alpha_j) ** (1.0 / 6.0)
    u = thole * r / alpha_eff
    return _s1(u)


def _dist(a: pygcmc.MCAtom, b: pygcmc.MCAtom) -> float:
    dx = a.x - b.x
    dy = a.y - b.y
    dz = a.z - b.z
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _spring_energy(p: pygcmc.DrudeParticle, parent: pygcmc.MCAtom, drude: pygcmc.MCAtom) -> float:
    dx = drude.x - parent.x
    dy = drude.y - parent.y
    dz = drude.z - parent.z
    return 0.5 * p.kSpring * (dx * dx + dy * dy + dz * dz)


def _expected_screened_coulomb_energy(
    p1: pygcmc.MCAtom,
    d1: pygcmc.MCAtom,
    p2: pygcmc.MCAtom,
    d2: pygcmc.MCAtom,
    alpha1: float,
    alpha2: float,
    thole: float,
    mode: pygcmc.TholeMode,
) -> float:
    k_e = pygcmc.DrudeConstants.ONE_4PI_EPS0

    def screened(qi: float, qj: float, r: float, involve_drude: bool) -> float:
        if r < 1e-12:
            return 0.0
        if mode == pygcmc.TholeMode.OpenMMCompat and not involve_drude:
            s = 1.0
        else:
            s = _s1_from_r(r, alpha1, alpha2, thole)
        return k_e * qi * qj * s / r

    e = 0.0
    e += screened(p1.charge, p2.charge, _dist(p1, p2), involve_drude=False)
    e += screened(p1.charge, d2.charge, _dist(p1, d2), involve_drude=True)
    e += screened(d1.charge, p2.charge, _dist(d1, p2), involve_drude=True)
    e += screened(d1.charge, d2.charge, _dist(d1, d2), involve_drude=True)
    return e


def _build_two_dipole_state() -> tuple[pygcmc.MCState, pygcmc.DrudeParticle, pygcmc.DrudeParticle]:
    """
    Build a fixed-geometry 2-dipole system in a large box (no PBC effects).

    Geometry is chosen so the parent-parent and cross terms are significant,
    making it impossible for a "Drude-Drude only" screening implementation to pass.
    """
    state = pygcmc.MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 4.0

    p1 = pygcmc.MCAtom()
    p1.x, p1.y, p1.z = 0.0, 0.0, 0.0
    p1.charge = 1.0
    p1.type = 0

    d1 = pygcmc.MCAtom()
    d1.x, d1.y, d1.z = 0.05, 0.0, 0.0
    d1.charge = -1.0
    d1.type = 1

    p2 = pygcmc.MCAtom()
    p2.x, p2.y, p2.z = 0.20, 0.00, 0.0
    p2.charge = 1.0
    p2.type = 0

    d2 = pygcmc.MCAtom()
    d2.x, d2.y, d2.z = 0.15, 0.02, 0.0
    d2.charge = -1.0
    d2.type = 1

    state.atoms = [p1, d1, p2, d2]
    state.activeAtomCount = 4

    # Define two "molecules" so intramolecular parent-drude interactions are excluded.
    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 2
    res1.active = True
    res1.type = 0

    res2 = pygcmc.MCResidue()
    res2.atomStart = 2
    res2.atomCount = 2
    res2.active = True
    res2.type = 1

    state.residues = [res1, res2]
    state.activeResidueCount = 2

    # Add two Drude particles matching the atom indices and charges.
    alpha = 0.001  # nm^3

    dp1 = pygcmc.DrudeParticle()
    dp1.drudeIndex = 1
    dp1.parentIndex = 0
    dp1.charge = d1.charge
    dp1.polarizability = alpha
    dp1.computeSpringConstants()

    dp2 = pygcmc.DrudeParticle()
    dp2.drudeIndex = 3
    dp2.parentIndex = 2
    dp2.charge = d2.charge
    dp2.polarizability = alpha
    dp2.computeSpringConstants()

    return state, dp1, dp2


def test_thole_screening_energy_standard_s1_uses_four_charge_pairs():
    state, dp1, dp2 = _build_two_dipole_state()

    pygcmc.DrudeComplete.clear()
    pygcmc.DrudeComplete.addParticle(dp1)
    pygcmc.DrudeComplete.addParticle(dp2)

    pair = pygcmc.ScreenedPair()
    pair.dipole1 = 0
    pair.dipole2 = 1
    pair.thole = 2.0
    pygcmc.DrudeComplete.addScreenedPair(pair)

    params = pygcmc.DrudeSCFParams()
    params.includeCoulombEnergy = True
    params.maxIterations = 0  # freeze geometry: test pure energy formula
    params.tholeMode = pygcmc.TholeMode.StandardS1
    pygcmc.DrudeComplete.setParameters(params)

    energy = pygcmc.DrudeComplete.calculateEnergy(state)

    p1, d1, p2, d2 = state.atoms
    spring = _spring_energy(dp1, p1, d1) + _spring_energy(dp2, p2, d2)
    expected_coulomb = _expected_screened_coulomb_energy(
        p1, d1, p2, d2, dp1.polarizability, dp2.polarizability, pair.thole, pygcmc.TholeMode.StandardS1
    )
    expected_total = spring + expected_coulomb

    assert abs(energy - expected_total) < 1e-8 * max(1.0, abs(expected_total))

    pygcmc.DrudeComplete.clear()


def test_thole_screening_energy_openmmcompat_screens_only_drude_involving_pairs():
    state, dp1, dp2 = _build_two_dipole_state()

    pygcmc.DrudeComplete.clear()
    pygcmc.DrudeComplete.addParticle(dp1)
    pygcmc.DrudeComplete.addParticle(dp2)

    pair = pygcmc.ScreenedPair()
    pair.dipole1 = 0
    pair.dipole2 = 1
    pair.thole = 2.0
    pygcmc.DrudeComplete.addScreenedPair(pair)

    params = pygcmc.DrudeSCFParams()
    params.includeCoulombEnergy = True
    params.maxIterations = 0  # freeze geometry: test pure energy formula
    params.tholeMode = pygcmc.TholeMode.OpenMMCompat
    pygcmc.DrudeComplete.setParameters(params)

    energy = pygcmc.DrudeComplete.calculateEnergy(state)

    p1, d1, p2, d2 = state.atoms
    spring = _spring_energy(dp1, p1, d1) + _spring_energy(dp2, p2, d2)
    expected_coulomb = _expected_screened_coulomb_energy(
        p1, d1, p2, d2, dp1.polarizability, dp2.polarizability, pair.thole, pygcmc.TholeMode.OpenMMCompat
    )
    expected_total = spring + expected_coulomb

    assert abs(energy - expected_total) < 1e-8 * max(1.0, abs(expected_total))

    pygcmc.DrudeComplete.clear()

