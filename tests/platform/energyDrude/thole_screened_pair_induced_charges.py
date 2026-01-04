"""
Regression: Thole screened-pair interactions must operate on induced Drude charges (±q_drude),
not on the full core atomic charges (which include permanent charge), matching OpenMM/CHARMM
DrudeForce semantics.

Additionally, because gcmc_cpu's DIRECT nonbonded backend excludes intra-residue Coulomb
interactions entirely, screened pairs within a residue must contribute their *full screened*
interaction energy (not just a correction term).

These tests are intentionally "non-cheating":
- Use an independent Python implementation of S1(u).
- Freeze SCF (maxIterations=0) so the geometry is fixed and the expected energy can be computed
  analytically.
"""

from __future__ import annotations

import math

import pygcmc
import pytest


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


def _expected_induced_four_charge_energy(
    p1: pygcmc.MCAtom,
    d1: pygcmc.MCAtom,
    p2: pygcmc.MCAtom,
    d2: pygcmc.MCAtom,
    *,
    qd1: float,
    qd2: float,
    alpha1: float,
    alpha2: float,
    thole: float,
    factor_mode: str,
) -> float:
    """
    Compute the screened-pair energy using induced charges only:
      parent carries -q_drude, drude carries +q_drude.

    factor_mode:
      - "screened": add full screened Coulomb (baseline absent, intra-residue in gcmc_cpu DIRECT).
      - "delta": add (screened - unscreened) correction (baseline present, inter-residue).
    """
    assert factor_mode in ("screened", "delta")
    k_e = pygcmc.DrudeConstants.ONE_4PI_EPS0

    # Induced charges for the dipole pair.
    q_p1 = -qd1
    q_d1 = qd1
    q_p2 = -qd2
    q_d2 = qd2

    def add(qi: float, ai: pygcmc.MCAtom, qj: float, aj: pygcmc.MCAtom) -> float:
        r = _dist(ai, aj)
        if r < 1e-12:
            return 0.0
        s = _s1_from_r(r, alpha1, alpha2, thole)
        base = k_e * qi * qj / r
        if factor_mode == "screened":
            return base * s
        return base * (s - 1.0)

    e = 0.0
    # 4 induced charge-charge interactions between dipoles.
    e += add(q_p1, p1, q_p2, p2)
    e += add(q_p1, p1, q_d2, d2)
    e += add(q_d1, d1, q_p2, p2)
    e += add(q_d1, d1, q_d2, d2)
    return e


def _build_two_dipoles(
    *, qd1: float, qperm1: float, qd2: float, qperm2: float, same_residue: bool
) -> tuple[pygcmc.MCState, pygcmc.DrudeParticle, pygcmc.DrudeParticle]:
    """
    Build a fixed-geometry 2-dipole system with nonzero permanent charge embedded in the core charge:
      q_core = q_perm - q_drude
      q_drude = qd

    The Drude screened-pair interaction must depend only on qd1/qd2, not on qperm1/qperm2.
    """
    state = pygcmc.MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 4.0

    # Geometry in nm (no PBC effects within this box).
    p1 = pygcmc.MCAtom()
    p1.x, p1.y, p1.z = 0.00, 0.00, 0.00
    p1.charge = qperm1 - qd1
    p1.type = 0

    d1 = pygcmc.MCAtom()
    d1.x, d1.y, d1.z = 0.05, 0.00, 0.00
    d1.charge = qd1
    d1.type = 1

    p2 = pygcmc.MCAtom()
    p2.x, p2.y, p2.z = 0.30, 0.00, 0.00
    p2.charge = qperm2 - qd2
    p2.type = 0

    d2 = pygcmc.MCAtom()
    d2.x, d2.y, d2.z = 0.25, 0.02, 0.00
    d2.charge = qd2
    d2.type = 1

    state.atoms = [p1, d1, p2, d2]
    state.activeAtomCount = 4

    if same_residue:
        res = pygcmc.MCResidue()
        res.atomStart = 0
        res.atomCount = 4
        res.active = True
        res.type = 0
        state.residues = [res]
        state.activeResidueCount = 1
    else:
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

    alpha = 0.001  # nm^3

    dp1 = pygcmc.DrudeParticle()
    dp1.parentIndex = 0
    dp1.drudeIndex = 1
    dp1.charge = qd1
    dp1.polarizability = alpha
    dp1.computeSpringConstants()

    dp2 = pygcmc.DrudeParticle()
    dp2.parentIndex = 2
    dp2.drudeIndex = 3
    dp2.charge = qd2
    dp2.polarizability = alpha
    dp2.computeSpringConstants()

    return state, dp1, dp2


def test_thole_induced_charge_semantics_inter_residue_is_core_charge_independent():
    """
    Inter-residue screened pair: DrudeComplete (includeCoulombEnergy=False) must return only the
    correction term (screened - unscreened) for the induced charge pair, independent of permanent
    charge embedded in the core charge.
    """
    qd1, qd2 = -1.2, -0.8
    alpha = 0.001
    thole = 2.6

    state_a, dp1_a, dp2_a = _build_two_dipoles(qd1=qd1, qperm1=-0.3, qd2=qd2, qperm2=0.1, same_residue=False)
    state_b, dp1_b, dp2_b = _build_two_dipoles(qd1=qd1, qperm1=0.9, qd2=qd2, qperm2=-0.7, same_residue=False)

    def run(state: pygcmc.MCState, dp1: pygcmc.DrudeParticle, dp2: pygcmc.DrudeParticle) -> float:
        pygcmc.DrudeComplete.clear()
        pygcmc.DrudeComplete.addParticle(dp1)
        pygcmc.DrudeComplete.addParticle(dp2)
        pair = pygcmc.ScreenedPair()
        pair.dipole1 = 0
        pair.dipole2 = 1
        pair.thole = thole
        pygcmc.DrudeComplete.addScreenedPair(pair)

        params = pygcmc.DrudeSCFParams()
        params.includeCoulombEnergy = False
        params.maxIterations = 0
        params.tholeMode = pygcmc.TholeMode.StandardS1
        pygcmc.DrudeComplete.setParameters(params)
        return pygcmc.DrudeComplete.calculateEnergy(state)

    e_a = run(state_a, dp1_a, dp2_a)
    e_b = run(state_b, dp1_b, dp2_b)
    assert abs(e_a - e_b) < 1e-10 * max(1.0, abs(e_a), abs(e_b))

    p1, d1, p2, d2 = state_a.atoms
    spring = _spring_energy(dp1_a, p1, d1) + _spring_energy(dp2_a, p2, d2)
    expected_thole_delta = _expected_induced_four_charge_energy(
        p1,
        d1,
        p2,
        d2,
        qd1=qd1,
        qd2=qd2,
        alpha1=alpha,
        alpha2=alpha,
        thole=thole,
        factor_mode="delta",
    )
    expected_total = spring + expected_thole_delta
    assert e_a == pytest.approx(expected_total, abs=2e-4, rel=1e-10)


def test_thole_induced_charge_semantics_intra_residue_adds_full_screened_energy():
    """
    Intra-residue screened pair: gcmc_cpu DIRECT excludes intra-residue Coulomb, so DrudeComplete
    must contribute the full screened induced-charge interaction (not just a delta term).
    """
    qd1, qd2 = -1.2, -0.8
    alpha = 0.001
    thole = 2.6

    state_a, dp1_a, dp2_a = _build_two_dipoles(qd1=qd1, qperm1=-0.3, qd2=qd2, qperm2=0.1, same_residue=True)
    state_b, dp1_b, dp2_b = _build_two_dipoles(qd1=qd1, qperm1=0.9, qd2=qd2, qperm2=-0.7, same_residue=True)

    def run(state: pygcmc.MCState, dp1: pygcmc.DrudeParticle, dp2: pygcmc.DrudeParticle) -> float:
        pygcmc.DrudeComplete.clear()
        pygcmc.DrudeComplete.addParticle(dp1)
        pygcmc.DrudeComplete.addParticle(dp2)
        pair = pygcmc.ScreenedPair()
        pair.dipole1 = 0
        pair.dipole2 = 1
        pair.thole = thole
        pygcmc.DrudeComplete.addScreenedPair(pair)

        params = pygcmc.DrudeSCFParams()
        params.includeCoulombEnergy = False
        params.maxIterations = 0
        params.tholeMode = pygcmc.TholeMode.StandardS1
        pygcmc.DrudeComplete.setParameters(params)
        return pygcmc.DrudeComplete.calculateEnergy(state)

    e_a = run(state_a, dp1_a, dp2_a)
    e_b = run(state_b, dp1_b, dp2_b)
    assert abs(e_a - e_b) < 1e-10 * max(1.0, abs(e_a), abs(e_b))

    p1, d1, p2, d2 = state_a.atoms
    spring = _spring_energy(dp1_a, p1, d1) + _spring_energy(dp2_a, p2, d2)
    expected_thole = _expected_induced_four_charge_energy(
        p1,
        d1,
        p2,
        d2,
        qd1=qd1,
        qd2=qd2,
        alpha1=alpha,
        alpha2=alpha,
        thole=thole,
        factor_mode="screened",
    )
    expected_total = spring + expected_thole
    assert e_a == pytest.approx(expected_total, abs=2e-4, rel=1e-10)
