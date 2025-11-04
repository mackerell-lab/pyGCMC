"""
GCMCEngine detailed balance tests for cavity bias + Lambda factors.

Mirrors the movementCPP cavity tests but exercises the engine pathway
to ensure insertion and deletion acceptance probabilities obey
the expected ratio on the same microstate.
"""

import math
import numpy as np
import pytest

import pygcmc
from pygcmc import movement

KB = 8.314e-3  # kJ/(mol*K)


def _build_simple_state(box_size: float = 3.0):
    """Create a minimal MCState and force field with a single fragment type."""
    state = pygcmc.MCState()
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.setTemperature(298.15)

    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.1]
    state.forcefield = ff

    reservoir = movement.FragmentReservoir()
    template = movement.FragmentTemplate()
    template.name = "test"
    template.typeId = 0

    atom = pygcmc.MCAtom()
    atom.x = atom.y = atom.z = 0.0
    atom.type = 0
    atom.charge = 0.0
    if hasattr(atom, "mass"):
        atom.mass = 18.0
    template.atoms = [atom]

    reservoir.addTemplate(template)
    return state, reservoir


@pytest.mark.parametrize(
    "use_cavity, lambda_nm",
    [
        (False, 0.3),  # Lambda only
        (True, 1.0),   # Cavity only
        (True, 0.3),   # Cavity + Lambda
    ],
)
def test_engine_cavity_lambda_detailed_balance(use_cavity: bool, lambda_nm: float):
    """Verify insertion/deletion ratios from GCMCEngine match theoretical expectations."""
    state, reservoir = _build_simple_state(box_size=3.0)

    engine = pygcmc.GCMCEngine()
    engine.initialize(state, reservoir)
    engine.setSeed(12345)
    engine.setTemperature(298.15)

    acceptance = pygcmc.GCMCAcceptance()
    acceptance.setTemperature(298.15)
    box_volume = float(np.prod(state.info.box))
    acceptance.setVolume(box_volume)

    mu = -5.0
    acceptance.setChemicalPotential(0, mu)
    acceptance.setThermalLambda(0, lambda_nm)
    acceptance.setSeed(12346)
    engine.setAcceptanceCalculator(acceptance)

    if use_cavity:
        engine.setConfigValue("useCavityBias", 1.0)
        cavity_mgr = pygcmc.CavityManager(2.0, 1.4)
        engine.setCavityManager(cavity_mgr)

    # Enable probability capture to retrieve detailed info
    engine.setConfigValue("storeProbabilities", 1.0)

    # Light equilibration to populate system with a few molecules
    for _ in range(400):
        if np.random.random() < 0.6:
            engine.attemptInsertion(0)
        else:
            engine.attemptDeletion(0)

    beta = 1.0 / (8.314e-3 * 298.15)
    lambda_cubed = lambda_nm ** 3

    ratios = []
    errors = []

    # Sample matched insertion/deletion attempts
    samples = 0
    attempts = 0
    max_attempts = 300
    while samples < 30 and attempts < max_attempts:
        attempts += 1
        n_before = reservoir.getActiveCount(0)
        if n_before == 0:
            engine.attemptInsertion(0)
            continue

        ins_result = engine.attemptInsertion(0)
        if not ins_result.accepted or ins_result.acceptanceProbability < 0:
            continue

        del_result = engine.attemptDeletion(0)
        if del_result.acceptanceProbability < 0:
            continue

        p_ins = ins_result.acceptanceProbability
        p_del = del_result.acceptanceProbability
        if p_ins <= 0 or p_del <= 0:
            continue

        v_eff = ins_result.effectiveVolume if ins_result.effectiveVolume > 0 else box_volume
        theory_ratio = math.exp(beta * mu) * v_eff / ((n_before + 1) * lambda_cubed)
        ratio = p_ins / p_del
        error_pct = abs(ratio - theory_ratio) / max(theory_ratio, 1e-12) * 100

        ratios.append(ratio)
        errors.append(error_pct)
        samples += 1

    assert samples >= 15, f"Insufficient matched samples collected (got {samples})"
    max_error = max(errors)
    mean_error = sum(errors) / len(errors)

    print(
        f"\nGCMCEngine cavity={use_cavity}, lambda={lambda_nm}: "
        f"samples={samples}, max_error={max_error:.3f}%, mean={mean_error:.3f}%"
    )

    assert max_error < 5.0, (
        f"Detailed balance error {max_error:.2f}% exceeds tolerance for "
        f"cavity={use_cavity}, lambda={lambda_nm}"
    )


def test_engine_acceptance_matches_formula():
    """Compare engine acceptance probabilities against closed-form expression."""
    state, reservoir = _build_simple_state(box_size=3.0)

    engine = pygcmc.GCMCEngine()
    engine.initialize(state, reservoir)
    engine.setSeed(77777)
    engine.setTemperature(298.15)

    acceptance = pygcmc.GCMCAcceptance()
    acceptance.setTemperature(298.15)
    volume = float(np.prod(state.info.box))
    acceptance.setVolume(volume)

    mu = -5.0
    lambda_nm = 0.6
    acceptance.setChemicalPotential(0, mu)
    acceptance.setThermalLambda(0, lambda_nm)
    acceptance.setSeed(77778)
    engine.setAcceptanceCalculator(acceptance)

    engine.setConfigValue("useCavityBias", 1.0)
    cavity_mgr = pygcmc.CavityManager(2.0, 1.4)
    engine.setCavityManager(cavity_mgr)
    engine.setConfigValue("storeProbabilities", 1.0)

    for _ in range(300):
        if np.random.random() < 0.7:
            engine.attemptInsertion(0)
        else:
            engine.attemptDeletion(0)

    beta = 1.0 / (KB * 298.15)
    activity = math.exp(beta * mu)
    lambda_cubed = lambda_nm ** 3
    proposal_bias = engine.getConfigValue("proposalBias")
    proposal_factor = 1.0
    if proposal_bias and proposal_bias > 0.0:
        proposal_factor = 1.0 / proposal_bias

    errors = []
    samples = 0
    max_attempts = 300

    for _ in range(max_attempts):
        n_before = reservoir.getActiveCount(0)
        if n_before == 0:
            engine.attemptInsertion(0)
            continue

        result = engine.attemptInsertion(0)
        p_actual = result.acceptanceProbability
        if p_actual < 0.0:
            continue

        v_eff = result.effectiveVolume if result.effectiveVolume > 0.0 else volume
        rosen = max(result.rosenbluthWeight, 1e-30)
        ratio = (
            proposal_factor
            * activity
            * v_eff
            / ((n_before + 1) * lambda_cubed)
        )
        ratio *= math.exp(-beta * result.deltaE) * rosen
        p_theory = 1.0 if ratio >= 1.0 else ratio

        errors.append(abs(p_actual - p_theory))
        samples += 1

        if samples >= 40:
            break

        # Keep population near steady-state to avoid runaway growth
        engine.attemptDeletion(0)

    assert samples >= 15, f"Collected only {samples} insertion samples"
    max_error = max(errors)
    mean_error = sum(errors) / len(errors)

    print(
        f"\nFormula verification: samples={samples}, "
        f"max_error={max_error:.2e}, mean_error={mean_error:.2e}"
    )

    assert max_error < 1e-4, f"Max error {max_error:.2e} exceeds tolerance 1e-4"
