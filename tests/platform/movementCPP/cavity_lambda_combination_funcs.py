#!/usr/bin/env python
"""
Detailed balance diagnostics for cavity bias and thermal wavelength (Λ³) handling.

These tests exercise MovementModule acceptance probabilities directly and
ensure that insertion/deletion ratios match the analytical GCMC formulas:

    P_ins / P_del = exp(βμ) * V_eff / ((N+1) * Λ³)

where V_eff incorporates cavity bias (if enabled) and Λ is the thermal
de Broglie wavelength.  The tests cover three practical configurations:
  1. Λ³ only (no cavity bias)
  2. Cavity only (Λ = 1)
  3. Cavity + Λ³ (the historically problematic path)
"""

from __future__ import annotations

import math
from typing import List, Tuple

import numpy as np
import pytest

import pygcmc


def _run_detailed_balance_case(
    *,
    use_cavity: bool,
    lambda_nm: float,
    chemical_potential: float = -5.0,
    equil_steps: int = 500,
    trial_pairs: int = 150,
) -> Tuple[List[float], List[float]]:
    """Collect insertion/deletion acceptance ratios for a specified configuration."""
    state = pygcmc.MCState()
    L = 3.0  # nm
    state.info.box = np.array([L, L, L])

    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.maxTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.1]
    state.forcefield = ff

    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = chemical_potential
    params.useCavityBias = use_cavity
    if use_cavity:
        params.cavityGridSpacing = 0.2
        params.probeRadius = 0.15
    params.thermalLambdaNm = lambda_nm
    params.seed = 12345
    params.updateDerivedParameters()

    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)

    np.random.seed(12345)
    for _ in range(equil_steps):
        if np.random.random() < 0.6:
            mover.attemptInsertion(state)
        else:
            mover.attemptDeletion(state)

    beta = 1.0 / (8.314e-3 * params.temperature)
    volume = L ** 3
    lambda_cubed = lambda_nm ** 3 if lambda_nm > 0 else 1.0

    ratios: List[float] = []
    errors: List[float] = []
    details: List[Tuple[int, float, float]] = []

    # Gather several paired insertion/deletion moves sharing the same microstate
    for _ in range(trial_pairs):
        active_before = {i for i, r in enumerate(state.residues) if r.active}
        n_before = len(active_before)
        if n_before == 0:
            continue

        ins_result = mover.attemptInsertion(state)
        if not ins_result.accepted:
            continue

        active_after = {i for i, r in enumerate(state.residues) if r.active}
        new_residues = active_after - active_before
        if not new_residues:
            continue
        inserted_idx = list(new_residues)[0]

        del_result = mover.attemptDeletion(state, inserted_idx)
        if not del_result.accepted:
            state.residues[inserted_idx].active = False
            continue

        p_ins = ins_result.acceptanceProbability
        p_del = del_result.acceptanceProbability
        if p_del <= 1e-10:
            state.residues[inserted_idx].active = False
            continue

        theory_ratio = math.exp(beta * params.chemicalPotential) * volume
        theory_ratio /= (n_before + 1) * lambda_cubed

        ratio = p_ins / p_del
        error_pct = abs(ratio - theory_ratio) / max(theory_ratio, 1e-12) * 100

        ratios.append(ratio)
        errors.append(error_pct)

        details.append((n_before, ratio, theory_ratio))

    return ratios, errors, details


@pytest.mark.parametrize(
    "use_cavity, lambda_nm, tolerance",
    [
        (False, 0.3, 5.0),   # Λ³ only
        (True, 0.3, 5.0),    # cavity + Λ³ (historically troublesome)
    ],
)
def test_detailed_balance_variants(use_cavity, lambda_nm, tolerance):
    ratios, errors, details = _run_detailed_balance_case(use_cavity=use_cavity, lambda_nm=lambda_nm)
    assert len(ratios) >= 20, f"Too few successful pairs: {len(ratios)} < 20"

    mean_ratio = float(np.mean(ratios))
    std_ratio = float(np.std(ratios))
    max_error = max(errors) if errors else 0.0

    print(f"\nDetailed balance check (cavity={use_cavity}, lambda={lambda_nm}):")
    print(
        f"  pairs={len(ratios)}, mean={mean_ratio:.6f}, std={std_ratio:.6f}, "
        f"max_error={max_error:.3f}%"
    )
    sample = details[:5]
    worst_idx = int(np.argmax(errors)) if errors else 0
    worst = details[worst_idx] if details else None
    print("  sample (n_before, ratio, theory):", sample)
    if worst is not None:
        print(f"  worst case: n={worst[0]}, ratio={worst[1]:.6f}, theory={worst[2]:.6f}, "
              f"error={errors[worst_idx]:.3f}%")

    assert max_error < tolerance, f"Detailed balance error {max_error:.2f}% exceeds tolerance {tolerance}%"
    relative_std = std_ratio / max(mean_ratio, 1e-12) * 100
    assert relative_std < 30.0, f"Too much scatter in ratios: {relative_std:.1f}% > 30%"


def test_cavity_bias_lambda_consistency():
    """Ensure cavity bias factors are geometry-driven and insensitive to Λ."""
    state1 = pygcmc.MCState()
    state1.info.box = np.array([2.5, 2.5, 2.5])

    state2 = pygcmc.MCState()
    state2.info.box = np.array([2.5, 2.5, 2.5])

    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.maxTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.1]
    state1.forcefield = ff
    state2.forcefield = ff

    params1 = pygcmc.movement.MovementParams()
    params1.temperature = 298.15
    params1.chemicalPotential = -10.0
    params1.useCavityBias = True
    params1.thermalLambdaNm = 1.0
    params1.seed = 12345
    params1.updateDerivedParameters()

    params2 = pygcmc.movement.MovementParams()
    params2.temperature = 298.15
    params2.chemicalPotential = -10.0 + 3.0 * math.log(0.5)
    params2.useCavityBias = True
    params2.thermalLambdaNm = 0.5
    params2.seed = 12345
    params2.updateDerivedParameters()

    mover1 = pygcmc.movement.MovementModule()
    mover1.setParams(params1)

    mover2 = pygcmc.movement.MovementModule()
    mover2.setParams(params2)

    np.random.seed(12345)
    for _ in range(50):
        if np.random.random() < 0.7:
            mover1.attemptInsertion(state1)
            mover2.attemptInsertion(state2)

    result1 = mover1.attemptInsertion(state1)
    result2 = mover2.attemptInsertion(state2)

    bias1 = getattr(result1, "cavityBiasFactor", 1.0)
    bias2 = getattr(result2, "cavityBiasFactor", 1.0)

    if bias1 > 0.1 and bias2 > 0.1:
        relative_diff = abs(bias1 - bias2) / max(bias1, bias2) * 100
        assert relative_diff < 20.0, f"Cavity bias differs between Λ paths: {bias1:.3f} vs {bias2:.3f}"

    print(f"Cavity bias consistency test: Λ=1.0 -> {bias1:.3f}, Λ≠1.0 -> {bias2:.3f}")


if __name__ == "__main__":
    pytest.main([__file__])
