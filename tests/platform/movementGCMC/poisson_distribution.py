"""
Poisson distribution check for ideal-gas GCMC (engine-level).

We validate <N> = z * V and Var(N) = <N> for a single-atom fragment with
zero interactions (DIRECT + cutoff path), using symmetric insertion/deletion
proposals.
"""

from __future__ import annotations

import math
import random
import statistics
import sys
from pathlib import Path

def _ensure_pygcmc_on_path() -> None:
    repo_root = Path(__file__).resolve().parents[3]
    bindings = repo_root / "build" / "modules" / "bindings"
    if bindings.exists():
        bindings_str = str(bindings)
        if bindings_str not in sys.path:
            sys.path.insert(0, bindings_str)


_ensure_pygcmc_on_path()
try:
    import pygcmc
except ModuleNotFoundError as exc:
    raise AssertionError(
        "pygcmc bindings not found; build with `cmake --build build --target pygcmc`"
    ) from exc


def test_ideal_gas_poisson_distribution():
    temperature = 300.0
    box_size_nm = 3.0
    volume = box_size_nm ** 3
    activity = 0.2  # nm^-3; expected mean = 0.2 * 27 = 5.4

    state = pygcmc.MCState()
    state.info.box = (box_size_nm, box_size_nm, box_size_nm)
    state.info.setTemperature(temperature)

    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljSigma = [0.3]
    ff.ljEps = [0.0]
    state.forcefield = ff

    tmpl = pygcmc.movement.FragmentTemplate()
    tmpl.typeId = 0
    atom = pygcmc.MCAtom()
    atom.type = 0
    atom.charge = 0.0
    tmpl.atoms = [atom]

    reservoir = pygcmc.movement.FragmentReservoir()
    reservoir.addTemplate(tmpl)

    engine = pygcmc.GCMCEngine()
    engine.initialize(state, reservoir)
    engine.setTemperature(temperature)
    engine.setSeed(2025)

    acceptance = pygcmc.GCMCAcceptance()
    acceptance.setTemperature(temperature)
    acceptance.setVolume(volume)
    acceptance.setActivity(0, activity)
    acceptance.setSeed(2026)
    engine.setAcceptanceCalculator(acceptance)

    n_steps = 10000
    burn_in = 1000
    stride = 5
    samples: list[int] = []

    rng = random.Random(2027)
    for step in range(n_steps):
        if rng.random() < 0.5:
            engine.attemptInsertion(0)
        else:
            engine.attemptDeletion(0)
        if step >= burn_in and step % stride == 0:
            samples.append(reservoir.getActiveCount(0))

    assert len(samples) >= 1000

    expected_mean = activity * volume
    sample_mean = statistics.mean(samples)
    sample_var = statistics.pvariance(samples)

    std_err = math.sqrt(expected_mean) / math.sqrt(len(samples))
    assert abs(sample_mean - expected_mean) < 4.0 * std_err + 0.1
    assert abs(sample_var - expected_mean) / expected_mean < 0.3
