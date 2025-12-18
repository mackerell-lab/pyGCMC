"""
P1 regression: strict_region_balance should prevent deletion fallback outside a region.

When a region constraint is active, deletion selection filters instances by whether their
reference position lies inside the region. The legacy fallback-to-all behavior can violate
detailed balance if the target distribution is defined only over the constrained region.
"""

from __future__ import annotations

import numpy as np

import pygcmc


def _build_engine_with_single_instance(*, strict: bool) -> tuple[pygcmc.GCMCEngine, pygcmc.movement.FragmentReservoir]:
    box_size = 3.0  # nm

    state = pygcmc.MCState()
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.setTemperature(298.15)

    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.1]
    state.forcefield = ff

    reservoir = pygcmc.movement.FragmentReservoir()
    tmpl = pygcmc.movement.FragmentTemplate()
    tmpl.name = "SOL"
    tmpl.typeId = 0

    atom = pygcmc.MCAtom()
    atom.x = atom.y = atom.z = 0.0
    atom.type = 0
    atom.charge = 0.0
    if hasattr(atom, "mass"):
        atom.mass = 18.0
    tmpl.atoms = [atom]
    reservoir.addTemplate(tmpl)

    engine = pygcmc.GCMCEngine()
    engine.initialize(state, reservoir)
    engine.setSeed(12345)
    engine.setTemperature(298.15)

    acceptance = pygcmc.GCMCAcceptance()
    acceptance.setTemperature(298.15)
    acceptance.setVolume(float(np.prod(state.info.box)))
    acceptance.setActivity(0, 1e-3)  # ensure deletion is in prob=1 regime for N=1
    engine.setAcceptanceCalculator(acceptance)

    engine.setConfigValue("strict_region_balance", 1.0 if strict else 0.0)

    # Create one active instance outside the constrained region and sync into MCState.
    outside = pygcmc.movement.Vector3(2.7, 2.7, 2.7)
    instance_id = reservoir.createInstance(0, outside)
    assert instance_id >= 0
    engine.synchronizeStateWithReservoir(instance_id, True)

    # Small sphere near origin (nm) which excludes the instance above.
    engine.setRegionConstraintFromSpec(
        "sphere 0.3 0.3 0.3 0.2",
        box_size,
        box_size,
        box_size,
    )

    assert reservoir.getActiveCount(0) == 1
    return engine, reservoir


def test_strict_region_balance_blocks_deletion_fallback_when_region_empty() -> None:
    engine_relaxed, reservoir_relaxed = _build_engine_with_single_instance(strict=False)
    relaxed = engine_relaxed.attemptDeletion(0)
    assert relaxed.accepted is True
    assert reservoir_relaxed.getActiveCount(0) == 0

    engine_strict, reservoir_strict = _build_engine_with_single_instance(strict=True)
    strict = engine_strict.attemptDeletion(0)
    assert strict.accepted is False
    assert reservoir_strict.getActiveCount(0) == 1
