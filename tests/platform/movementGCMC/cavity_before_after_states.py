"""
Ensure cavity bias is evaluated in the correct state for insertion/deletion moves.

Insertion should use the pre-insertion (before) cavity volume fraction,
while deletion should use the post-deletion (after) cavity fraction.
"""

import numpy as np
import pytest

import pygcmc
from pygcmc import movement


def _set_up_engine(box_size: float = 3.5):
    """Create a small engine with cavity bias enabled."""
    state = pygcmc.MCState()
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.setTemperature(300.0)

    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.2]
    ff.ljSigma = [0.3]
    state.forcefield = ff

    reservoir = movement.FragmentReservoir()
    template = movement.FragmentTemplate()
    template.name = "mono"
    template.typeId = 0

    atom = pygcmc.MCAtom()
    atom.type = 0
    atom.charge = 0.0
    atom.x = atom.y = atom.z = 0.0
    template.atoms = [atom]
    reservoir.addTemplate(template)

    engine = pygcmc.GCMCEngine()
    engine.initialize(state, reservoir)
    engine.setSeed(1234)
    engine.setTemperature(300.0)

    acceptance = pygcmc.GCMCAcceptance()
    acceptance.setTemperature(300.0)
    acceptance.setVolume(box_size ** 3)
    acceptance.setActivity(0, 500.0)  # encourage insertions
    engine.setAcceptanceCalculator(acceptance)

    cavity_mgr = pygcmc.CavityManager(1.0, 0.25)
    engine.setCavityManager(cavity_mgr)
    engine.setConfigValue("useCavityBias", 1.0)
    engine.setConfigValue("storeProbabilities", 1.0)

    return engine, state, reservoir, cavity_mgr


def test_insertion_uses_before_state_cavity():
    engine, state, reservoir, cavity_mgr = _set_up_engine()

    # Populate system a bit to make cavity meaningful
    for _ in range(20):
        engine.attemptInsertion(0)

    assert reservoir.getActiveCount(0) > 0, "Failed to populate system"

    before_fraction = cavity_mgr.getCavityVolumeFraction(state)
    result = engine.attemptInsertion(0)

    assert result.acceptanceProbability >= 0.0, "Probabilities are not being stored"
    assert before_fraction > 0.0
    assert before_fraction <= 1.0 + 1e-9

    assert pytest.approx(before_fraction, rel=1e-6) == result.cavityBiasComponent


def test_deletion_uses_after_state_cavity():
    engine, state, reservoir, cavity_mgr = _set_up_engine()

    # Ensure we have several molecules in the box first
    for _ in range(30):
        engine.attemptInsertion(0)

    assert reservoir.getActiveCount(0) > 5, "System too sparse for deletion test"

    # Attempt deletions until one is accepted (to observe after-state directly)
    for _ in range(40):
        result = engine.attemptDeletion(0)
        if not result.accepted:
            continue

        after_fraction = cavity_mgr.getCavityVolumeFraction(state)
        assert result.acceptanceProbability >= 0.0
        assert after_fraction > 0.0
        assert after_fraction <= 1.0 + 1e-9
        assert pytest.approx(after_fraction, rel=1e-6) == result.cavityBiasComponent
        break
    else:
        pytest.skip("No deletion was accepted; cannot verify after-state cavity.")
