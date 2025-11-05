"""
Verify that CBMC Rosenbluth weights (W/K) remain invariant w.r.t. the number of
trials K when all trial energies are identical.
"""

import math
import numpy as np

import pytest

import pygcmc
from pygcmc import movement

KB = 8.314e-3  # kJ/(mol*K)


def _build_zero_energy_state(box_size: float = 3.0):
    state = pygcmc.MCState()
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.setTemperature(298.15)

    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]   # No interactions => all trial energies identical (0)
    ff.ljSigma = [0.1]
    state.forcefield = ff

    reservoir = movement.FragmentReservoir()
    template = movement.FragmentTemplate()
    template.name = "neutral"
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


def _setup_engine(cbmc_trials: int):
    state, reservoir = _build_zero_energy_state()

    engine = pygcmc.GCMCEngine()
    engine.initialize(state, reservoir)
    engine.setSeed(424242 + cbmc_trials)
    engine.setTemperature(298.15)

    acceptance = pygcmc.GCMCAcceptance()
    acceptance.setTemperature(298.15)
    volume = float(np.prod(state.info.box))
    acceptance.setVolume(volume)
    acceptance.setChemicalPotential(0, -6.0)
    acceptance.setThermalLambda(0, 1.0)
    acceptance.setSeed(424200 + cbmc_trials)
    engine.setAcceptanceCalculator(acceptance)

    engine.setConfigValue("storeProbabilities", 0.0)

    if cbmc_trials > 1:
        engine.setConfigValue("useConfBias", 1.0)
    engine.setCBMCTrialsPerType([cbmc_trials])
    else:
        engine.setConfigValue("useConfBias", 0.0)

    return engine, reservoir


@pytest.mark.parametrize("cbmc_trials", [2, 8])
def test_cbmc_rosenbluth_weight_invariant(cbmc_trials: int):
    """When trial energies are identical, W/K should be exactly 1 for any K."""
    engine, reservoir = _setup_engine(cbmc_trials)

    insertion_weight = None
    for _ in range(50):
        result = engine.attemptInsertion(0)
        if result.rosenbluthWeight > 0.0:
            insertion_weight = result.rosenbluthWeight
        if result.accepted:
            break

    assert insertion_weight is not None, f"No insertion move recorded for K={cbmc_trials}"
    assert abs(insertion_weight - 1.0) < 1e-12, (
        f"Rosenbluth weight deviated from unity (insertion, K={cbmc_trials}): "
        f"{insertion_weight}"
    )

    assert (
        reservoir.getActiveCount(0) > 0
    ), "Insertion never succeeded; cannot check deletion invariance"

    deletion_weight = None
    for _ in range(50):
        result = engine.attemptDeletion(0)
        if result.rosenbluthWeight > 0.0:
            deletion_weight = result.rosenbluthWeight
            break

    assert deletion_weight is not None, f"No deletion move recorded for K={cbmc_trials}"
    assert abs(deletion_weight - 1.0) < 1e-12, (
        f"Rosenbluth weight deviated from unity (deletion, K={cbmc_trials}): "
        f"{deletion_weight}"
    )
