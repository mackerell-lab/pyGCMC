"""
Regression: ΔU from GCMC insertion/deletion must match unique-pairs total energy differences.

This locks in the "unique total energy" convention introduced to avoid accidental ×2
contamination when using total-energy difference paths.
"""

from __future__ import annotations

import math

import pytest
import pygcmc


def _make_minimal_state(*, box_nm: float, cutoff_nm: float) -> pygcmc.MCState:
    state = pygcmc.MCState()

    info = pygcmc.MCInfo()
    info.box = [float(box_nm), float(box_nm), float(box_nm)]
    info.cutoff = float(cutoff_nm)
    info.setTemperature(300.0)
    info.max_residues = 16
    info.max_atoms = 16
    info.volume = float(box_nm**3)
    state.info = info

    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.ljSigma = [0.30]  # nm
    ff.ljEps = [0.20]  # kJ/mol
    state.forcefield = ff

    state.residues = []
    state.atoms = []
    state.activeResidueCount = 0
    state.activeAtomCount = 0
    return state


def _make_single_atom_template(*, name: str, type_id: int) -> pygcmc.movement.FragmentTemplate:
    tmpl = pygcmc.movement.FragmentTemplate()
    tmpl.name = name
    tmpl.typeId = type_id

    atom = pygcmc.MCAtom()
    atom.type = type_id
    atom.charge = 0.0
    atom.x = 0.0
    atom.y = 0.0
    atom.z = 0.0
    tmpl.atoms = [atom]

    tmpl.isRigid = True
    tmpl.allowRotation = False
    tmpl.allowTranslation = True
    return tmpl


def _total_energy_unique_pairs_pbc_cutoff(state: pygcmc.MCState) -> float:
    pygcmc.computeSystemEnergyPBCCutoff(state)
    return pygcmc.getTotalEnergyUniquePairs(state, pygcmc.EnergyMethod.DIRECT)


def test_insert_then_delete_matches_unique_total_energy_diff() -> None:
    type_id = 0
    box_nm = 1.0
    cutoff_nm = 1.2  # ensure two random points always interact under MIC

    state = _make_minimal_state(box_nm=box_nm, cutoff_nm=cutoff_nm)

    reservoir = pygcmc.movement.FragmentReservoir()
    reservoir.addTemplate(_make_single_atom_template(name="X", type_id=type_id))

    engine = pygcmc.GCMCEngine()
    engine.initialize(state, reservoir)
    engine.setSeed(1234)
    engine.setConfigValue("useCavityBias", 0.0)
    engine.setConfigValue("useConfBias", 0.0)

    acc = pygcmc.GCMCAcceptance()
    acc.setTemperature(300.0)
    acc.setVolume(box_nm**3)
    acc.setThermalLambda(type_id, 1.0)
    engine.setAcceptanceCalculator(acc)

    # Insert first molecule (no pair interactions yet).
    acc.setActivity(type_id, 1e300)
    ins1 = engine.attemptInsertion(type_id)
    assert ins1.accepted
    e1 = _total_energy_unique_pairs_pbc_cutoff(state)
    assert math.isfinite(e1)

    # Insert second molecule and verify ΔU matches total energy diff.
    acc.setActivity(type_id, 1e300)
    ins2 = engine.attemptInsertion(type_id)
    assert ins2.accepted
    e2 = _total_energy_unique_pairs_pbc_cutoff(state)
    assert math.isfinite(e2)
    assert ins2.deltaE == pytest.approx(e2 - e1, rel=0, abs=1e-6)

    # Force deletion acceptance and verify ΔU matches total energy diff + cycle closes.
    acc.setActivity(type_id, 1e-300)
    dele = engine.attemptDeletion(type_id)
    assert dele.accepted
    e3 = _total_energy_unique_pairs_pbc_cutoff(state)
    assert math.isfinite(e3)
    assert dele.deltaE == pytest.approx(e3 - e2, rel=0, abs=1e-6)
    assert e3 == pytest.approx(e1, rel=0, abs=1e-6)
