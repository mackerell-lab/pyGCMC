"""
Total energy convention tests: unique-pairs vs double-counted residue sums.

Goals:
- Lock in the existing per-residue "one-to-all" storage convention for DIRECT energies:
  sum(residue energies) is double-counted for pair terms.
- Validate the dedicated unique-pairs total energy helper used by GCMC non-DIRECT ΔU paths.
"""

import pytest
import pygcmc


def _make_two_residue_system(distance=0.3):
    state = pygcmc.MCState()
    state.info.box = [4.0, 4.0, 4.0]
    state.info.cutoff = 1.2

    atom1 = pygcmc.MCAtom()
    atom1.x, atom1.y, atom1.z = 0.0, 0.0, 0.0
    atom1.charge = 1.0
    atom1.type = 0

    atom2 = pygcmc.MCAtom()
    atom2.x, atom2.y, atom2.z = float(distance), 0.0, 0.0
    atom2.charge = -1.0
    atom2.type = 0

    state.atoms = [atom1, atom2]
    state.activeAtomCount = 2

    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 1
    res1.active = True

    res2 = pygcmc.MCResidue()
    res2.atomStart = 1
    res2.atomCount = 1
    res2.active = True

    state.residues = [res1, res2]
    state.activeResidueCount = 2

    state.forcefield.numTotalTypes = 1
    state.forcefield.ljSigma = [0.3]
    state.forcefield.ljEps = [0.5]

    return state


def test_total_energy_unique_pairs_direct_matches_half_residue_sum():
    state = _make_two_residue_system(distance=0.3)

    energy_unique_from_compute = pygcmc.computeSystemEnergyCutoff(state)

    residue_sum_double = sum(
        (res.energy_vdw + res.energy_elec) for res in state.residues if res.active
    )
    expected_unique = 0.5 * residue_sum_double

    assert energy_unique_from_compute == pytest.approx(expected_unique, rel=0, abs=1e-10)
    assert residue_sum_double == pytest.approx(2.0 * energy_unique_from_compute, rel=0, abs=1e-10)

    energy_unique_cpp = pygcmc.getTotalEnergyUniquePairs(state, pygcmc.EnergyMethod.DIRECT)
    assert energy_unique_cpp == pytest.approx(expected_unique, rel=0, abs=1e-10)


def test_total_energy_unique_pairs_ewald_uses_ewald_total_plus_half_vdw():
    state = _make_two_residue_system(distance=0.6)
    pygcmc.initializeEwaldParameters(state.info.cutoff, state.info.box)

    elec_total, vdw_sum, _ = pygcmc.computeSystemEnergyEwald(state)

    # Ewald bindings currently return VdW as a double-counted residue sum.
    expected_unique = elec_total + 0.5 * vdw_sum

    energy_unique_cpp = pygcmc.getTotalEnergyUniquePairs(state, pygcmc.EnergyMethod.EWALD)
    assert energy_unique_cpp == pytest.approx(expected_unique, rel=0, abs=1e-8)

