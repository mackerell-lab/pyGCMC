"""
Test hard sphere exclusion in systems with repulsive interactions
"""
import math
import pygcmc


def _lj_energy(eps, sigma, r):
    sr = sigma / r
    sr6 = sr ** 6
    sr12 = sr6 * sr6
    return 4.0 * eps * (sr12 - sr6)


def test_hard_sphere_repulsion_documentation():
    """Repulsive LJ should strongly penalize close inter-molecular distances."""
    state = pygcmc.MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 2.0

    eps = 5.0
    sigma = 0.4

    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.maxTypes = 1
    ff.ljEps = [eps]
    ff.ljSigma = [sigma]
    state.forcefield = ff

    atom1 = pygcmc.MCAtom()
    atom1.x = 0.0
    atom1.y = 0.0
    atom1.z = 0.0
    atom1.charge = 0.0
    atom1.type = 0

    atom2 = pygcmc.MCAtom()
    atom2.y = 0.0
    atom2.z = 0.0
    atom2.charge = 0.0
    atom2.type = 0

    state.atoms = [atom1, atom2]
    state.activeAtomCount = 2

    res1 = pygcmc.MCResidue()
    res1.active = True
    res1.type = 0
    res1.atomStart = 0
    res1.atomCount = 1

    res2 = pygcmc.MCResidue()
    res2.active = True
    res2.type = 0
    res2.atomStart = 1
    res2.atomCount = 1

    state.residues = [res1, res2]
    state.activeResidueCount = 2

    r_close = 0.2
    r_far = 0.6

    state.atoms[1].x = r_close
    energy_close = pygcmc.computeSystemEnergyCutoff(state)

    state.atoms[1].x = r_far
    energy_far = pygcmc.computeSystemEnergyCutoff(state)

    expected_close = _lj_energy(eps, sigma, r_close)
    expected_far = _lj_energy(eps, sigma, r_far)

    assert energy_close > energy_far
    assert math.isfinite(energy_close)
    assert math.isfinite(energy_far)
    assert abs(energy_close - expected_close) / expected_close < 1e-5
    assert abs(energy_far - expected_far) / max(abs(expected_far), 1e-6) < 1e-5


if __name__ == "__main__":
    test_hard_sphere_repulsion_documentation()
    print("Test completed")
