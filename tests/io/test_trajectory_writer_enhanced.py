#!/usr/bin/env python
"""
Enhanced tests for trajectory/data writers
Focus on correctness of counts, units, and select numeric values.
"""

import os
import sys
import tempfile
import math
import pytest

# Add build path for pygcmc module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build'))

try:
    import pygcmc
    PYGCMC_AVAILABLE = True
except ImportError:
    PYGCMC_AVAILABLE = False
    pygcmc = None


def _build_min_state():
    state = pygcmc.MCState()
    state.info.box = (5.0, 5.0, 5.0)
    state.info.setTemperature(300.0)

    # Minimal forcefield
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljSigma = [0.315]
    ff.ljEps = [0.0]
    state.forcefield = ff

    # Two atoms only
    a1 = pygcmc.MCAtom(); a1.type = 0; a1.x=a1.y=a1.z=1.0
    a2 = pygcmc.MCAtom(); a2.type = 0; a2.x=a2.y=a2.z=2.0
    state.atoms = [a1, a2]

    # One residue declaring more atoms than exist to test XYZ count robustness
    res = pygcmc.MCResidue()
    res.active = True
    res.atomStart = 0
    res.atomCount = 10  # intentionally larger than actual
    res.atoms = [a1, a2]
    res.energy_vdw = -1.0
    res.energy_elec = -2.0
    state.residues = [res]
    state.activeAtomCount = 2
    state.activeResidueCount = 1

    return state


@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="PyGCMC not available")
def test_xyz_counts_and_units():
    state = _build_min_state()
    with tempfile.NamedTemporaryFile(suffix='.xyz', delete=False) as f:
        filename = f.name
    try:
        w = pygcmc.TrajectoryWriter(filename, pygcmc.TrajectoryWriter.Format.XYZ)
        w.write_frame(state, 0)
        w.close()
        with open(filename, 'r') as fh:
            lines = [l.rstrip('\n') for l in fh]
        # First line must equal actual number of valid atoms (2)
        assert lines[0].strip() == '2'
        # There should be exactly 2 atom lines after the header and comment lines
        atom_lines = lines[2:]
        assert len(atom_lines) == 2
        # Coordinates converted nm->Å (1.0->10.0, 2.0->20.0)
        coords1 = atom_lines[0].split()[1:]
        assert all(abs(float(v) - 10.0) < 1e-6 for v in coords1)
        coords2 = atom_lines[1].split()[1:]
        assert all(abs(float(v) - 20.0) < 1e-6 for v in coords2)
    finally:
        if os.path.exists(filename):
            os.remove(filename)


@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="PyGCMC not available")
def test_dat_energy_and_volume():
    state = _build_min_state()
    with tempfile.NamedTemporaryFile(suffix='.dat', delete=False) as f:
        filename = f.name
    try:
        w = pygcmc.TrajectoryWriter(filename, pygcmc.TrajectoryWriter.Format.DAT)
        w.write_frame(state, 0)
        w.close()
        with open(filename, 'r') as fh:
            parts = fh.readline().split()
        # frame, n_molecules, energy, volume
        assert len(parts) == 4
        assert int(parts[1]) == 1
        assert abs(float(parts[2]) - (-3.0)) < 1e-9
        assert abs(float(parts[3]) - 125.0) < 1e-9
    finally:
        if os.path.exists(filename):
            os.remove(filename)


@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="PyGCMC not available")
def test_top_temperature_line():
    state = _build_min_state()
    with tempfile.NamedTemporaryFile(suffix='.top', delete=False) as f:
        filename = f.name
    try:
        w = pygcmc.TrajectoryWriter(filename, pygcmc.TrajectoryWriter.Format.TOP)
        w.write_topology(state)
        w.close()
        with open(filename, 'r') as fh:
            content = fh.read()
        # Extract temperature value
        for line in content.splitlines():
            if line.startswith('# Temperature:'):
                temp_val = float(line.split(':', 1)[1].split()[0])
                assert abs(temp_val - 300.0) < 1e-6
                break
        else:
            pytest.fail('Temperature line not found')
    finally:
        if os.path.exists(filename):
            os.remove(filename)

