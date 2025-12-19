"""Regression for strict GROMACS pairtypes (1-4 only) in direct cutoff energy."""

import math

import pytest
import pygcmc


def _write_text(path, content):
    path.write_text(content.strip() + "\n")


def test_pairtypes_14_override_applies_in_direct_cutoff(tmp_path):
    work = tmp_path / "pairtypes_14_direct"
    work.mkdir(parents=True, exist_ok=True)

    pdb_path = work / "mol.pdb"
    _write_text(
        pdb_path,
        """
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1
ATOM      1  A1  MOL A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  B1  MOL A   1       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  B2  MOL A   1       2.000   0.000   0.000  1.00  0.00           C
ATOM      4  A2  MOL A   1       3.000   0.000   0.000  1.00  0.00           C
END
""",
    )

    top_path = work / "mol.top"
    _write_text(
        top_path,
        """
[ moleculetype ]
MOL  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   A   1   MOL  A1  1  0.000  12.011
2   B   1   MOL  B1  1  0.000  12.011
3   B   1   MOL  B2  1  0.000  12.011
4   A   1   MOL  A2  1  0.000  12.011

[ bonds ]
1 2 1
2 3 1
3 4 1

[ dihedrals ]
1 2 3 4 1

[ system ]
Pairtypes14

[ molecules ]
MOL 1
""",
    )

    structure = pygcmc.PDBParser.parse_file(str(pdb_path))
    topology = pygcmc.TOPParser.parse_file(str(top_path))
    molecular = pygcmc.MolecularSystem().combine(structure, topology)

    mc_system = pygcmc.MonteCarloSystem()
    info = pygcmc.MCInfo()
    info.max_residues = 10
    info.max_atoms = 10
    mc_system.initialize(info)
    mc_system.initialize_from_molecular(molecular)

    state = mc_system.get_state_mutable()
    state.info.cutoff = 2.0

    assert state.isPair14(0, 3)

    type_a = state.atomTypes.get_or_add_type("A")
    type_b = state.atomTypes.get_or_add_type("B")
    num_types = max(type_a, type_b) + 1

    state.forcefield.numTotalTypes = num_types
    state.forcefield.numMovementTypes = num_types
    state.forcefield.ljSigma = [0.0] * (num_types * num_types)
    state.forcefield.ljEps = [0.0] * (num_types * num_types)
    state.forcefield.ljMatrixInitialized = True

    sigma_nm = 0.28
    eps_kj = 1.2
    state.forcefield.setPairtype14(type_a, type_a, sigma_nm, eps_kj)

    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 0
    movement_info.activeCount = 1
    movement_info.totalCount = 1
    movement_info.resName = "MOL"
    state.movementResidues = [movement_info]

    pygcmc.computeMovementEnergyCutoff(state)

    atom0 = state.atoms[0]
    atom3 = state.atoms[3]
    dx = atom0.x - atom3.x
    dy = atom0.y - atom3.y
    dz = atom0.z - atom3.z
    r_nm = math.sqrt(dx * dx + dy * dy + dz * dz)
    sr = sigma_nm / r_nm
    expected = 4.0 * eps_kj * (sr ** 12 - sr ** 6)

    assert state.residues[0].energy_vdw == pytest.approx(expected, rel=1e-6, abs=1e-6)

    state.clearPair14()
    pygcmc.computeMovementEnergyCutoff(state)
    assert state.residues[0].energy_vdw == pytest.approx(0.0, abs=1e-8)
