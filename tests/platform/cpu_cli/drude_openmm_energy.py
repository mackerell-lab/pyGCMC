"""
gcmc_cpu Drude energy closure + OpenMM cross-check.

This is a strict, file-driven regression:
- Runs the compiled gcmc_cpu CLI on a minimal Drude PSF system.
- Reads energy from *_statistics.dat (no stdout/stderr parsing).
- Verifies total energy closes against an analytic decomposition (Coulomb + Drude spring).
- Cross-checks the same configuration against OpenMM (Reference platform) energy components.
"""

from __future__ import annotations

import math
import subprocess
import warnings
import importlib.util
from pathlib import Path

import pytest

# OpenMM (SWIG) on Python 3.12 emits noisy DeprecationWarnings like:
# "builtin type SwigPyPacked has no __module__ attribute".
# Filter locally so this test remains quiet without relying on pytest.ini.
warnings.filterwarnings(
    "ignore",
    category=DeprecationWarning,
    message=r"builtin type SwigPyPacked has no __module__ attribute",
)
warnings.filterwarnings(
    "ignore",
    category=DeprecationWarning,
    message=r"builtin type SwigPyObject has no __module__ attribute",
)
warnings.filterwarnings(
    "ignore",
    category=DeprecationWarning,
    message=r"builtin type swigvarlink has no __module__ attribute",
)

HAS_OPENMM = importlib.util.find_spec("openmm") is not None
openmm = None
unit = None


COULOMB_KJ_NM_PER_MOL_E2 = 138.935456  # must match DrudeConstants.ONE_4PI_EPS0


def _require_openmm():
    global openmm, unit
    if openmm is None or unit is None:
        import openmm as _openmm
        import openmm.unit as _unit

        openmm = _openmm
        unit = _unit
    return openmm, unit


def _write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content.strip() + "\n")


def _pdb_atom_xyz(pdb_path: Path, *, resname: str, atom_name: str, resid: int) -> tuple[float, float, float]:
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        parts = line.split()
        if len(parts) < 9:
            continue
        if parts[2].upper() != atom_name.upper():
            continue
        if parts[3].upper() != resname.upper():
            continue
        if int(parts[5]) != resid:
            continue
        return (float(parts[6]), float(parts[7]), float(parts[8]))  # Å
    raise AssertionError(f"Atom not found in {pdb_path}: {resname}{resid} {atom_name}")


def _dist_nm(a_a: tuple[float, float, float], b_a: tuple[float, float, float]) -> float:
    dx = (a_a[0] - b_a[0]) * 0.1
    dy = (a_a[1] - b_a[1]) * 0.1
    dz = (a_a[2] - b_a[2]) * 0.1
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _read_last_stats_energy_kj(stats_path: Path) -> float:
    lines = [ln for ln in stats_path.read_text().splitlines() if ln and not ln.startswith("#")]
    assert lines, f"No data lines in statistics file: {stats_path}"
    last = lines[-1].split()
    assert len(last) >= 2, f"Malformed statistics line: {lines[-1]}"
    return float(last[1])


def _setup_minimal_drude_external_field_system(work: Path) -> tuple[Path, Path, Path, Path, Path]:
    """
    Minimal 2-residue system:
      RES: parent + Drude (bonded)
      FIX: external point charge (fixed framework)

    External units: Å / kcal decks are supported; PSF/PDB use Å for coordinates.
    """
    psf = work / "system.psf"
    pdb = work / "system.pdb"
    prm = work / "ff.prm"
    frag_itp = work / "res.itp"
    monomerdir = work / "monomerdir"
    monomerdir.mkdir(parents=True, exist_ok=True)
    frag_pdb = monomerdir / "RES.pdb"

    parent_type = "P"
    drude_type = "D"
    ext_type = "Q"

    alpha_a3 = 1.0  # Å^3 -> 0.001 nm^3
    thole = 2.6
    drude_dx_a = 0.005  # Å, tiny non-zero initial offset

    _write_text(
        psf,
        f"""
PSF EXT DRUDE

       1 !NTITLE
REMARKS minimal Drude energy/openmm regression

       3 !NATOM
       1 SYS  1 RES  P1   {parent_type}   1.000000  12.0110  0  {thole:.3f}  {alpha_a3:.3f}
       2 SYS  1 RES  DP1  {drude_type}   -1.000000   0.4000  0  0.000  0.000
       3 SYS  2 FIX  Q1   {ext_type}      1.000000  12.0110  0  0.000  0.000

       1 !NBOND: bonds
       1       2
""",
    )

    cx, cy, cz = 25.0, 25.0, 25.0
    ext_dx_a = 3.0  # Å
    _write_text(
        pdb,
        f"""
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1
ATOM      1  P1  RES A   1    {cx:8.3f}{cy:8.3f}{cz:8.3f}  1.00  0.00           C
ATOM      2 DP1  RES A   1    {cx + drude_dx_a:8.3f}{cy:8.3f}{cz:8.3f}  1.00  0.00           D
ATOM      3  Q1  FIX A   2    {cx + ext_dx_a:8.3f}{cy:8.3f}{cz:8.3f}  1.00  0.00           Q
END
""",
    )

    _write_text(
        frag_pdb,
        f"""
ATOM      1  P1  RES A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2 DP1  RES A   1       {drude_dx_a:6.3f}   0.000   0.000  1.00  0.00           D
END
""",
    )

    _write_text(
        prm,
        f"""
* minimal PRM for Drude energy/openmm regression
*

NONBONDED nbxmod 5 atom cdiel -
cutnb 999.0 ctonnb 0.0 ctofnb 0.0 eps 1.0 e14fac 1.0 wmin 1.5
{parent_type:8s}  0   0.0   0.0
{drude_type:8s}   0   0.0   0.0
{ext_type:8s}     0   0.0   0.0

ALPHA
{parent_type:8s}  {alpha_a3:.3f}  {thole:.3f}
END
""",
    )

    _write_text(
        frag_itp,
        f"""
[ moleculetype ]
RES     3

[ atoms ]
; nr  type   resnr  residue  atom  cgnr  charge    mass
1   {parent_type}  1   RES  P1   1   1.000000  12.0110
2   {drude_type}   1   RES  DP1  1  -1.000000   0.4000

[ bonds ]
1  2  1
""",
    )

    return psf, pdb, prm, frag_itp, monomerdir


def _openmm_energy_components_kj(positions_a: list[tuple[float, float, float]]) -> tuple[float, float, float]:
    assert HAS_OPENMM
    openmm_mod, unit_mod = _require_openmm()
    positions_nm = [(x * 0.1, y * 0.1, z * 0.1) for (x, y, z) in positions_a]

    system = openmm_mod.System()
    system.addParticle(12.011 * unit_mod.dalton)  # parent
    system.addParticle(0.400 * unit_mod.dalton)  # drude
    system.addParticle(12.011 * unit_mod.dalton)  # external

    drude_force = openmm_mod.DrudeForce()
    drude_force.setForceGroup(1)
    drude_force.addParticle(1, 0, -1, -1, -1, -1.0, 0.001, 1.0, 1.0)  # q=-1, alpha=0.001 nm^3
    system.addForce(drude_force)

    nb = openmm_mod.NonbondedForce()
    nb.setForceGroup(2)
    nb.setNonbondedMethod(openmm_mod.NonbondedForce.NoCutoff)
    sigma = 0.1 * unit_mod.nanometer
    epsilon = 0.0 * unit_mod.kilojoule_per_mole
    nb.addParticle(1.0 * unit_mod.elementary_charge, sigma, epsilon)   # parent
    nb.addParticle(-1.0 * unit_mod.elementary_charge, sigma, epsilon)  # drude
    nb.addParticle(1.0 * unit_mod.elementary_charge, sigma, epsilon)   # external
    nb.addException(0, 1, 0.0, 1.0, 0.0)  # exclude parent-drude nonbonded (bonded in PSF)
    system.addForce(nb)

    platform = openmm_mod.Platform.getPlatformByName("Reference")
    integrator = openmm_mod.VerletIntegrator(0.001 * unit_mod.picoseconds)
    context = openmm_mod.Context(system, integrator, platform)
    context.setPositions(positions_nm * unit_mod.nanometer)

    total = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit_mod.kilojoule_per_mole)
    drude_e = context.getState(getEnergy=True, groups=1 << 1).getPotentialEnergy().value_in_unit(unit_mod.kilojoule_per_mole)
    nb_e = context.getState(getEnergy=True, groups=1 << 2).getPotentialEnergy().value_in_unit(unit_mod.kilojoule_per_mole)
    return total, nb_e, drude_e


@pytest.mark.skipif(not HAS_OPENMM, reason="OpenMM not available")
def test_gcmc_cpu_drude_energy_matches_openmm_and_analytic_closure(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "drude_openmm_energy"
    work.mkdir(parents=True, exist_ok=True)

    psf, pdb, prm, frag_itp, monomerdir = _setup_minimal_drude_external_field_system(work)

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    inp = work / "run.inp"
    _write_text(
        inp,
        f"""
energy_method:direct
top:{psf}
pdb:{pdb}
par:{prm}
monomerdir:{monomerdir}
fragitp:{frag_itp}
fragname:RES
fragconc:55.0
fragmuex:0.0

box_size:50.0 50.0 50.0
cutoff:12.0
temperature:300.0

moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:0 0 1 0
max_translation:0.0
max_rotation:0.0
""",
    )

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--seed", "123", "--stats-interval", "1"],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=40,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    out_pdb = Path(f"{out_prefix}_final.pdb")
    stats = Path(f"{out_prefix}_statistics.dat")
    assert out_pdb.exists()
    assert stats.exists()

    energy_gcmc = _read_last_stats_energy_kj(stats)

    parent = _pdb_atom_xyz(out_pdb, resname="RES", atom_name="P1", resid=1)
    drude = _pdb_atom_xyz(out_pdb, resname="RES", atom_name="DP1", resid=1)
    ext = _pdb_atom_xyz(out_pdb, resname="FIX", atom_name="Q1", resid=2)

    # Analytic closure: total = Coulomb(parent-ext + drude-ext) + harmonic spring.
    q_parent, q_drude, q_ext = 1.0, -1.0, 1.0
    alpha_nm3 = 0.001
    r_pe = _dist_nm(parent, ext)
    r_de = _dist_nm(drude, ext)
    disp_pd = _dist_nm(parent, drude)

    coulomb = COULOMB_KJ_NM_PER_MOL_E2 * q_parent * q_ext / r_pe + COULOMB_KJ_NM_PER_MOL_E2 * q_drude * q_ext / r_de
    k_spring = COULOMB_KJ_NM_PER_MOL_E2 * (q_drude * q_drude) / alpha_nm3
    spring = 0.5 * k_spring * disp_pd * disp_pd
    expected_total = coulomb + spring

    # Energy is computed from full-precision internal coordinates, while the analytic closure
    # is reconstructed from the (rounded) output PDB; allow a small absolute tolerance.
    assert energy_gcmc == pytest.approx(expected_total, abs=5e-3, rel=1e-6)

    omm_total, omm_nb, omm_drude = _openmm_energy_components_kj([parent, drude, ext])
    assert omm_total == pytest.approx(energy_gcmc, abs=5e-3, rel=1e-6)
    assert omm_nb == pytest.approx(coulomb, abs=5e-3, rel=1e-6)
    assert omm_drude == pytest.approx(spring, abs=5e-3, rel=1e-6)
