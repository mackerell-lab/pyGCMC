"""
Quantify and lock down the C/D/E energy-method differences for gcmc_cpu.

Modes (gcmc_cpu INP `energy_method:`):
- C: pgp_full        (no mesh-self; "more physical"/translation-invariant baseline)
- D: pme             (full PME; includes mesh-self from the discrete grid)
- E: pgp_full_pme    (PGP(full) + mesh-self; must match D per-move ΔU)

We keep these tests end-to-end and file-driven:
- Run `gcmc_cpu` with `--dump-accept` JSONL.
- Compare per-move `deltaU` sequences and derived statistics (no stdout/stderr matching).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from cpu_cli.inp_units_compat_helpers import _run_gcmc_cpu, _write_inp


def _write_text(path: Path, content: str) -> None:
    path.write_text(content.strip() + "\n")


def _load_translation_delta_us(path: Path, *, species: str) -> list[float]:
    want_species = species.strip().upper()
    deltas: list[float] = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        if str(rec.get("move", "")).strip().lower() != "translation":
            continue
        if str(rec.get("species", "")).strip().upper() != want_species:
            continue
        deltas.append(float(rec["deltaU"]))
    return deltas


def _mean_abs(values: list[float]) -> float:
    if not values:
        return 0.0
    return sum(abs(v) for v in values) / float(len(values))


def _setup_neutral_dipole_system(
    work: Path,
    *,
    dipole_charge: float,
    include_fixed_host: bool,
) -> tuple[Path, Path, Path, Path]:
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    top = work / "sys.top"
    par = work / "par.itp"
    frag = work / "frag.itp"
    # Coordinate template for the fragment (required for multi-atom moves).
    # `SimulationInputBuilder` will auto-detect monomer PDB as:
    #   monomerdir/<frag>.pdb or <fragitp_dir>/<frag>.pdb
    # and `FragmentLibrary` will center it to COM and store per-atom relative coords.
    frag_pdb = work / "frg.pdb"

    host_pdb = ""
    host_top = ""
    host_molecules = ""
    host_types = ""
    if include_fixed_host:
        host_pdb = """
ATOM      1  HP  HST A   1       5.000   5.000   5.000  1.00  0.00           C
ATOM      2  HM  HSN A   2      45.000  45.000  45.000  1.00  0.00           C
"""
        host_top = """
[ moleculetype ]
HST  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   HP    1      HST      HP    1     1.000  12.011

[ moleculetype ]
HSN  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   HM    1      HSN      HM    1    -1.000  12.011
"""
        host_molecules = """
HST  1
HSN  1
"""
        host_types = """
HP      0       12.011 0.000   A      0.300   0.000
HM      0       12.011 0.000   A      0.300   0.000
"""

    _write_text(
        pdb,
        f"""
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1
{host_pdb.rstrip()}
ATOM      3  QP  FRG B   1      25.000  25.000  25.000  1.00  0.00           C
ATOM      4  QM  FRG B   1      26.000  25.000  25.000  1.00  0.00           C
END
""",
    )

    q = float(dipole_charge)
    _write_text(
        top,
        f"""
[ defaults ]
1 2 yes 1.0 1.0

{host_top.rstrip()}

[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   QP    1      FRG      QP    1     {q:.6f}  12.011
2   QM    1      FRG      QM    1    {-q:.6f}  12.011

[ system ]
NeutralDipole

[ molecules ]
{host_molecules.rstrip()}
FRG  1
""",
    )

    _write_text(
        par,
        f"""
[ defaults ]
1 2 yes 1.0 1.0

[ atomtypes ]
; name  at.num  mass   charge  ptype  sigma   epsilon
QP      0       12.011 0.000   A      0.300   0.000
QM      0       12.011 0.000   A      0.300   0.000
{host_types.rstrip()}
""",
    )

    _write_text(
        frag,
        f"""
[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   QP    1      FRG      QP    1     {q:.6f}  12.011
2   QM    1      FRG      QM    1    {-q:.6f}  12.011
""",
    )

    # Provide a minimal FRG coordinate PDB so translation/rotation preserve geometry.
    _write_text(
        frag_pdb,
        """
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1
ATOM      1  QP  FRG A   1       0.000   0.000   0.000  1.00  0.00           Q
ATOM      2  QM  FRG A   1       1.000   0.000   0.000  1.00  0.00           Q
END
""",
    )

    return pdb, top, par, frag


def _run_translation_delta_us(
    gcmc_cpu: str,
    *,
    work: Path,
    method: str,
    pdb: Path,
    top: Path,
    par: Path,
    frag: Path,
    seed: int,
    mcsteps: int,
) -> list[float]:
    run_dir = work / method
    run_dir.mkdir(parents=True, exist_ok=True)

    out_prefix = run_dir / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = run_dir / "out" / "accept.jsonl"

    inp = run_dir / "run.inp"
    _write_inp(
        inp,
        f"""
random_seed:{seed}
energy_method:{method}
par:{par}
fragitp:{frag}
fragname:FRG
fragconc:0.0
fragmuex:0.0

pdb:{pdb}
top:{top}
box_size:50.0 50.0 50.0
cutoff:6.0
temperature:1e30
moves_per_step:1
mcsteps:{mcsteps}
nprint:{mcsteps}
mc_move_prob:0 0 1 0
max_translation:5.0
max_rotation:0.0

use_cavity_bias:no
use_conf_bias:no
""",
    )

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=run_dir,
        inp=inp,
        out_prefix=out_prefix,
        extra_args=[
            "--dump-accept",
            str(accept_log),
            "--no-stats",
            "--print-freq",
            "1000000000",
            "--traj-freq",
            "1000000000",
            "--checkpoint-freq",
            "0",
        ],
        timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    deltas = _load_translation_delta_us(accept_log, species="FRG")
    assert len(deltas) == mcsteps
    return deltas


def test_cde_mesh_self_gap_is_background_independent_for_neutral_dipole_translation(gcmc_cpu, temp_dir):
    """
    For a rigid neutral dipole:
    - The per-move gap (D - C) is the mesh-self term Δ(0.5 ρᵀKρ) and must be independent
      of any fixed background, given identical move proposals.

    Notes:
    - Mode E is designed to match the discrete PME target (Mode D) by adding the PME mesh-self term.
      For multi-atom fragments there can still be a tiny residual (D-E) from intra-residue real-space
      electrostatics that Mode C/E intentionally do not treat as part of the interaction ΔU.
      We therefore lock down that E is *much closer* to D than C is, while keeping the hard,
      theory-backed (D-C) contracts strict.
    - Strict per-move D==E is covered separately in `cpu_cli/pgp_pme_equivalence.py` using a 1-atom
      charged fragment, which removes intra-residue pair contributions.
    """
    base = Path(temp_dir) / "pgp_cde_quantification" / "neutral_dipole_background_independence"
    base.mkdir(parents=True, exist_ok=True)

    mcsteps = 25
    seed = 424242
    dipole_charge = 2.0  # amplify mesh-self effect while staying neutral overall

    vac = base / "vacuum"
    vac_files = _setup_neutral_dipole_system(vac, dipole_charge=dipole_charge, include_fixed_host=False)
    host = base / "with_host"
    host_files = _setup_neutral_dipole_system(host, dipole_charge=dipole_charge, include_fixed_host=True)

    def _run_all(system_dir: Path, files: tuple[Path, Path, Path, Path]) -> tuple[list[float], list[float], list[float]]:
        pdb, top, par, frag = files
        deltas_c = _run_translation_delta_us(
            gcmc_cpu,
            work=system_dir,
            method="pgp_full",
            pdb=pdb,
            top=top,
            par=par,
            frag=frag,
            seed=seed,
            mcsteps=mcsteps,
        )
        deltas_d = _run_translation_delta_us(
            gcmc_cpu,
            work=system_dir,
            method="pme",
            pdb=pdb,
            top=top,
            par=par,
            frag=frag,
            seed=seed,
            mcsteps=mcsteps,
        )
        deltas_e = _run_translation_delta_us(
            gcmc_cpu,
            work=system_dir,
            method="pgp_full_pme",
            pdb=pdb,
            top=top,
            par=par,
            frag=frag,
            seed=seed,
            mcsteps=mcsteps,
        )
        return deltas_c, deltas_d, deltas_e

    vac_c, vac_d, vac_e = _run_all(vac, vac_files)
    host_c, host_d, host_e = _run_all(host, host_files)

    def _max_abs_diff(a: list[float], b: list[float]) -> float:
        return max(abs(x - y) for x, y in zip(a, b))

    vac_gap = [d - c for d, c in zip(vac_d, vac_c)]
    host_gap = [d - c for d, c in zip(host_d, host_c)]

    assert host_gap == pytest.approx(vac_gap, rel=1e-10, abs=1e-6)
    assert max(abs(x) for x in vac_gap) > 1e-4

    # E should dramatically reduce the PME mesh-self mismatch relative to C.
    # (Exact D==E is validated elsewhere in a single-atom charged fragment case.)
    vac_de = _max_abs_diff(vac_d, vac_e)
    host_de = _max_abs_diff(host_d, host_e)
    vac_dc = _max_abs_diff(vac_d, vac_c)
    host_dc = _max_abs_diff(host_d, host_c)

    assert vac_de < 0.1 * vac_dc
    assert host_de < 0.1 * host_dc
    assert vac_de < 0.02
    assert host_de < 0.02


def test_mesh_self_gap_scales_quadratically_with_charge_for_neutral_dipole(gcmc_cpu, temp_dir):
    """
    Theory: mesh-self term is quadratic in charge density: E_mesh_self ∝ ρᵀKρ.

    So if we scale all dipole charges by factor s, the (D - C) mesh-self gap
    should scale by s^2 (here s=2 => ~4x).
    """
    base = Path(temp_dir) / "pgp_cde_quantification" / "neutral_dipole_charge_scaling"
    base.mkdir(parents=True, exist_ok=True)

    mcsteps = 25
    seed = 13579

    sys_q1 = base / "q1"
    files_q1 = _setup_neutral_dipole_system(sys_q1, dipole_charge=1.0, include_fixed_host=False)
    sys_q2 = base / "q2"
    files_q2 = _setup_neutral_dipole_system(sys_q2, dipole_charge=2.0, include_fixed_host=False)

    def _gap_mean_abs(system_dir: Path, files: tuple[Path, Path, Path, Path]) -> float:
        pdb, top, par, frag = files
        deltas_c = _run_translation_delta_us(
            gcmc_cpu,
            work=system_dir,
            method="pgp_full",
            pdb=pdb,
            top=top,
            par=par,
            frag=frag,
            seed=seed,
            mcsteps=mcsteps,
        )
        deltas_d = _run_translation_delta_us(
            gcmc_cpu,
            work=system_dir,
            method="pme",
            pdb=pdb,
            top=top,
            par=par,
            frag=frag,
            seed=seed,
            mcsteps=mcsteps,
        )
        gap = [d - c for d, c in zip(deltas_d, deltas_c)]
        return _mean_abs(gap)

    gap1 = _gap_mean_abs(sys_q1, files_q1)
    gap2 = _gap_mean_abs(sys_q2, files_q2)

    assert gap1 > 1e-6
    assert gap2 / gap1 == pytest.approx(4.0, rel=0.05, abs=0.0)
