"""
CBMC boundary/robustness regressions on a multi-atom fragment (no stdout/stderr matching).

Contracts:
- cbmcTrials=1 must reduce to an unbiased choice: rosenbluthWeight == 1
- cbmcTrials>1 must expose finite trial energies and close the Rosenbluth terms
- inserted fragment geometry must remain rigid (pairwise distances invariant vs template)
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from .inp_units_compat_helpers import _first_accept_record, _run_gcmc_cpu, _write_inp


def _write_text(path: Path, content: str) -> None:
    path.write_text(content.strip() + "\n")


def _parse_first_residue_coords(pdb_path: Path, *, resname: str) -> list[tuple[str, float, float, float]]:
    want = resname.strip().upper()
    atoms: list[tuple[str, float, float, float]] = []
    first_resid = None
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[17:20].strip().upper() != want:
            continue
        resid = int(line[22:26])
        atom_name = line[12:16].strip()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        if first_resid is None:
            first_resid = resid
        if resid != first_resid:
            break
        atoms.append((atom_name, x, y, z))
    return atoms


def _pairwise_distances(coords: list[tuple[str, float, float, float]]) -> list[float]:
    pts = [(x, y, z) for _, x, y, z in coords]
    dists: list[float] = []
    for i in range(len(pts)):
        xi, yi, zi = pts[i]
        for j in range(i + 1, len(pts)):
            xj, yj, zj = pts[j]
            dx = xj - xi
            dy = yj - yi
            dz = zj - zi
            dists.append(math.sqrt(dx * dx + dy * dy + dz * dz))
    dists.sort()
    return dists


def _expected_cbmc_log_w_over_k(*, energies_kj_mol: list[float], beta: float) -> float:
    if not energies_kj_mol:
        raise AssertionError("No CBMC trial energies provided")
    min_e = min(energies_kj_mol)
    sum_scaled = 0.0
    for e in energies_kj_mol:
        sum_scaled += math.exp(-beta * (e - min_e))
    avg_scaled = sum_scaled / float(len(energies_kj_mol))
    return math.log(max(avg_scaled, 1e-30)) - beta * min_e


def _clamp(x: float, lo: float, hi: float) -> float:
    return min(max(x, lo), hi)


def _write_multiatom_cbmc_system(work: Path) -> dict[str, Path]:
    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
ATOM      1  C   MOL A   1      15.000  15.000  15.000  1.00  0.00           C
END
""",
    )

    top = work / "sys.top"
    _write_text(
        top,
        """
[ defaults ]
1 2 yes 0.5 0.8333

[ moleculetype ]
MOL  2

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   C     1      MOL      C     1     0.500   12.011

[ system ]
Minimal

[ molecules ]
MOL  1
""",
    )

    par = work / "par.itp"
    _write_text(
        par,
        """
[ defaults ]
1 2 yes 0.5 0.8333

[ atomtypes ]
; name  at.num  mass   charge  ptype  sigma   epsilon
C       0       12.011 0.000   A      0.320   0.000
X       0       1.000  0.000   A      0.280   0.000
""",
    )

    frag = work / "frag.itp"
    _write_text(
        frag,
        """
[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   X     1      FRG      X1    1     -0.100  1.000
2   X     1      FRG      X2    1     -0.100  1.000
3   X     1      FRG      X3    1     -0.100  1.000
4   X     1      FRG      X4    1     -0.100  1.000
5   X     1      FRG      X5    1     -0.100  1.000
6   X     1      FRG      X6    1     -0.100  1.000
""",
    )

    frag_pdb = work / "frag.pdb"
    _write_text(
        frag_pdb,
        """
ATOM      1  X1  FRG A   1      -2.500   0.000   0.000  1.00  0.00           C
ATOM      2  X2  FRG A   1      -1.500   0.000   0.000  1.00  0.00           C
ATOM      3  X3  FRG A   1      -0.500   0.000   0.000  1.00  0.00           C
ATOM      4  X4  FRG A   1       0.500   0.000   0.000  1.00  0.00           C
ATOM      5  X5  FRG A   1       1.500   0.000   0.000  1.00  0.00           C
ATOM      6  X6  FRG A   1       2.500   0.000   0.000  1.00  0.00           C
END
""",
    )

    return {"pdb": pdb, "top": top, "par": par, "frag": frag, "frag_pdb": frag_pdb}


def test_cbmc_trials_one_has_unity_rosenbluth_and_preserves_geometry(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "cbmc_multiatom" / "k1"
    work.mkdir(parents=True, exist_ok=True)

    files = _write_multiatom_cbmc_system(work)

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "run.inp"
    _write_inp(
        inp,
        f"""
random_seed:12345
par:{files["par"]}
fragitp:{files["frag"]}
fragname:FRG
fragconc:55.0
fragmuex:0.0

pdb:{files["pdb"]}
top:{files["top"]}
box_size:30.0 30.0 30.0
gcmc_region:box 10.0 10.0 10.0 20.0 20.0 20.0
cutoff:12.0
temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1

use_cavity_bias:no
use_conf_bias:yes
fragconf:1

attempt_prob_ins:1.0
attempt_prob_del:0.0
attempt_prob_trn:0.0
attempt_prob_rot:0.0
""",
    )

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work,
        inp=inp,
        out_prefix=out_prefix,
        extra_args=["--dump-accept", str(accept_log)],
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="insertion", species="FRG")
    assert int(rec.get("cbmcTrials", 0)) == 1

    beta = float(rec["beta"])
    selected = float(rec["cbmcSelectedEnergy"])
    energies = [float(x) for x in rec.get("cbmcTrialEnergies", [])]
    # Implementation detail: for K=1 we allow omitting the trial-energy vector,
    # but the selected-energy / logWOverK / rosenbluthWeight semantics must still close.
    assert len(energies) in (0, 1)
    if energies:
        assert energies[0] == pytest.approx(selected, abs=1e-12)

    # K=1 => logWOverK = -βE, and rosenbluthWeight == exp(logWOverK + βE) == 1 exactly.
    assert float(rec["cbmcLogWOverK"]) == pytest.approx(-beta * selected, abs=1e-12)
    assert float(rec["rosenbluthWeight"]) == pytest.approx(1.0, abs=1e-12)

    out_pdb = Path(f"{out_prefix}_final.pdb")
    assert out_pdb.exists()

    inserted = _parse_first_residue_coords(out_pdb, resname="FRG")
    template = _parse_first_residue_coords(files["frag_pdb"], resname="FRG")
    assert len(inserted) == len(template) == 6

    inserted_d = _pairwise_distances(inserted)
    template_d = _pairwise_distances(template)
    assert len(inserted_d) == len(template_d)
    for i in range(len(template_d)):
        assert inserted_d[i] == pytest.approx(template_d[i], abs=0.03)


def test_cbmc_trials_multiatom_closes_rosenbluth_terms(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "cbmc_multiatom" / "k8"
    work.mkdir(parents=True, exist_ok=True)

    files = _write_multiatom_cbmc_system(work)

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "run.inp"
    _write_inp(
        inp,
        f"""
random_seed:54321
par:{files["par"]}
fragitp:{files["frag"]}
fragname:FRG
fragconc:1.0
fragmuex:0.0

pdb:{files["pdb"]}
top:{files["top"]}
box_size:30.0 30.0 30.0
gcmc_region:box 0.0 0.0 0.0 30.0 30.0 30.0
cutoff:15.0
temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1

use_cavity_bias:no
use_conf_bias:yes
fragconf:8

attempt_prob_ins:1.0
attempt_prob_del:0.0
attempt_prob_trn:0.0
attempt_prob_rot:0.0
""",
    )

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work,
        inp=inp,
        out_prefix=out_prefix,
        extra_args=["--dump-accept", str(accept_log)],
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="insertion", species="FRG")
    cbmc_trials = int(rec.get("cbmcTrials", 0))
    assert 1 < cbmc_trials <= 8

    energies = [float(x) for x in rec.get("cbmcTrialEnergies", [])]
    assert len(energies) == cbmc_trials
    assert all(math.isfinite(e) for e in energies)
    assert max(energies) - min(energies) > 1e-6, "CBMC trial energies unexpectedly degenerate"

    beta = float(rec["beta"])
    expected_log_w_over_k = _expected_cbmc_log_w_over_k(energies_kj_mol=energies, beta=beta)
    assert float(rec["cbmcLogWOverK"]) == pytest.approx(expected_log_w_over_k, abs=1e-12)

    log_rosen = expected_log_w_over_k + beta * float(rec["cbmcSelectedEnergy"])
    log_rosen = _clamp(log_rosen, math.log(1e-30), 700.0)
    expected_rosen = math.exp(log_rosen)
    assert float(rec["rosenbluthWeight"]) == pytest.approx(expected_rosen, rel=1e-12, abs=1e-12)
