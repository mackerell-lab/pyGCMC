"""
CBMC diagnostics contract for gcmc_cpu (no stdout/stderr matching).

We export CBMC trial energies in --dump-accept and verify they close:
- cbmcLogWOverK against the logged energies
- rosenbluthWeight against exp(cbmcLogWOverK + beta * cbmcSelectedEnergy)
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from cpu_cli.inp_units_compat_helpers import _first_accept_record, _run_gcmc_cpu, _write_inp


def _write_text(path: Path, content: str) -> None:
    path.write_text(content.strip() + "\n")


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


def test_dump_accept_exposes_cbmc_trial_energies_and_closes_rosenbluth_terms(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "cbmc_trial_energies"
    work.mkdir(parents=True, exist_ok=True)

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
1   X     1      FRG      X     1     -0.250  1.000
""",
    )

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
random_seed:12345
par:{par}
fragitp:{frag}
fragname:FRG
fragconc:1.0
fragmuex:0.0

pdb:{pdb}
top:{top}
box_size:30.0 30.0 30.0
gcmc_region:box 14.0 14.0 14.0 16.0 16.0 16.0
cutoff:6.0
temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:1 0 0 0

use_cavity_bias:no
use_conf_bias:yes
fragconf:5
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
    assert int(rec.get("cbmcTrials", 0)) > 1

    energies = [float(x) for x in rec.get("cbmcTrialEnergies", [])]
    assert len(energies) == int(rec["cbmcTrials"])
    assert all(math.isfinite(e) for e in energies)

    beta = float(rec["beta"])
    expected_log_w_over_k = _expected_cbmc_log_w_over_k(energies_kj_mol=energies, beta=beta)
    assert float(rec["cbmcLogWOverK"]) == pytest.approx(expected_log_w_over_k, abs=1e-12)

    log_rosen = expected_log_w_over_k + beta * float(rec["cbmcSelectedEnergy"])
    log_rosen = _clamp(log_rosen, math.log(1e-30), 700.0)
    expected_rosen = math.exp(log_rosen)
    assert float(rec["rosenbluthWeight"]) == pytest.approx(expected_rosen, rel=1e-12, abs=1e-12)

