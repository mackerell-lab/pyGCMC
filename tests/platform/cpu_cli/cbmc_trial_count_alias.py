"""
CBMC trial-count compatibility: num_conf_bias_trial (gcmc_opencl-style) must drive cbmcTrials
when fragconf is not provided.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from cpu_cli.inp_units_compat_helpers import _first_accept_record, _run_gcmc_cpu, _write_inp


def _write_text(path: Path, content: str) -> None:
    path.write_text(content.strip() + "\n")


def test_num_conf_bias_trial_sets_cbmc_trials_when_fragconf_missing(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "cbmc_trials_num_conf_bias_trial"
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
    cbmc_trials = 7
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
cutoff:6.0
temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:1 0 0 0

use_cavity_bias:no
use_conf_bias:yes
num_conf_bias_trial:{cbmc_trials}
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
    assert accept_log.exists()

    rec = _first_accept_record(accept_log, move="insertion", species="FRG")
    assert int(rec["cbmcTrials"]) == cbmc_trials
    energies = rec.get("cbmcTrialEnergies", [])
    assert isinstance(energies, list)
    assert len(energies) == cbmc_trials

