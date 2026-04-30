"""
Pybind CPU GCMC interface tests.

These tests exercise the C++ simulation through the installed/imported Python
extension, not through the standalone gcmc_cpu executable.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import pygcmc

from cpu_cli.inp_units_compat_helpers import _write_inp

PACKAGE_DIR = Path(__file__).resolve().parents[3] / "python"
if PACKAGE_DIR.exists() and str(PACKAGE_DIR) not in sys.path:
    sys.path.insert(0, str(PACKAGE_DIR))

from pygcmc_tools.paper_validation import PAPER_MODES, run_paper_modes


def _write_text(path: Path, content: str) -> None:
    path.write_text(content.strip() + "\n")


def _write_minimal_charged_translation_case(work: Path, *, mode: str = "pgp_full_pme", steps: int = 3) -> Path:
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1
ATOM      1  Q   HST A   1      10.000  20.000  20.000  1.00  0.00           C
ATOM      2  Q   FRG B   1      35.000  20.000  20.000  1.00  0.00           C
END
""",
    )

    top = work / "sys.top"
    _write_text(
        top,
        """
[ defaults ]
1 2 yes 1.0 1.0

[ moleculetype ]
HST  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   Q     1      HST      Q     1     1.000   12.011

[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   Q     1      FRG      Q     1    -1.000   12.011

[ system ]
Minimal

[ molecules ]
HST  1
FRG  1
""",
    )

    par = work / "par.itp"
    _write_text(
        par,
        """
[ defaults ]
1 2 yes 1.0 1.0

[ atomtypes ]
; name  at.num  mass   charge  ptype  sigma   epsilon
Q       0       12.011 0.000   A      0.300   0.000
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
1   Q     1      FRG      Q     1    -1.000  12.011
""",
    )

    inp = work / "run.inp"
    _write_inp(
        inp,
        f"""
random_seed:12345
energy_method:{mode}
par:{par}
fragitp:{frag}
fragname:FRG
fragconc:0.0
fragmuex:-1000.0

pdb:{pdb}
top:{top}
box_size:40.0 40.0 40.0
cutoff:6.0
temperature:1e30
moves_per_step:1
mcsteps:{steps}
nprint:{steps}
mc_move_prob:0 0 1 0
max_translation:1.0
max_rotation:0.0

use_cavity_bias:no
use_conf_bias:no
""",
    )
    return inp


def test_pybind_run_gcmc_cpu_executes_pgp_full_pme(temp_dir):
    work = Path(temp_dir) / "pybind_gcmc_cpu" / "pgp_full_pme"
    inp = _write_minimal_charged_translation_case(work, mode="pgp_full_pme", steps=3)

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "accept.jsonl"
    params_json = work / "out" / "params.json"

    config = pygcmc.GCMCCPUConfig()
    config.inputFile = str(inp)
    config.outputPrefix = str(out_prefix)
    config.randomSeed = 12345
    config.printFrequency = 1000000000
    config.trajectoryFrequency = 1000000000
    config.checkpointFrequency = 0
    config.enableStatistics = False

    result = pygcmc.run_gcmc_cpu(
        config,
        dumpAccept=str(accept_log),
        dumpParams=str(params_json),
    )

    assert result["returncode"] == 0
    assert result["initialized"] is True
    assert result["ran"] is True
    assert Path(f"{out_prefix}_final.pdb").exists()
    assert json.loads(params_json.read_text())["basic"]["energy_method"] == "pgp_full_pme"

    records = [json.loads(line) for line in accept_log.read_text().splitlines() if line.strip()]
    translation_deltas = [
        float(rec["deltaU"])
        for rec in records
        if str(rec.get("move", "")).lower() == "translation"
    ]
    assert len(translation_deltas) == 3
    assert all(math.isfinite(delta) for delta in translation_deltas)
    assert result["statistics"]["totalSteps"] == 3


def test_installed_paper_modes_cover_two_pgp_methods_and_pme(temp_dir):
    work = Path(temp_dir) / "paper_modes" / "input"
    inp = _write_minimal_charged_translation_case(work, mode="pgp_full", steps=3)

    summary = run_paper_modes(
        inp,
        Path(temp_dir) / "paper_modes" / "runs",
        modes=PAPER_MODES,
        seed=12345,
        steps=3,
        move_prob="0 0 1 0",
        use_conf_bias="no",
        use_cavity_bias="no",
        backend="pybind",
    )

    assert tuple(summary["modes"]) == PAPER_MODES
    assert set(summary["runs"]) == {"pgp_full", "pme", "pgp_full_pme"}
    assert summary["record_counts"] == {"pgp_full": 3, "pme": 3, "pgp_full_pme": 3}
    assert {"E-C", "D-C", "E-D"} <= set(summary["pairs"])
    assert summary["pairs"]["E-D"]["all"]["max_abs"] < 1.0e-8
    assert Path(summary["summary_path"]).exists()
