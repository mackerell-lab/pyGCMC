"""
Mode E (pgp_full_pme) and Mode D (pme) must be energy-consistent.

For gcmc_cpu we require that `energy_method:pme` reproduces the same per-move ΔU
as `energy_method:pgp_full_pme` when running the same random seed and move proposals.
This prevents backend-dependent sampling bias in μVT simulations.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from cpu_cli.inp_units_compat_helpers import _run_gcmc_cpu, _write_inp


def _write_text(path: Path, content: str) -> None:
    path.write_text(content.strip() + "\n")


def _load_translation_delta_us(path: Path) -> list[float]:
    deltas: list[float] = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        if str(rec.get("move", "")).strip().lower() != "translation":
            continue
        if str(rec.get("species", "")).strip().upper() != "FRG":
            continue
        deltas.append(float(rec["deltaU"]))
    return deltas


def _load_move_delta_us(path: Path, *, move: str, species: str) -> tuple[list[float], list[bool]]:
    deltas: list[float] = []
    accepted: list[bool] = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        if str(rec.get("move", "")).strip().lower() != move.strip().lower():
            continue
        if str(rec.get("species", "")).strip().upper() != species.strip().upper():
            continue
        deltas.append(float(rec["deltaU"]))
        accepted.append(bool(rec["accepted"]))
    return deltas, accepted


def test_pgp_full_pme_and_pme_translation_deltaU_match(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "pgp_pme_equivalence" / "translation_deltaU_match"
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

    def _run(method: str) -> list[float]:
        run_dir = work / method
        run_dir.mkdir(parents=True, exist_ok=True)

        inp = run_dir / "run.inp"
        _write_inp(
            inp,
            f"""
random_seed:12345
energy_method:{method}
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
mcsteps:20
nprint:20
mc_move_prob:0 0 1 0
max_translation:1.0
max_rotation:0.0

use_cavity_bias:no
use_conf_bias:no
""",
        )

        out_prefix = run_dir / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        accept_log = run_dir / "out" / "accept.jsonl"

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

        deltas = _load_translation_delta_us(accept_log)
        assert len(deltas) == 20
        return deltas

    deltas_pgp = _run("pgp_full_pme")
    deltas_pme = _run("pme")

    assert deltas_pme == pytest.approx(deltas_pgp, rel=1e-10, abs=1e-6)


def test_pgp_full_pme_and_pme_insertion_deltaU_match_neutral_lj(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "pgp_pme_equivalence" / "insertion_deltaU_match_rejected"
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
1   Q     1      HST      Q     1     0.000   12.011

[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   Q     1      FRG      Q     1     0.000   12.011

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
Q       0       12.011 0.000   A      0.300   0.200
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
1   Q     1      FRG      Q     1     0.000  12.011
""",
    )

    def _run(method: str) -> list[float]:
        run_dir = work / method
        run_dir.mkdir(parents=True, exist_ok=True)

        inp = run_dir / "run.inp"
        _write_inp(
            inp,
            f"""
random_seed:12345
energy_method:{method}
par:{par}
fragitp:{frag}
fragname:FRG
fragconc:1.0
fragmuex:0.0

pdb:{pdb}
top:{top}
box_size:40.0 40.0 40.0
cutoff:6.0
temperature:1e30
moves_per_step:1
mcsteps:20
nprint:20
mc_move_prob:1 0 0 0
max_translation:0.0
max_rotation:0.0

use_cavity_bias:no
use_conf_bias:no
""",
        )

        out_prefix = run_dir / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        accept_log = run_dir / "out" / "accept.jsonl"

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

        deltas, accepted = _load_move_delta_us(accept_log, move="insertion", species="FRG")
        assert len(deltas) == 20
        assert accepted == [True] * 20
        return deltas

    deltas_pgp = _run("pgp_full_pme")
    deltas_pme = _run("pme")

    assert deltas_pme == pytest.approx(deltas_pgp, rel=1e-10, abs=1e-6)


def test_pgp_full_pme_and_pme_deletion_deltaU_match_neutral_lj(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "pgp_pme_equivalence" / "deletion_deltaU_match_rejected"
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
1   Q     1      HST      Q     1     0.000   12.011

[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   Q     1      FRG      Q     1     0.000   12.011

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
Q       0       12.011 0.000   A      0.300   0.200
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
1   Q     1      FRG      Q     1     0.000  12.011
""",
    )

    def _run(method: str) -> list[float]:
        run_dir = work / method
        run_dir.mkdir(parents=True, exist_ok=True)

        inp = run_dir / "run.inp"
        _write_inp(
            inp,
            f"""
random_seed:12345
energy_method:{method}
par:{par}
fragitp:{frag}
fragname:FRG
fragconc:1e-4
fragmuex:0.0

pdb:{pdb}
top:{top}
box_size:40.0 40.0 40.0
cutoff:6.0
temperature:1e30
moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:0 1 0 0
max_translation:0.0
max_rotation:0.0

use_cavity_bias:no
use_conf_bias:no
""",
        )

        out_prefix = run_dir / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        accept_log = run_dir / "out" / "accept.jsonl"

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

        deltas, accepted = _load_move_delta_us(accept_log, move="deletion", species="FRG")
        assert len(deltas) == 1
        assert accepted == [True]
        return deltas

    deltas_pgp = _run("pgp_full_pme")
    deltas_pme = _run("pme")

    assert deltas_pme == pytest.approx(deltas_pgp, rel=1e-10, abs=1e-6)


def test_pgp_full_pme_and_pme_insertion_deltaU_match_charged_coulomb(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "pgp_pme_equivalence" / "insertion_deltaU_match_charged"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1
ATOM      1  QN  HST A   1      10.000  20.000  20.000  1.00  0.00           C
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
1   QN    1      HST      QN    1    -1.000   12.011

[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   QP    1      FRG      QP    1     1.000   12.011

[ system ]
ChargedInsertion

[ molecules ]
HST  1
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
QP      0       12.011 0.000   A      0.300   0.000
QN      0       12.011 0.000   A      0.300   0.000
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
1   QP    1      FRG      QP    1     1.000  12.011
""",
    )

    def _run(method: str) -> tuple[list[float], list[bool]]:
        run_dir = work / method
        run_dir.mkdir(parents=True, exist_ok=True)

        inp = run_dir / "run.inp"
        _write_inp(
            inp,
            f"""
random_seed:12345
energy_method:{method}
par:{par}
fragitp:{frag}
fragname:FRG
fragconc:1.0
fragmuex:0.0

pdb:{pdb}
top:{top}
box_size:40.0 40.0 40.0
cutoff:6.0
temperature:1e30
moves_per_step:1
mcsteps:5
nprint:5
mc_move_prob:1 0 0 0
max_translation:0.0
max_rotation:0.0

use_cavity_bias:no
use_conf_bias:no
""",
        )

        out_prefix = run_dir / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        accept_log = run_dir / "out" / "accept.jsonl"

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

        deltas, accepted = _load_move_delta_us(accept_log, move="insertion", species="FRG")
        assert len(deltas) == 5
        return deltas, accepted

    deltas_pgp, accepted_pgp = _run("pgp_full_pme")
    deltas_pme, accepted_pme = _run("pme")

    assert accepted_pme == accepted_pgp
    assert deltas_pme == pytest.approx(deltas_pgp, rel=1e-10, abs=1e-6)


def test_pgp_full_pme_and_pme_deletion_deltaU_match_charged_coulomb(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "pgp_pme_equivalence" / "deletion_deltaU_match_charged"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1
ATOM      1  QN  HST A   1      10.000  20.000  20.000  1.00  0.00           C
ATOM      2  QP  FRG B   1      35.000  20.000  20.000  1.00  0.00           C
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
1   QN    1      HST      QN    1    -1.000   12.011

[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   QP    1      FRG      QP    1     1.000   12.011

[ system ]
ChargedDeletion

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
QP      0       12.011 0.000   A      0.300   0.000
QN      0       12.011 0.000   A      0.300   0.000
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
1   QP    1      FRG      QP    1     1.000  12.011
""",
    )

    def _run(method: str) -> list[float]:
        run_dir = work / method
        run_dir.mkdir(parents=True, exist_ok=True)

        inp = run_dir / "run.inp"
        _write_inp(
            inp,
            f"""
random_seed:12345
energy_method:{method}
par:{par}
fragitp:{frag}
fragname:FRG
fragconc:1e-4
fragmuex:0.0

pdb:{pdb}
top:{top}
box_size:40.0 40.0 40.0
cutoff:6.0
temperature:1e30
moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:0 1 0 0
max_translation:0.0
max_rotation:0.0

use_cavity_bias:no
use_conf_bias:no
""",
        )

        out_prefix = run_dir / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        accept_log = run_dir / "out" / "accept.jsonl"

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

        deltas, accepted = _load_move_delta_us(accept_log, move="deletion", species="FRG")
        assert len(deltas) == 1
        assert accepted == [True]
        return deltas

    deltas_pgp = _run("pgp_full_pme")
    deltas_pme = _run("pme")

    assert deltas_pme == pytest.approx(deltas_pgp, rel=1e-10, abs=1e-6)


def test_pgp_full_pme_and_pme_rotation_deltaU_match_neutral_lj(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "pgp_pme_equivalence" / "rotation_deltaU_match_neutral_lj"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1
ATOM      1  QH  HST A   1      12.000  20.000  20.000  1.00  0.00           C
ATOM      2  QP  FRG B   1      32.000  20.000  20.000  1.00  0.00           C
ATOM      3  QN  FRG B   1      33.000  20.000  20.000  1.00  0.00           C
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
1   QH    1      HST      QH    1     0.000   12.011

[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   QP    1      FRG      QP    1     0.000   12.011
2   QN    1      FRG      QN    2     0.000   12.011

[ system ]
NeutralLJRotation

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
QH      0       12.011 0.000   A      0.300   0.200
QP      0       12.011 0.000   A      0.300   0.200
QN      0       12.011 0.000   A      0.300   0.200
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
1   QP    1      FRG      QP    1     0.000  12.011
2   QN    1      FRG      QN    2     0.000  12.011
""",
    )

    def _run(method: str) -> tuple[list[float], list[bool]]:
        run_dir = work / method
        run_dir.mkdir(parents=True, exist_ok=True)

        inp = run_dir / "run.inp"
        _write_inp(
            inp,
            f"""
random_seed:12345
energy_method:{method}
par:{par}
fragitp:{frag}
fragname:FRG
fragconc:0.0
fragmuex:-1000.0

pdb:{pdb}
top:{top}
box_size:40.0 40.0 40.0
cutoff:6.0
temperature:300.0
moves_per_step:1
mcsteps:20
nprint:20
mc_move_prob:0 0 0 1
max_translation:0.0
max_rotation:0.3

use_cavity_bias:no
use_conf_bias:no
""",
        )

        out_prefix = run_dir / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        accept_log = run_dir / "out" / "accept.jsonl"

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

        deltas, accepted = _load_move_delta_us(accept_log, move="rotation", species="FRG")
        assert len(deltas) == 20
        return deltas, accepted

    deltas_pgp, accepted_pgp = _run("pgp_full_pme")
    deltas_pme, accepted_pme = _run("pme")

    assert accepted_pme == accepted_pgp
    assert deltas_pme == pytest.approx(deltas_pgp, rel=1e-10, abs=1e-6)
