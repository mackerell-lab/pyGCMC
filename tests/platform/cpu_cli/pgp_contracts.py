"""
Additional PGP mode contracts for gcmc_cpu (A=direct, B=pgp_host, C=pgp_full).

These tests are end-to-end and do not depend on stdout/stderr log text.
They assert behavior via --dump-accept (JSONL) and --dump-params (JSON).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from cpu_cli.inp_units_compat_helpers import _first_accept_record, _run_gcmc_cpu, _write_inp


def _write_text(path: Path, content: str) -> None:
    path.write_text(content.strip() + "\n")


def test_pgp_host_does_not_include_nonfixed_long_range_beyond_cutoff(gcmc_cpu, temp_dir):
    """
    Mode B (pgp_host) must not add long-range Coulomb interactions from non-fixed residues.

    Construct:
    - One fixed residue (HST) with charge 0 to satisfy pgp_host's "must have fixed" requirement.
    - One non-fixed charged FRG at r > cutoff from the insertion region.
    Expectation:
    - deltaU for inserting another FRG is ~0 (no LJ, no Coulomb within cutoff, and host is neutral).
    """
    work = Path(temp_dir) / "pgp_host_no_nonfixed_long_range"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1
ATOM      1  Q0  HST A   1      10.000  10.000  10.000  1.00  0.00           C
ATOM      2  QP  FRG A   2      10.000  10.000  10.000  1.00  0.00           C
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
1   Q0    1      HST      Q0    1     0.000   12.011

[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   QP    1      FRG      QP    1     1.000   12.011

[ system ]
HostPlusOneGuest

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
Q0      0       12.011 0.000   A      0.300   0.000
QP      0       12.011 0.000   A      0.300   0.000
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

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "accept.jsonl"

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
random_seed:12345
energy_method:pgp_host
par:{par}
fragitp:{frag}
fragname:FRG
fragconc:1.0
fragmuex:0.0

pdb:{pdb}
top:{top}
box_size:50.0 50.0 50.0
gcmc_region:box 19.5 9.5 9.5 20.5 10.5 10.5
cutoff:6.0
temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:1 0 0 0

use_cavity_bias:no
use_conf_bias:no
""",
    )

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work,
        inp=inp,
        out_prefix=out_prefix,
        extra_args=["--dump-accept", str(accept_log)],
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="insertion", species="FRG")
    assert float(rec["deltaU"]) == pytest.approx(0.0, abs=1e-6)


def test_pgp_full_translation_self_exclusion_keeps_deltaU_zero_in_one_particle_system(gcmc_cpu, temp_dir):
    """
    Mode C (pgp_full) must exclude the moved residue from the background grid.

    With only one charged particle in the entire system:
    - Energy should be purely self-energy (no pair interactions).
    - Under continuous Ewald, translating the particle must not change energy (ΔU = 0).
    - In standard PME (Mode D), a finite mesh introduces a small mesh-self fluctuation; we
      therefore require Mode E (pgp_full_pme) to match PME, while Mode C stays "more physical".
    """
    work = Path(temp_dir) / "pgp_full_translation_self_exclusion"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1
ATOM      1  QP  FRG A   1      10.000  10.000  10.000  1.00  0.00           C
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
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   QP    1      FRG      QP    1     1.000   12.011

[ system ]
SingleCharged

[ molecules ]
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

        out_prefix = run_dir / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        accept_log = run_dir / "out" / "accept.jsonl"

        inp = run_dir / "test.inp"
        _write_inp(
            inp,
            f"""
random_seed:42
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
max_translation:5.0
temperature:300.0
moves_per_step:1
mcsteps:10
nprint:1
mc_move_prob:0 0 1 0

use_cavity_bias:no
use_conf_bias:no
""",
        )

        result = _run_gcmc_cpu(
            gcmc_cpu,
            workdir=run_dir,
            inp=inp,
            out_prefix=out_prefix,
            extra_args=["--dump-accept", str(accept_log)],
            timeout=60,
        )
        assert result.returncode == 0, result.stdout + result.stderr

        records = [
            json.loads(line)
            for line in accept_log.read_text().splitlines()
            if line.strip()
        ]
        deltas = [
            float(r["deltaU"])
            for r in records
            if str(r.get("move", "")).strip().lower() == "translation"
            and str(r.get("species", "")).strip().upper() == "FRG"
        ]
        assert deltas, f"No translation/FRG records found in {accept_log}"
        return deltas

    deltas_pgp = _run("pgp_full")
    deltas_pme = _run("pme")
    deltas_pgp_pme = _run("pgp_full_pme")

    assert deltas_pgp == pytest.approx([0.0] * len(deltas_pgp), abs=1e-6)
    assert deltas_pgp_pme == pytest.approx(deltas_pme, rel=1e-10, abs=1e-6)
    assert max(abs(x) for x in deltas_pme) > 1e-4


def test_pgp_host_cbmc_trial_energies_include_long_range_host_coulomb(gcmc_cpu, temp_dir):
    """
    Under CBMC insertion, Mode B must evaluate trial energies using the PGP(host-only) backend.

    Setup:
    - One fixed charged host (+1) at (10,10,10) Å.
    - Insert a charged fragment (-1) in a region at r > cutoff, so DIRECT+cutoff gives 0.
    Expectation:
    - direct: cbmcTrialEnergies ~ 0 and deltaU ~ 0
    - pgp_host: cbmcTrialEnergies contain significant negative values and deltaU < 0
    """
    work = Path(temp_dir) / "pgp_host_cbmc_trial_energies"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "host.pdb"
    _write_text(
        pdb,
        """
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1
ATOM      1  QP  HST A   1      10.000  10.000  10.000  1.00  0.00           C
END
""",
    )

    top = work / "host.top"
    _write_text(
        top,
        """
[ defaults ]
1 2 yes 1.0 1.0

[ moleculetype ]
HST  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   QP    1      HST      QP    1     1.000   12.011

[ system ]
HostOnly

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
QM      0       12.011 0.000   A      0.300   0.000
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
1   QM    1      FRG      QM    1     -1.000  12.011
""",
    )

    def _run(method: str) -> dict:
        out_prefix = work / method / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        accept_log = work / method / "out" / "accept.jsonl"
        inp = work / method / "test.inp"
        inp.parent.mkdir(parents=True, exist_ok=True)

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
box_size:50.0 50.0 50.0
gcmc_region:box 19.5 9.5 9.5 20.5 10.5 10.5
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
            workdir=work / method,
            inp=inp,
            out_prefix=out_prefix,
            extra_args=["--dump-accept", str(accept_log)],
            timeout=60,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        return _first_accept_record(accept_log, move="insertion", species="FRG")

    rec_direct = _run("direct")
    assert int(rec_direct.get("cbmcTrials", 0)) >= 5
    energies_direct = [float(x) for x in rec_direct.get("cbmcTrialEnergies", [])]
    assert len(energies_direct) == int(rec_direct["cbmcTrials"])
    selected_direct = float(rec_direct.get("cbmcSelectedEnergy", 0.0))
    assert float(rec_direct["deltaU"]) == pytest.approx(selected_direct, rel=1e-6, abs=1e-6)
    assert min(abs(selected_direct - e) for e in energies_direct) <= 1e-6
    assert max(abs(e) for e in energies_direct) < 1e-6
    assert abs(float(rec_direct["deltaU"])) < 1e-6

    rec_pgp = _run("pgp_host")
    assert int(rec_pgp.get("cbmcTrials", 0)) >= 5
    energies_pgp = [float(x) for x in rec_pgp.get("cbmcTrialEnergies", [])]
    assert len(energies_pgp) == int(rec_pgp["cbmcTrials"])
    selected_pgp = float(rec_pgp.get("cbmcSelectedEnergy", 0.0))
    assert float(rec_pgp["deltaU"]) == pytest.approx(selected_pgp, rel=1e-6, abs=1e-6)
    assert min(abs(selected_pgp - e) for e in energies_pgp) <= 1e-6
    assert min(energies_pgp) < -1.0
    assert float(rec_pgp["deltaU"]) < -1.0
