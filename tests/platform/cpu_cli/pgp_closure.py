"""
PGP closure/consistency regressions for gcmc_cpu.

These tests are end-to-end (run the compiled gcmc_cpu binary) and avoid
stdout/stderr log text matching. Assertions use --dump-accept JSONL.
"""

from __future__ import annotations

import statistics
from pathlib import Path

import pytest

from cpu_cli.inp_units_compat_helpers import _first_accept_record, _run_gcmc_cpu, _write_inp
from cpu_cli.physical_contracts import _write_charged_single_atom_system


def _write_text(path: Path, content: str) -> None:
    path.write_text(content.strip() + "\n")


def _run_pgp_full_cbmc_with_background_charge(
    gcmc_cpu: str, *, workdir: Path, background_charge: float
) -> dict:
    workdir.mkdir(parents=True, exist_ok=True)

    pdb = workdir / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1
ATOM      1  Q0  HST A   1       5.000   5.000   5.000  1.00  0.00           C
ATOM      2  QP  FRG A   2      10.000  10.000  10.000  1.00  0.00           C
END
""",
    )

    top = workdir / "sys.top"
    _write_text(
        top,
        f"""
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
1   QP    1      FRG      QP    1     {background_charge:.6f}   12.011

[ system ]
PgpFullCbmcBackgroundCharge

[ molecules ]
HST  1
FRG  1
""",
    )

    par = workdir / "par.itp"
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

    frag = workdir / "frag.itp"
    _write_text(
        frag,
        """
[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   QP    1      FRG      QP    1     1.000   12.011
""",
    )

    out_prefix = workdir / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = workdir / "out" / "accept.jsonl"

    inp = workdir / "run.inp"
    _write_inp(
        inp,
        f"""
random_seed:12345
energy_method:pgp_full
par:{par}
fragitp:{frag}
fragname:FRG
fragconc:1.0
fragmuex:0.0

pdb:{pdb}
top:{top}
box_size:50.0 50.0 50.0
gcmc_region:box 29.5 9.5 9.5 30.5 10.5 10.5
cutoff:15.0
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
        workdir=workdir,
        inp=inp,
        out_prefix=out_prefix,
        extra_args=["--dump-accept", str(accept_log)],
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="insertion", species="FRG")
    cbmc_trials = int(rec.get("cbmcTrials", 0))
    assert cbmc_trials >= 5
    energies = [float(x) for x in rec.get("cbmcTrialEnergies", [])]
    assert len(energies) == cbmc_trials

    selected = float(rec.get("cbmcSelectedEnergy", 0.0))
    delta_u = float(rec["deltaU"])
    assert delta_u == pytest.approx(selected, rel=1e-6, abs=1e-6)
    assert min(abs(selected - e) for e in energies) <= 1e-6

    return {"deltaU": delta_u, "cbmcSelectedEnergy": selected, "cbmcTrialEnergies": energies}


def test_pgp_full_cbmc_trial_energies_shift_with_nonfixed_background_charge(gcmc_cpu, temp_dir):
    """
    Mode C (pgp_full) must include non-fixed background charges in the reciprocal grid.

    We run the same CBMC insertion twice with identical RNG seed/geometry, but with the
    existing (non-fixed) FRG background charge set to 0 vs +1 in the topology.
    Self energy is identical; any energy shift must come from background electrostatics.
    """
    base = Path(temp_dir) / "pgp_full" / "cbmc_background_charge_shift"

    rec_neutral = _run_pgp_full_cbmc_with_background_charge(
        gcmc_cpu, workdir=base / "bg_q0", background_charge=0.0
    )
    rec_charged = _run_pgp_full_cbmc_with_background_charge(
        gcmc_cpu, workdir=base / "bg_qp", background_charge=1.0
    )

    mean_neutral = statistics.mean(rec_neutral["cbmcTrialEnergies"])
    mean_charged = statistics.mean(rec_charged["cbmcTrialEnergies"])

    shift_trials = mean_charged - mean_neutral
    shift_selected = rec_charged["deltaU"] - rec_neutral["deltaU"]

    # Order-of-magnitude contract: switching a background charge 0 -> +1 must produce a clear
    # repulsive shift, but we do not hard-code a vacuum 1/r estimate because periodic Ewald/PGP
    # includes image/background effects and CBMC selects among trial positions.
    assert shift_trials > 10.0
    assert shift_trials < 100.0
    assert shift_selected > 10.0
    assert shift_selected < 100.0


def test_pgp_full_deletion_deltaU_is_negative_of_insertion_for_charged_system(gcmc_cpu, temp_dir):
    """
    Mode C (pgp_full): for the same configuration, deletion must report deltaU = -insertion deltaU.

    This is a strong end-to-end closure check that catches:
    - grid background inclusion/exclusion mistakes
    - missing/extra self-energy terms
    - per-move backend mismatches between insertion and deletion
    """
    work = Path(temp_dir) / "pgp_full" / "charged_deltaU_symmetry"
    work.mkdir(parents=True, exist_ok=True)

    files = _write_charged_single_atom_system(
        work,
        mol_charge=0.400,
        frag_charge=-0.200,
        mol_sigma_nm=0.320,
        mol_eps_kj_mol=0.000,
        frag_sigma_nm=0.280,
        frag_eps_kj_mol=0.000,
    )

    ins_prefix = work / "ins" / "gcmc"
    ins_prefix.parent.mkdir(parents=True, exist_ok=True)
    ins_accept = work / "ins" / "acceptance.jsonl"
    ins_inp = work / "ins" / "test.inp"
    _write_inp(
        ins_inp,
        f"""
random_seed:8642
energy_method:pgp_full
par:{files["par"]}
fragitp:{files["frag"]}
fragname:FRG
fragconc:0.0
fragmuex:30.0

pdb:{files["pdb"]}
top:{files["top_ins"]}
box_size:30.0 30.0 30.0
gcmc_region:box 18.0 15.0 15.0 19.0 16.0 16.0
cutoff:6.0
temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:1 0 0 0
""",
    )

    ins_result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work / "ins",
        inp=ins_inp,
        out_prefix=ins_prefix,
        extra_args=["--dump-accept", str(ins_accept)],
        timeout=60,
    )
    assert ins_result.returncode == 0, ins_result.stdout + ins_result.stderr

    ins_rec = _first_accept_record(ins_accept, move="insertion", species="FRG")
    assert bool(ins_rec.get("accepted")) is True

    deltaU_ins = float(ins_rec["deltaU"])
    assert abs(deltaU_ins) > 5.0

    final_pdb = Path(f"{ins_prefix}_final.pdb")
    assert final_pdb.exists()

    del_prefix = work / "del" / "gcmc"
    del_prefix.parent.mkdir(parents=True, exist_ok=True)
    del_accept = work / "del" / "acceptance.jsonl"
    del_inp = work / "del" / "test.inp"
    _write_inp(
        del_inp,
        f"""
random_seed:97531
energy_method:pgp_full
par:{files["par"]}
fragitp:{files["frag"]}
fragname:FRG
fragconc:0.0
fragmuex:30.0

pdb:{final_pdb}
top:{files["top_del"]}
box_size:30.0 30.0 30.0
gcmc_region:box 18.0 15.0 15.0 19.0 16.0 16.0
cutoff:6.0
temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:0 1 0 0
""",
    )

    del_result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work / "del",
        inp=del_inp,
        out_prefix=del_prefix,
        extra_args=["--dump-accept", str(del_accept)],
        timeout=60,
    )
    assert del_result.returncode == 0, del_result.stdout + del_result.stderr

    del_rec = _first_accept_record(del_accept, move="deletion", species="FRG")
    assert int(del_rec.get("nBefore", 0)) == 1

    deltaU_del = float(del_rec["deltaU"])
    assert deltaU_del == pytest.approx(-deltaU_ins, rel=1e-4, abs=1e-2)
