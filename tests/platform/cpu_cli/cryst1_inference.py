"""
CRYST1-less PDB handling in gcmc_gpu-style decks.

Some legacy/example inputs omit the PDB CRYST1 record. In gcmc_gpu/opencl-style INP,
`box_size` commonly denotes the *GCMC region* size, while `sys_center` can be used to
infer the periodic box when CRYST1 is missing.

This test is file-driven and validates behavior via the output PDB CRYST1 record.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from .inp_units_compat_helpers import _read_cryst1_box_angstrom, _run_gcmc_cpu, _write_inp


def test_missing_cryst1_uses_sys_center_to_infer_periodic_box_in_gcmc_gpu_units(
    gcmc_cpu, temp_dir
):
    work = Path(temp_dir) / "cryst1_missing_infer_from_sys_center"
    work.mkdir(parents=True, exist_ok=True)

    # PDB intentionally omits CRYST1.
    pdb = work / "one_water_no_cryst1.pdb"
    pdb.write_text(
        "\n".join(
            [
                "ATOM      1  O   SOL A   1       0.000   0.000   0.000  1.00  0.00           O",
                "ATOM      2  H1  SOL A   1       0.957   0.000   0.000  1.00  0.00           H",
                "ATOM      3  H2  SOL A   1      -0.239   0.927   0.000  1.00  0.00           H",
                "END",
            ]
        )
        + "\n"
    )

    # Minimal GROMACS topology for the SOL residue.
    top = work / "one_water.top"
    top.write_text(
        "\n".join(
            [
                "[ defaults ]",
                "1 2",
                "",
                "[ atomtypes ]",
                "O   8  15.9994  0.0  A  3.15061e-01  6.36386e-01",
                "H   1   1.008   0.0  A  0.00000e+00  0.00000e+00",
                "",
                "[ moleculetype ]",
                "SOL  2",
                "",
                "[ atoms ]",
                "; nr  type  resnr  residue  atom  cgnr  charge  mass",
                "1   O    1   SOL  O    1   -0.834   15.9994",
                "2   H    1   SOL  H1   1    0.417    1.008",
                "3   H    1   SOL  H2   1    0.417    1.008",
                "",
                "[ system ]",
                "OneWater",
                "",
                "[ molecules ]",
                "SOL 1",
                "",
            ]
        )
        + "\n"
    )

    # In gcmc_gpu-style INP units (default inp_units:auto), coordinates and lengths are in Å.
    # `box_size` denotes the *region* size; with missing CRYST1, use `sys_center` to infer the periodic box.
    box_size_a = 10.0
    sys_center_a = 20.0
    expected_box_a = 2.0 * sys_center_a

    inp = work / "run.inp"
    _write_inp(
        inp,
        f"""
pdb:{pdb}
top:{top}
op_pdb:out.pdb
op_top:out.top

box_size:{box_size_a} {box_size_a} {box_size_a}
sys_center:{sys_center_a} {sys_center_a} {sys_center_a}
gc_center:{sys_center_a} {sys_center_a} {sys_center_a}

temperature:300.0
cutoff:8.0

fragname:SOL
fragconc:55.0
fragmuex:0.0

mcsteps:0
nprint:1
""",
    )

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work,
        inp=inp,
        out_prefix=out_prefix,
        extra_args=["--seed", "123"],
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()
    cryst = _read_cryst1_box_angstrom(final_pdb)
    assert cryst == pytest.approx((expected_box_a, expected_box_a, expected_box_a), abs=1e-3)

