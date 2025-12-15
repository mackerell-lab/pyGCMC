import csv
import os
import subprocess
from pathlib import Path

import pytest


GCMC_CPU_PATH = Path(__file__).parent.parent.parent / "build" / "bin" / "gcmc_cpu"


def test_itp_nonbonded_builds_lj_matrix(tmp_path: Path):
    ff_file = tmp_path / "ffnonbonded.itp"
    ff_file.write_text(
        """[ atomtypes ]
; name  bond_type  mass  charge  ptype  sigma     epsilon
A       0          1.0   0.0     A      0.200     0.500
B       0          1.0   0.0     A      0.400     2.000

[ nonbond_params ]
; i  j  func  sigma  epsilon
A  B  1  0.123  9.870
"""
    )

    inp_file = tmp_path / "test.inp"
    out_prefix = tmp_path / "out"
    lj_csv = Path(str(out_prefix) + "_lj.csv")
    frag_itp = tmp_path / "WAT.itp"
    frag_itp.write_text(
        """[ moleculetype ]
WAT 3

[ atoms ]
1 A 1 WAT O 1 0.0 1.0
2 B 1 WAT H1 2 0.0 1.0
3 B 1 WAT H2 3 0.0 1.0
"""
    )

    inp_file.write_text(
        f"""par:{ff_file}
fragitp:{frag_itp}

fragname:WAT
fragconc:55.0
fragmuex:0.0

box_size:5.0 5.0 5.0
cutoff:1.0
temperature:300
mcsteps:1
nprint:1000
eqsteps:0

op_top:{tmp_path}/out.top
op_pdb:{tmp_path}/out.pdb
"""
    )

    env = os.environ.copy()
    env["GCMC_DUMP_LJ"] = "1"

    result = subprocess.run(
        [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--seed", "42", "--prefix", str(out_prefix)],
        cwd=str(tmp_path),
        capture_output=True,
        text=True,
        timeout=20,
        env=env,
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    assert lj_csv.exists(), "Expected LJ CSV export file was not created"

    with lj_csv.open(newline="") as f:
        rows = list(csv.DictReader(f))

    # For two atom types, we expect a 2x2 matrix = 4 entries.
    assert len(rows) == 4

    matrix = {}
    for r in rows:
        key = (int(r["i"]), int(r["j"]))
        matrix[key] = (float(r["sigma_ij"]), float(r["eps_ij"]))

    # Type indices are deterministic: std::map order -> A=0, B=1
    assert matrix[(0, 0)][0] == pytest.approx(0.200)
    assert matrix[(0, 0)][1] == pytest.approx(0.500)
    assert matrix[(1, 1)][0] == pytest.approx(0.400)
    assert matrix[(1, 1)][1] == pytest.approx(2.000)

    # nonbond_params overrides A-B mixing rule
    assert matrix[(0, 1)][0] == pytest.approx(0.123)
    assert matrix[(0, 1)][1] == pytest.approx(9.870)
    assert matrix[(1, 0)][0] == pytest.approx(0.123)
    assert matrix[(1, 0)][1] == pytest.approx(9.870)


def test_itp_pairtypes_overrides_mixing(tmp_path: Path):
    ff_file = tmp_path / "ffnonbonded.itp"
    ff_file.write_text(
        """[ atomtypes ]
; name  bond_type  mass  charge  ptype  sigma     epsilon
A       0          1.0   0.0     A      0.200     0.500
B       0          1.0   0.0     A      0.400     2.000

[ pairtypes ]
; i  j  func  sigma  epsilon
A  B  1  0.321  8.760
"""
    )

    inp_file = tmp_path / "test.inp"
    out_prefix = tmp_path / "out"
    lj_csv = Path(str(out_prefix) + "_lj.csv")
    frag_itp = tmp_path / "WAT.itp"
    frag_itp.write_text(
        """[ moleculetype ]
WAT 3

[ atoms ]
1 A 1 WAT O 1 0.0 1.0
2 B 1 WAT H1 2 0.0 1.0
3 B 1 WAT H2 3 0.0 1.0
"""
    )

    inp_file.write_text(
        f"""par:{ff_file}
fragitp:{frag_itp}

fragname:WAT
fragconc:55.0
fragmuex:0.0

box_size:5.0 5.0 5.0
cutoff:1.0
temperature:300
mcsteps:1
nprint:1000
eqsteps:0

op_top:{tmp_path}/out.top
op_pdb:{tmp_path}/out.pdb
"""
    )

    env = os.environ.copy()
    env["GCMC_DUMP_LJ"] = "1"

    result = subprocess.run(
        [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--seed", "42", "--prefix", str(out_prefix)],
        cwd=str(tmp_path),
        capture_output=True,
        text=True,
        timeout=20,
        env=env,
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    assert lj_csv.exists(), "Expected LJ CSV export file was not created"

    with lj_csv.open(newline="") as f:
        rows = list(csv.DictReader(f))

    matrix = {(int(r["i"]), int(r["j"])): (float(r["sigma_ij"]), float(r["eps_ij"])) for r in rows}

    # Type indices are deterministic: std::map order -> A=0, B=1
    assert matrix[(0, 1)][0] == pytest.approx(0.321)
    assert matrix[(0, 1)][1] == pytest.approx(8.760)
    assert matrix[(1, 0)][0] == pytest.approx(0.321)
    assert matrix[(1, 0)][1] == pytest.approx(8.760)
