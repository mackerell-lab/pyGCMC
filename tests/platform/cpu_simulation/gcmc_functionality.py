"""
Functional validation tests for gcmc_cpu
测试 gcmc_cpu 的基本功能是否正确运行

Tests:
1. mc_move_prob parameter parsing and functionality
2. Statistics DAT file output
3. Basic GCMC simulation runs
4. INP parameter handling
5. Numerical correctness (energy, MC statistics, weights)
"""

import json
import pytest
import numpy as np
import subprocess
from pathlib import Path

# Try to import pygcmc bindings for energy verification
try:
    import pygcmc
    import pygcmc.platform.cpu as cpu_platform
    HAS_PYGCMC = True
except ImportError:
    HAS_PYGCMC = False

# Path to gcmc_cpu executable
GCMC_CPU_PATH = Path(__file__).parent.parent.parent.parent / "build" / "bin" / "gcmc_cpu"

MOVE_TYPES = ("insertion", "deletion", "translation", "rotation")


def _load_accept_records(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def _load_params(path: Path) -> dict:
    return json.loads(path.read_text())


def _parse_itp_atomtypes(path: Path) -> dict[str, tuple[float, float]]:
    section = ""
    atomtypes: dict[str, tuple[float, float]] = {}
    for raw in path.read_text().splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and "]" in line:
            section = line[1:line.index("]")].strip().lower()
            continue
        if section == "atomtypes":
            tokens = line.split()
            if len(tokens) >= 7:
                atomtypes[tokens[0]] = (float(tokens[-2]), float(tokens[-1]))
    if not atomtypes:
        raise AssertionError(f"No atomtypes parsed from {path}")
    return atomtypes


def _expected_cdf(weights: list[float]) -> list[float]:
    if all(w <= 0 for w in weights):
        weights = [1.0, 1.0, 1.0, 1.0]
    total = sum(weights)
    if total <= 0:
        total = 1.0
    c0 = weights[0] / total
    c1 = c0 + weights[1] / total
    c2 = c1 + weights[2] / total
    return [c0, c1, c2, 1.0]


def _fragment_index(params: dict, name: str) -> int:
    names = params.get("fragment", {}).get("names", [])
    target = name.strip().upper()
    for idx, frag_name in enumerate(names):
        if str(frag_name).strip().upper() == target:
            return idx
    if not names:
        return 0
    raise AssertionError(f"Fragment name {name} not found in params: {names}")


def _fragment_cdf(params: dict, name: str) -> list[float]:
    idx = _fragment_index(params, name)
    cdf_list = params.get("fragment", {}).get("move_cdf", [])
    if idx >= len(cdf_list):
        raise AssertionError(f"move_cdf missing index {idx} for fragment {name}")
    return [float(x) for x in cdf_list[idx]]


def _move_counts(records: list[dict], *, species: str | None = None) -> dict[str, int]:
    counts = {move: 0 for move in MOVE_TYPES}
    for rec in records:
        if species is not None:
            if str(rec.get("species", "")).strip().upper() != species.strip().upper():
                continue
        move = str(rec.get("move", "")).strip().lower()
        if move in counts:
            counts[move] += 1
    return counts


def _move_fractions(counts: dict[str, int]) -> dict[str, float]:
    total = sum(counts.values())
    if total <= 0:
        raise AssertionError("No move attempts found in acceptance log")
    return {move: counts[move] / total for move in MOVE_TYPES}


def _assert_move_fractions_close(actual: dict[str, float], expected: dict[str, float], *, tol: float) -> None:
    for move in MOVE_TYPES:
        assert abs(actual[move] - expected[move]) <= tol, (
            f"{move} fraction {actual[move]:.3f} differs from expected {expected[move]:.3f}"
        )


def _run_with_accept_log(inp_file: Path, tmp_path: Path, *, seed: str = "42") -> tuple[subprocess.CompletedProcess, Path, Path]:
    accept_log = tmp_path / "accept.jsonl"
    params_json = tmp_path / "params.json"
    out_prefix = tmp_path / "out"
    result = subprocess.run(
        [
            str(GCMC_CPU_PATH),
            "--inp",
            str(inp_file),
            "--seed",
            seed,
            "--prefix",
            str(out_prefix),
            "--dump-accept",
            str(accept_log),
            "--dump-params",
            str(params_json),
        ],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=str(tmp_path),
    )
    return result, accept_log, out_prefix


class TestMcMoveProb:
    """Test mc_move_prob parameter functionality"""

    @staticmethod
    def create_minimal_system(tmpdir, mc_move_prob_values=None, mcsteps=100):
        """Create minimal water system for testing mc_move_prob"""

        # Minimal PDB
        pdb_file = tmpdir / "test.pdb"
        pdb_content = """CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1
ATOM      1  O   WAT     1      10.000  10.000  10.000  1.00  0.00
ATOM      2  H1  WAT     1      10.757  10.586  10.000  1.00  0.00
ATOM      3  H2  WAT     1       9.243  10.586  10.000  1.00  0.00
END
"""
        pdb_file.write_text(pdb_content)

        # Minimal TOP
        top_file = tmpdir / "test.top"
        top_content = """[ defaults ]
1 2 yes 0.5 0.8333

[ atomtypes ]
O   8   15.9994  -0.834  A  3.15061e-01  6.36386e-01
H   1   1.008     0.417  A  0.00000e+00  0.00000e+00

[ moleculetype ]
WAT   3

[ atoms ]
1  O   1  WAT  O   1  -0.834  15.9994
2  H   1  WAT  H1  2   0.417   1.008
3  H   1  WAT  H2  3   0.417   1.008

[ bonds ]
1  2  1  0.09572  502416.0
1  3  1  0.09572  502416.0

[ angles ]
2  1  3  1  104.52  628.02

[ system ]
Test Water

[ molecules ]
WAT  1
"""
        top_file.write_text(top_content)

        # Atomtypes file
        atp_file = tmpdir / "atomtypes.atp"
        atp_content = """O   15.9994
H    1.008
"""
        atp_file.write_text(atp_content)

        # Force field parameter file
        ff_file = tmpdir / "ffnonbonded.itp"
        ff_content = """[ atomtypes ]
O   8   15.9994  -0.834  A  3.15061e-01  6.36386e-01
H   1   1.008     0.417  A  0.00000e+00  0.00000e+00
"""
        ff_file.write_text(ff_content)

        # INP file
        inp_file = tmpdir / "test.inp"

        # Build mc_move_prob line if provided
        mc_move_line = ""
        if mc_move_prob_values:
            mc_move_line = f"mc_move_prob: {' '.join(map(str, mc_move_prob_values))}"

        inp_content = f"""par:{ff_file}
atomtypes:{atp_file}
top:{top_file}
pdb:{pdb_file}
protitp:{top_file}

fragname: water
fragconc: 55.0
fragmuex: -5.60

box_size: 20.0 20.0 20.0
cutoff: 10.0
temperature: 300
mcsteps: {mcsteps}
nprint: 20
eqsteps: 10

{mc_move_line}

op_top: {tmpdir}/output.top
op_pdb: {tmpdir}/output.pdb
"""
        inp_file.write_text(inp_content)

        return inp_file

    def test_mc_move_prob_equal_probabilities(self, tmp_path):
        """Test mc_move_prob with equal probabilities (0.25 each)"""
        inp_file = self.create_minimal_system(tmp_path, [0.25, 0.25, 0.25, 0.25], mcsteps=2000)
        result, accept_log, _ = _run_with_accept_log(inp_file, tmp_path)

        assert result.returncode == 0, f"gcmc_cpu failed: {result.stderr}"
        assert accept_log.exists(), "Acceptance log missing"
        params_json = tmp_path / "params.json"
        assert params_json.exists(), "dump-params output missing"

        params = _load_params(params_json)
        actual_cdf = _fragment_cdf(params, "water")
        expected_cdf = _expected_cdf([0.25, 0.25, 0.25, 0.25])
        assert actual_cdf == pytest.approx(expected_cdf, abs=1e-6)

    def test_mc_move_prob_biased_insertion(self, tmp_path):
        """Test mc_move_prob with bias toward insertion (0.7, 0.1, 0.1, 0.1)"""
        inp_file = self.create_minimal_system(tmp_path, [0.7, 0.1, 0.1, 0.1], mcsteps=2000)
        result, accept_log, _ = _run_with_accept_log(inp_file, tmp_path)

        assert result.returncode == 0, f"gcmc_cpu failed: {result.stderr}"
        assert accept_log.exists(), "Acceptance log missing"
        params_json = tmp_path / "params.json"
        assert params_json.exists(), "dump-params output missing"

        params = _load_params(params_json)
        actual_cdf = _fragment_cdf(params, "water")
        expected_cdf = _expected_cdf([0.7, 0.1, 0.1, 0.1])
        assert actual_cdf == pytest.approx(expected_cdf, abs=1e-6)

    def test_mc_move_prob_before_fragname(self, tmp_path):
        """Test that mc_move_prob works even when placed before fragname (time-order independent)"""

        # Create files with both water and methanol topologies
        pdb_file, top_file, atp_file, ff_file = self._create_multi_fragment_files(tmp_path)

        # INP with mc_move_prob BEFORE fragname
        inp_file = tmp_path / "test.inp"
        inp_content = f"""par:{ff_file}
atomtypes:{atp_file}
top:{top_file}
pdb:{pdb_file}
protitp:{top_file}

mc_move_prob: 0.4 0.3 0.2 0.1

fragname: water methanol
fragconc: 55.0 1.0
fragmuex: -5.60 -3.0

box_size: 20.0 20.0 20.0
cutoff: 10.0
temperature: 300
mcsteps: 2000
nprint: 200
eqsteps: 5

op_top: {tmp_path}/output.top
op_pdb: {tmp_path}/output.pdb
"""
        inp_file.write_text(inp_content)

        result, accept_log, _ = _run_with_accept_log(inp_file, tmp_path)

        assert result.returncode == 0, f"gcmc_cpu failed: {result.stderr}"

        assert accept_log.exists(), "Acceptance log missing"
        params_json = tmp_path / "params.json"
        assert params_json.exists(), "dump-params output missing"
        records = _load_accept_records(accept_log)

        seen_species = {str(r.get("species", "")).strip().upper() for r in records}
        assert {"WATER", "METHANOL"}.issubset(seen_species)

        params = _load_params(params_json)
        expected_cdf = _expected_cdf([0.4, 0.3, 0.2, 0.1])
        assert _fragment_cdf(params, "water") == pytest.approx(expected_cdf, abs=1e-6)
        assert _fragment_cdf(params, "methanol") == pytest.approx(expected_cdf, abs=1e-6)

    def _create_basic_files(self, tmpdir):
        """Helper to create basic simulation files"""
        pdb_file = tmpdir / "test.pdb"
        pdb_file.write_text("""CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1
ATOM      1  O   WAT     1      10.000  10.000  10.000  1.00  0.00
ATOM      2  H1  WAT     1      10.757  10.586  10.000  1.00  0.00
ATOM      3  H2  WAT     1       9.243  10.586  10.000  1.00  0.00
END
""")

        top_file = tmpdir / "test.top"
        top_file.write_text("""[ defaults ]
1 2 yes 0.5 0.8333
[ atomtypes ]
O   8   15.9994  -0.834  A  3.15061e-01  6.36386e-01
H   1   1.008     0.417  A  0.00000e+00  0.00000e+00
[ moleculetype ]
WAT   3
[ atoms ]
1  O   1  WAT  O   1  -0.834  15.9994
2  H   1  WAT  H1  2   0.417   1.008
3  H   1  WAT  H2  3   0.417   1.008
[ system ]
Test
[ molecules ]
WAT  1
""")

        atp_file = tmpdir / "atomtypes.atp"
        atp_file.write_text("O   15.9994\nH    1.008\n")

        ff_file = tmpdir / "ffnonbonded.itp"
        ff_file.write_text("""[ atomtypes ]
O   8   15.9994  -0.834  A  3.15061e-01  6.36386e-01
H   1   1.008     0.417  A  0.00000e+00  0.00000e+00
""")

        return pdb_file, top_file, atp_file, ff_file

    def _create_multi_fragment_files(self, tmpdir):
        """Helper to create simulation files with water and methanol topologies"""
        pdb_file = tmpdir / "test.pdb"
        pdb_file.write_text("""CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1
ATOM      1  O   WAT     1      10.000  10.000  10.000  1.00  0.00
ATOM      2  H1  WAT     1      10.757  10.586  10.000  1.00  0.00
ATOM      3  H2  WAT     1       9.243  10.586  10.000  1.00  0.00
END
""")

        top_file = tmpdir / "test.top"
        top_file.write_text("""[ defaults ]
1 2 yes 0.5 0.8333
[ atomtypes ]
O   8   15.9994  -0.834  A  3.15061e-01  6.36386e-01
H   1   1.008     0.417  A  0.00000e+00  0.00000e+00
C   6   12.011    0.145  A  3.50000e-01  2.76144e-01

[ moleculetype ]
WAT   3
[ atoms ]
1  O   1  WAT  O   1  -0.834  15.9994
2  H   1  WAT  H1  2   0.417   1.008
3  H   1  WAT  H2  3   0.417   1.008
[ bonds ]
1  2  1  0.09572  502416.0
1  3  1  0.09572  502416.0
[ angles ]
2  1  3  1  104.52  628.02

[ moleculetype ]
MEOH  6
[ atoms ]
1  C   1  MEOH  C1  1   0.145  12.011
2  O   1  MEOH  O2  2  -0.683  15.9994
3  H   1  MEOH  H3  3   0.040   1.008
4  H   1  MEOH  H4  4   0.040   1.008
5  H   1  MEOH  H5  5   0.040   1.008
6  H   1  MEOH  H6  6   0.418   1.008
[ bonds ]
1  2  1  0.143  267776.0
1  3  1  0.109  284512.0
1  4  1  0.109  284512.0
1  5  1  0.109  284512.0
2  6  1  0.0945 462750.4
[ angles ]
2  1  3  1  109.5  292.88
2  1  4  1  109.5  292.88
2  1  5  1  109.5  292.88
3  1  4  1  109.5  276.14
3  1  5  1  109.5  276.14
4  1  5  1  109.5  276.14
1  2  6  1  108.5  460.24

[ system ]
Test Multi-fragment
[ molecules ]
WAT  1
""")

        atp_file = tmpdir / "atomtypes.atp"
        atp_file.write_text("O   15.9994\nH    1.008\nC   12.011\n")

        ff_file = tmpdir / "ffnonbonded.itp"
        ff_file.write_text("""[ atomtypes ]
O   8   15.9994  -0.834  A  3.15061e-01  6.36386e-01
H   1   1.008     0.417  A  0.00000e+00  0.00000e+00
C   6   12.011    0.145  A  3.50000e-01  2.76144e-01
""")

        return pdb_file, top_file, atp_file, ff_file


class TestStatisticsOutput:
    """Test statistics DAT file output"""

    def test_statistics_dat_created(self, tmp_path):
        """Test that _statistics.dat file is created"""

        # Use helper from TestMcMoveProb
        test_helper = TestMcMoveProb()
        inp_file = test_helper.create_minimal_system(tmp_path)

        # Run gcmc_cpu in tmp_path directory
        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file)],
            capture_output=True,
            text=True,
            timeout=30,
            cwd=str(tmp_path)  # Run in tmp_path to ensure files are created there
        )

        assert result.returncode == 0, f"gcmc_cpu failed: {result.stderr}"

        # Check that statistics DAT file was created
        # Default prefix is "gcmc", so file should be gcmc_statistics.dat
        dat_files = list(tmp_path.glob("*statistics.dat"))
        assert len(dat_files) > 0, f"No statistics DAT file created. Files: {list(tmp_path.iterdir())}"

        dat_file = dat_files[0]
        dat_content = dat_file.read_text()

        # Check header
        assert "# GCMC Statistics File" in dat_content
        assert "# Generated by PyGCMC" in dat_content
        assert "Step" in dat_content
        assert "Energy" in dat_content
        assert "N_total" in dat_content
        assert "Accept_rate" in dat_content

        # Check that data lines exist
        lines = dat_content.strip().split('\n')
        data_lines = [l for l in lines if not l.startswith('#')]
        assert len(data_lines) > 0, "No data lines in DAT file"

    def test_statistics_dat_format(self, tmp_path):
        """Test that statistics DAT file has correct format"""

        test_helper = TestMcMoveProb()
        inp_file = test_helper.create_minimal_system(tmp_path)

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file)],
            capture_output=True,
            text=True,
            timeout=30,
            cwd=str(tmp_path)
        )

        assert result.returncode == 0

        dat_files = list(tmp_path.glob("*statistics.dat"))
        assert len(dat_files) > 0, f"No statistics DAT file. Files: {list(tmp_path.iterdir())}"
        dat_file = dat_files[0]
        dat_content = dat_file.read_text()

        # Parse first data line
        lines = dat_content.strip().split('\n')
        data_lines = [l for l in lines if not l.startswith('#') and l.strip()]

        if data_lines:
            first_data = data_lines[0].split()
            # Should have: Step Energy N_total Accept_rate InsAtt InsAcc DelAtt DelAcc TrnAtt TrnAcc RotAtt RotAcc
            assert len(first_data) >= 12, f"Data line has {len(first_data)} columns, expected at least 12"

            # Check that step is numeric
            step = int(first_data[0])
            assert step > 0, "Step should be positive"

            # Check that energy is numeric
            energy = float(first_data[1])
            # Energy can be negative or positive

            # Check that N_total is numeric
            n_total = int(first_data[2])
            assert n_total >= 0, "N_total should be non-negative"


class TestBasicGCMC:
    """Test basic GCMC simulation runs"""

    def test_minimal_water_simulation(self, tmp_path):
        """Test that a minimal water GCMC simulation completes successfully"""

        test_helper = TestMcMoveProb()
        inp_file = test_helper.create_minimal_system(tmp_path, [0.25, 0.25, 0.25, 0.25])

        # Run simulation
        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file)],
            capture_output=True,
            text=True,
            timeout=30,
            cwd=str(tmp_path)
        )

        # Should complete without error
        assert result.returncode == 0, f"Simulation failed: {result.stderr}"

        # Output files should exist (default prefix is "gcmc")
        top_files = list(tmp_path.glob("*.top"))
        pdb_files = list(tmp_path.glob("*.pdb"))
        assert len(top_files) > 0, f"Output TOP file not created. Files: {list(tmp_path.iterdir())}"
        assert len(pdb_files) > 0, f"Output PDB file not created. Files: {list(tmp_path.iterdir())}"

        # Read and validate statistics.dat file
        dat_files = list(tmp_path.glob("*statistics.dat"))
        assert len(dat_files) > 0, "Statistics DAT file not created"

        dat_file = dat_files[0]
        dat_content = dat_file.read_text()
        lines = dat_content.strip().split('\n')
        data_lines = [l for l in lines if not l.startswith('#') and l.strip()]

        assert len(data_lines) > 0, "No data lines in statistics DAT file"

        # Parse last data line (final statistics)
        last_line = data_lines[-1].split()
        assert len(last_line) >= 12, f"Statistics line has {len(last_line)} columns, expected at least 12"

        # Validate acceptance rate (should be between 0 and 1)
        accept_rate = float(last_line[3])
        assert 0.0 <= accept_rate <= 1.0, f"Acceptance rate {accept_rate} out of valid range [0, 1]"

        # Validate molecule count (should be non-negative)
        n_total = int(last_line[2])
        assert n_total >= 0, f"Molecule count {n_total} should be non-negative"

        # Validate energy (should be numeric, can be positive or negative)
        energy = float(last_line[1])
        assert isinstance(energy, float), "Energy should be numeric"

        # Validate attempt/accept counts are reasonable
        ins_att, ins_acc = int(last_line[4]), int(last_line[5])
        del_att, del_acc = int(last_line[6]), int(last_line[7])
        trn_att, trn_acc = int(last_line[8]), int(last_line[9])
        rot_att, rot_acc = int(last_line[10]), int(last_line[11])

        # Accepted should never exceed attempted
        assert ins_acc <= ins_att, f"Insert accepted {ins_acc} > attempted {ins_att}"
        assert del_acc <= del_att, f"Delete accepted {del_acc} > attempted {del_att}"
        assert trn_acc <= trn_att, f"Translate accepted {trn_acc} > attempted {trn_att}"
        assert rot_acc <= rot_att, f"Rotate accepted {rot_acc} > attempted {rot_att}"

        # Total attempts should be around mcsteps (100)
        total_attempts = ins_att + del_att + trn_att + rot_att
        assert total_attempts >= 50, f"Total attempts {total_attempts} seems too low for 100 MC steps"

    def test_multiple_fragments(self, tmp_path):
        """Test GCMC with multiple fragment types"""

        # Create files with actual methanol topology
        test_helper = TestMcMoveProb()
        pdb_file, top_file, atp_file, ff_file = test_helper._create_multi_fragment_files(tmp_path)

        # INP with two fragments
        inp_file = tmp_path / "test.inp"
        inp_content = f"""par:{ff_file}
atomtypes:{atp_file}
top:{top_file}
pdb:{pdb_file}
protitp:{top_file}

fragname: water methanol
fragconc: 1.0 1.0
fragmuex: -5.60 -3.0
mctime: 10 1

box_size: 25.0 25.0 25.0
cutoff: 10.0
temperature: 300
mcsteps: 2000
nprint: 200
eqsteps: 5

mc_move_prob: 0.3 0.3 0.2 0.2

op_top: {tmp_path}/output.top
op_pdb: {tmp_path}/output.pdb
"""
        inp_file.write_text(inp_content)

        result, accept_log, _ = _run_with_accept_log(inp_file, tmp_path)

        # Should complete
        assert result.returncode == 0, f"Multi-fragment simulation failed: {result.stderr}"

        assert accept_log.exists(), "Acceptance log missing"
        records = _load_accept_records(accept_log)

        seen_species = {str(r.get("species", "")).strip().upper() for r in records}
        assert {"WATER", "METHANOL"}.issubset(seen_species)

        counts = _move_counts(records)
        for move in MOVE_TYPES:
            assert counts[move] > 0, f"Expected at least one {move} attempt"

    def test_energy_calculation(self, tmp_path):
        """Test that energy is calculated and reported"""

        test_helper = TestMcMoveProb()
        inp_file = test_helper.create_minimal_system(tmp_path)

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file)],
            capture_output=True,
            text=True,
            timeout=30,
            cwd=str(tmp_path)
        )

        assert result.returncode == 0

        dat_files = list(tmp_path.glob("*statistics.dat"))
        assert dat_files, "Statistics DAT file not created"
        data_lines = [
            line for line in dat_files[0].read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ]
        assert data_lines, "No data lines in statistics DAT file"
        energy = float(data_lines[-1].split()[1])
        assert isinstance(energy, float)


class TestINPParameters:
    """Test INP parameter parsing"""

    def test_mctime_parsing(self, tmp_path):
        """Test that mctime parameter is parsed correctly"""

        # Create minimal files
        test_helper = TestMcMoveProb()
        pdb_file, top_file, atp_file, ff_file = test_helper._create_basic_files(tmp_path)

        inp_file = tmp_path / "test.inp"
        inp_content = f"""par:{ff_file}
atomtypes:{atp_file}
top:{top_file}
pdb:{pdb_file}
protitp:{top_file}

fragname: water methanol ethanol
fragconc: 55.0 1.0 0.5
fragmuex: -5.60 -3.0 -2.5
mctime: 10 1 0.5

box_size: 20.0 20.0 20.0
cutoff: 10.0
temperature: 300
mcsteps: 30
nprint: 15
eqsteps: 5

op_top: {tmp_path}/output.top
op_pdb: {tmp_path}/output.pdb
"""
        inp_file.write_text(inp_content)

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file)],
            capture_output=True,
            text=True,
            timeout=30,
            cwd=str(tmp_path)
        )

        # Should parse without error
        assert result.returncode == 0, f"Failed to parse mctime: {result.stderr}"

        # mctime should be mentioned or processed
        # (exact output depends on implementation)

    def test_invalid_mc_move_prob(self, tmp_path):
        """Test handling of invalid mc_move_prob (wrong number of values)"""

        test_helper = TestMcMoveProb()
        pdb_file, top_file, atp_file, ff_file = test_helper._create_basic_files(tmp_path)

        inp_file = tmp_path / "test.inp"
        inp_content = f"""par:{ff_file}
atomtypes:{atp_file}
top:{top_file}
pdb:{pdb_file}
protitp:{top_file}

fragname: water
fragconc: 55.0
fragmuex: 5.0

mc_move_prob: 0.5 0.5

box_size: 20.0 20.0 20.0
cutoff: 10.0
temperature: 300
mcsteps: 2000
nprint: 200
eqsteps: 0

op_top: {tmp_path}/output.top
op_pdb: {tmp_path}/output.pdb
"""
        inp_file.write_text(inp_content)
        result, accept_log, _ = _run_with_accept_log(inp_file, tmp_path)
        assert result.returncode == 0, f"Simulation failed: {result.stderr}"
        assert accept_log.exists(), "Acceptance log missing"
        params_json = tmp_path / "params.json"
        assert params_json.exists(), "dump-params output missing"

        params = _load_params(params_json)
        expected_cdf = _expected_cdf([0.25, 0.25, 0.25, 0.25])
        assert _fragment_cdf(params, "water") == pytest.approx(expected_cdf, abs=1e-6)


class TestAdvancedFeatures:
    """Test advanced GCMC features like cavity bias and CBMC"""

    def test_cavity_bias_parameter(self, tmp_path):
        """
        Test that use_cavity_bias parameter is recognized and enables cavity bias.
        Should print enabling message in setupEngine.

        Reference: GCMCSimulation.cpp:1213 (cavity bias setup)
        """

        test_helper = TestMcMoveProb()
        pdb_file, top_file, atp_file, ff_file = test_helper._create_basic_files(tmp_path)

        inp_file = tmp_path / "test.inp"
        inp_content = f"""par:{ff_file}
atomtypes:{atp_file}
top:{top_file}
pdb:{pdb_file}
protitp:{top_file}

fragname: water
fragconc: 55.0
fragmuex: -5.60

box_size: 20.0 20.0 20.0
cutoff: 10.0
temperature: 300
mcsteps: 50
nprint: 25
eqsteps: 5

mc_move_prob: 0.25 0.25 0.25 0.25
use_cavity_bias: yes

op_top: {tmp_path}/output.top
op_pdb: {tmp_path}/output.pdb
"""
        inp_file.write_text(inp_content)

        params_json = tmp_path / "params.json"
        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--dump-params", str(params_json)],
            capture_output=True,
            text=True,
            timeout=30,
            cwd=str(tmp_path)
        )

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"
        assert params_json.exists(), "dump-params output missing"
        params = json.loads(params_json.read_text())
        assert params.get("bias", {}).get("use_cavity_bias") is True

    def test_configurational_bias_parameter(self, tmp_path):
        """
        Test that use_conf_bias (CBMC) parameter is recognized and enables CBMC.
        Should print enabling message in setupEngine.

        Reference: GCMCSimulation.cpp:1152 (CBMC setup)
        """

        test_helper = TestMcMoveProb()
        pdb_file, top_file, atp_file, ff_file = test_helper._create_basic_files(tmp_path)

        inp_file = tmp_path / "test.inp"
        inp_content = f"""par:{ff_file}
atomtypes:{atp_file}
top:{top_file}
pdb:{pdb_file}
protitp:{top_file}

fragname: water
fragconc: 55.0
fragmuex: -5.60

box_size: 20.0 20.0 20.0
cutoff: 10.0
temperature: 300
mcsteps: 50
nprint: 25
eqsteps: 5

mc_move_prob: 0.25 0.25 0.25 0.25
use_conf_bias: yes
fragconf: 5

op_top: {tmp_path}/output.top
op_pdb: {tmp_path}/output.pdb
"""
        inp_file.write_text(inp_content)

        params_json = tmp_path / "params.json"
        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--dump-params", str(params_json)],
            capture_output=True,
            text=True,
            timeout=30,
            cwd=str(tmp_path)
        )

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"
        assert params_json.exists(), "dump-params output missing"
        params = json.loads(params_json.read_text())
        assert params.get("bias", {}).get("use_conf_bias") is True

    def test_combined_advanced_features(self, tmp_path):
        """
        Test using both cavity bias and CBMC together.
        Both features should be acknowledged in output if enabled.

        Reference: GCMCSimulation.cpp:1152, 1213 (setup messages)
        """

        test_helper = TestMcMoveProb()
        pdb_file, top_file, atp_file, ff_file = test_helper._create_basic_files(tmp_path)

        inp_file = tmp_path / "test.inp"
        inp_content = f"""par:{ff_file}
atomtypes:{atp_file}
top:{top_file}
pdb:{pdb_file}
protitp:{top_file}

fragname: water
fragconc: 55.0
fragmuex: -5.60

box_size: 20.0 20.0 20.0
cutoff: 10.0
temperature: 300
mcsteps: 50
nprint: 25
eqsteps: 5

mc_move_prob: 0.25 0.25 0.25 0.25
use_cavity_bias: yes
use_conf_bias: yes
fragconf: 5

op_top: {tmp_path}/output.top
op_pdb: {tmp_path}/output.pdb
"""
        inp_file.write_text(inp_content)

        params_json = tmp_path / "params.json"
        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--dump-params", str(params_json)],
            capture_output=True,
            text=True,
            timeout=30,
            cwd=str(tmp_path)
        )

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"
        assert params_json.exists(), "dump-params output missing"
        params = json.loads(params_json.read_text())
        assert params.get("bias", {}).get("use_cavity_bias") is True
        assert params.get("bias", {}).get("use_conf_bias") is True


class TestNumericalCorrectness:
    """Test numerical correctness of energy calculations and MC statistics"""

    @pytest.mark.skipif(not HAS_PYGCMC, reason="pygcmc bindings not available")
    def test_energy_baseline_comparison(self, tmp_path):
        """
        Verify O(N²) energy implementation by comparing:
        1. Energy from computeSystemEnergy (direct calculation)
        2. Energy from statistics.dat (simulation output)

        Reference: SimulationBasicBindings.cpp:70
        """
        test_helper = TestMcMoveProb()
        pdb_file, top_file, atp_file, ff_file = test_helper._create_basic_files(tmp_path)

        pdb_file.write_text(
            """CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1
ATOM      1  O   WAT     1      10.000  10.000  10.000  1.00  0.00
ATOM      2  H1  WAT     1      10.757  10.586  10.000  1.00  0.00
ATOM      3  H2  WAT     1       9.243  10.586  10.000  1.00  0.00
ATOM      4  O   WAT     2      15.000  10.000  10.000  1.00  0.00
ATOM      5  H1  WAT     2      15.757  10.586  10.000  1.00  0.00
ATOM      6  H2  WAT     2      14.243  10.586  10.000  1.00  0.00
END
"""
        )
        top_file.write_text(
            """[ defaults ]
1 2 yes 0.5 0.8333
[ atomtypes ]
O   8   15.9994  -0.834  A  3.15061e-01  6.36386e-01
H   1   1.008     0.417  A  0.00000e+00  0.00000e+00
[ moleculetype ]
WAT   3
[ atoms ]
1  O   1  WAT  O   1  -0.834  15.9994
2  H   1  WAT  H1  2   0.417   1.008
3  H   1  WAT  H2  3   0.417   1.008
[ system ]
Test
[ molecules ]
WAT  2
"""
        )

        inp_file = tmp_path / "energy.inp"
        inp_file.write_text(
            f"""par:{ff_file}
atomtypes:{atp_file}
top:{top_file}
pdb:{pdb_file}
protitp:{top_file}

fragname: WAT
fragconc: 55.0
fragmuex: -5.60

box_size: 20.0 20.0 20.0
cutoff: 10.0
temperature: 300
mcsteps: 1
nprint: 100
eqsteps: 0

mc_move_prob: 0.0 0.0 1.0 0.0

op_top: {tmp_path}/output.top
op_pdb: {tmp_path}/output.pdb
"""
        )

        params_json = tmp_path / "params.json"
        out_prefix = tmp_path / "energy"

        result = subprocess.run(
            [
                str(GCMC_CPU_PATH),
                "--inp",
                str(inp_file),
                "--seed",
                "123",
                "--stats-interval",
                "1",
                "--dump-params",
                str(params_json),
                "--prefix",
                str(out_prefix),
            ],
            capture_output=True,
            text=True,
            timeout=30,
            cwd=str(tmp_path)
        )

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"
        assert params_json.exists(), "dump-params output missing"

        # Read statistics.dat to get energy from simulation
        dat_files = list(tmp_path.glob("*statistics.dat"))
        assert len(dat_files) > 0, "Statistics DAT file not created"

        dat_content = dat_files[0].read_text()
        lines = dat_content.strip().split('\n')
        data_lines = [l for l in lines if not l.startswith('#') and l.strip()]
        assert data_lines, "Statistics DAT file contains no data rows"

        # Get energy from the last data line (final snapshot)
        last_line = data_lines[-1].split()
        sim_energy = float(last_line[1])  # Energy column

        params = _load_params(params_json)
        space = params.get("space", {})

        final_pdb = Path(str(out_prefix) + "_final.pdb")
        assert final_pdb.exists(), "Final PDB output missing"

        structure = pygcmc.PDBParser.parse_file(str(final_pdb))
        topology = pygcmc.TOPParser.parse_file(str(top_file))
        molecular = pygcmc.MolecularSystem().combine(structure, topology)

        mc_system = pygcmc.MonteCarloSystem()
        info = pygcmc.MCInfo()
        info.max_residues = 100
        info.max_atoms = 1000
        mc_system.initialize(info)
        mc_system.initialize_from_molecular(molecular)

        state = mc_system.get_state_mutable()
        if "cutoff_nm" in space:
            state.info.cutoff = float(space["cutoff_nm"])
        if "box_size_nm" in space:
            state.info.box = [float(x) for x in space["box_size_nm"]]

        atomtypes = _parse_itp_atomtypes(ff_file)
        for name in sorted(atomtypes):
            state.atomTypes.get_or_add_type(name)

        num_types = len(state.atomTypes.atomTypes)
        state.forcefield.numTotalTypes = num_types
        state.forcefield.numMovementTypes = num_types
        state.forcefield.mixingRule = pygcmc.MCForceField.MixingRule.LorentzBerthelot
        sigma_types = [0.0] * num_types
        eps_types = [0.0] * num_types
        for name in sorted(atomtypes):
            sigma, eps = atomtypes[name]
            idx = state.atomTypes.get_or_add_type(name)
            if 0 <= idx < num_types:
                sigma_types[idx] = float(sigma)
                eps_types[idx] = float(eps)

        state.forcefield.ljSigmaType = sigma_types
        state.forcefield.ljEpsType = eps_types
        state.forcefield.rebuildLJMatrix()

        cpu_platform.computeSystemEnergyPBCCutoff(state)
        computed_energy = pygcmc.getTotalEnergyUniquePairs(state, pygcmc.EnergyMethod.DIRECT)

        assert not np.isnan(sim_energy), "Simulation energy is NaN"
        assert not np.isinf(sim_energy), "Simulation energy is Inf"
        assert not np.isnan(computed_energy), "Computed energy is NaN"
        assert not np.isinf(computed_energy), "Computed energy is Inf"

        assert computed_energy == pytest.approx(sim_energy, rel=5e-2, abs=2e-1)

    def test_mc_move_probability_distribution(self, tmp_path):
        """
        Test MC move sampling with fixed seed for reproducibility.
        Verify that InsAtt/DelAtt/TrnAtt/RotAtt ratios are reasonable.

        Note: Actual ratios may differ from mc_move_prob due to:
        1. System state constraints (can't delete/translate/rotate if N=0)
        2. Fragment selection when multiple fragments exist
        3. Move type selection may be per-fragment

        Reference: GCMCSimulation.cpp:1757 (writeStatisticsDAT)
        """
        test_helper = TestMcMoveProb()
        pdb_file, top_file, atp_file, ff_file = test_helper._create_basic_files(tmp_path)

        # Use equal probabilities to avoid state-dependent biases
        mc_probs = [0.25, 0.25, 0.25, 0.25]  # Insert, Delete, Translate, Rotate

        inp_file = tmp_path / "test.inp"
        inp_content = f"""par:{ff_file}
atomtypes:{atp_file}
top:{top_file}
pdb:{pdb_file}
protitp:{top_file}

fragname: water
fragconc: 55.0
fragmuex: -5.60

box_size: 20.0 20.0 20.0
cutoff: 10.0
temperature: 300
mcsteps: 1000
nprint: 200
eqsteps: 0
seed: 42

mc_move_prob: {' '.join(map(str, mc_probs))}

op_top: {tmp_path}/output.top
op_pdb: {tmp_path}/output.pdb
"""
        inp_file.write_text(inp_content)

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file)],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(tmp_path)
        )

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"

        # Read final statistics
        dat_files = list(tmp_path.glob("*statistics.dat"))
        assert len(dat_files) > 0, "Statistics DAT file not created"

        dat_content = dat_files[0].read_text()
        lines = dat_content.strip().split('\n')
        data_lines = [l for l in lines if not l.startswith('#') and l.strip()]

        # Parse last line for cumulative statistics
        last_line = data_lines[-1].split()
        ins_att = int(last_line[4])
        del_att = int(last_line[6])
        trn_att = int(last_line[8])
        rot_att = int(last_line[10])

        total_att = ins_att + del_att + trn_att + rot_att
        assert total_att > 0, "No MC moves were attempted"

        # Print diagnostic information
        print(f"\nMove attempts: Ins={ins_att}, Del={del_att}, Trn={trn_att}, Rot={rot_att}")
        print(f"Total attempts: {total_att}")

        # Calculate observed ratios
        obs_ratios = np.array([ins_att, del_att, trn_att, rot_att]) / total_att
        expected_ratios = np.array(mc_probs)

        print(f"Observed ratios: {obs_ratios}")
        print(f"Expected ratios: {expected_ratios}")

        # Verify basic properties:
        # 1. All move types should be attempted at least once
        assert ins_att > 0, "Insert moves were never attempted"
        # Note: Delete/Translate/Rotate may be 0 if system never had molecules

        # 2. Total attempts should match mcsteps
        # Each MC step may involve multiple fragment moves, so total >= mcsteps
        assert total_att >= 1000, f"Total attempts {total_att} < mcsteps (1000)"

        # 3. Ratios should be reasonable (not all 100% one type)
        # With mc_move_prob = [0.25, 0.25, 0.25, 0.25], no single move should dominate completely
        for i, ratio in enumerate(obs_ratios):
            assert ratio < 0.95, f"Move type {i} ratio {ratio:.3f} is too dominant (>95%)"

        # Note: Reproducibility with seed is NOT tested here because:
        # - The seed parameter may not be fully implemented
        # - There may be non-deterministic factors in the simulation
        # - This would require investigation of the RNG implementation

    def test_multi_fragment_weight_verification(self, tmp_path):
        """
        Verify mctime weights affect fragment selection probability.
        Check that Move probability cdf (per fragment) appears for all fragments
        and that probabilities sum to 1.0.

        Reference: InpParserGCMC.cpp:41 (mctime), GCMCSimulation.cpp:1006 (CDF print)
        """
        test_helper = TestMcMoveProb()
        pdb_file, top_file, atp_file, ff_file = test_helper._create_multi_fragment_files(tmp_path)

        # Set different mctime weights: water gets 5x more selection probability than methanol
        inp_file = tmp_path / "test.inp"
        inp_content = f"""par:{ff_file}
atomtypes:{atp_file}
top:{top_file}
pdb:{pdb_file}
protitp:{top_file}

fragname: water methanol
fragconc: 55.0 1.0
fragmuex: -5.60 -3.0
mctime: 5.0 1.0

box_size: 20.0 20.0 20.0
cutoff: 10.0
temperature: 300
mcsteps: 2000
nprint: 200
eqsteps: 5

mc_move_prob: 0.25 0.25 0.25 0.25

op_top: {tmp_path}/output.top
op_pdb: {tmp_path}/output.pdb
"""
        inp_file.write_text(inp_content)

        result, accept_log, _ = _run_with_accept_log(inp_file, tmp_path)

        assert result.returncode == 0, f"Multi-fragment simulation failed: {result.stderr}"
        assert accept_log.exists(), "Acceptance log missing"
        params_json = tmp_path / "params.json"
        assert params_json.exists(), "dump-params output missing"

        params = _load_params(params_json)
        names = params.get("fragment", {}).get("names", [])
        probs = params.get("fragment", {}).get("selection_prob", [])
        assert names and probs, "Fragment selection probabilities missing in params"

        idx_water = _fragment_index(params, "water")
        idx_meoh = _fragment_index(params, "methanol")
        expected_water_prob = 5.0 / (5.0 + 1.0)
        observed_water_prob = float(probs[idx_water])
        observed_meoh_prob = float(probs[idx_meoh])
        assert observed_water_prob == pytest.approx(expected_water_prob, abs=1e-6)
        assert observed_meoh_prob == pytest.approx(1.0 - expected_water_prob, abs=1e-6)


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v", "-s"])
