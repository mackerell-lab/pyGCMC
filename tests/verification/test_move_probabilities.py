"""
Test mc_move_prob parsing and CDF generation.

Verifies P0.2: mc_move_prob生效验证
- Global mc_move_prob is parsed correctly
- Single-element vector is broadcasted to all fragments
- CDF is correctly normalized
"""

import json
import pytest
import subprocess
from pathlib import Path

# Path to gcmc_cpu executable
GCMC_CPU_PATH = Path(__file__).parent.parent.parent / "build" / "bin" / "gcmc_cpu"


MOVE_TYPES = ("insertion", "deletion", "translation", "rotation")


def _load_params(path: Path) -> dict:
    return json.loads(path.read_text())


def _load_accept_records(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


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


def _requested_move_counts(records: list[dict]) -> dict[str, int]:
    counts = {move: 0 for move in MOVE_TYPES}
    for rec in records:
        move = str(rec.get("requestedMove", "")).strip().lower()
        if move in counts:
            counts[move] += 1
    return counts


def _requested_move_fractions(counts: dict[str, int]) -> dict[str, float]:
    total = sum(counts.values())
    if total <= 0:
        raise AssertionError("No requestedMove records found in acceptance log")
    return {move: counts[move] / total for move in MOVE_TYPES}


def _run_with_params(inp_file: Path, tmp_path: Path) -> tuple[subprocess.CompletedProcess, Path]:
    params_json = tmp_path / "params.json"
    result = subprocess.run(
        [
            str(GCMC_CPU_PATH),
            "--inp",
            str(inp_file),
            "--seed",
            "42",
            "--prefix",
            str(tmp_path / "out"),
            "--dump-params",
            str(params_json),
        ],
        capture_output=True,
        text=True,
        timeout=60,
        cwd=str(tmp_path),
    )
    return result, params_json


class TestMoveProbabilities:
    """Test move probability parsing and distribution via acceptance log"""

    @staticmethod
    def create_minimal_system(
        tmpdir,
        mc_move_prob_line=None,
        fragment_specific=False,
        mcsteps=2000,
        fragmuex_value=-5.60,
    ):
        """Create minimal water system with mc_move_prob setting"""

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

        # Atomtypes
        atp_file = tmpdir / "atomtypes.atp"
        atp_file.write_text("O   15.9994\nH    1.008\n")

        # Force field
        ff_file = tmpdir / "ffnonbonded.itp"
        ff_file.write_text("""[ atomtypes ]
O   8   15.9994  -0.834  A  3.15061e-01  6.36386e-01
H   1   1.008     0.417  A  0.00000e+00  0.00000e+00
""")

        # INP file
        inp_file = tmpdir / "test.inp"
        inp_content = f"""par:{ff_file}
atomtypes:{atp_file}
top:{top_file}
pdb:{pdb_file}
protitp:{top_file}

fragname: water
fragconc: 55.0
fragmuex: {fragmuex_value}

box_size: 20.0 20.0 20.0
cutoff: 10.0
temperature: 300
mcsteps: {mcsteps}
nprint: 200
eqsteps: 0

"""
        # Add mc_move_prob or fragment-specific lines
        if mc_move_prob_line:
            inp_content += f"mc_move_prob: {mc_move_prob_line}\n"

        if fragment_specific:
            # Add per-fragment attempt probabilities
            inp_content += "attempt_prob_ins: 5.0\n"
            inp_content += "attempt_prob_del: 4.0\n"
            inp_content += "attempt_prob_trn: 3.0\n"
            inp_content += "attempt_prob_rot: 2.0\n"

        inp_content += f"""
seed: 42

op_top: {tmpdir}/output.top
op_pdb: {tmpdir}/output.pdb
"""
        inp_file.write_text(inp_content)

        return inp_file

    def test_global_mc_move_prob_cdf(self, tmp_path):
        """
        P0.2验收测试：全局mc_move_prob设置，验证CDF正确归一化

        设置 mc_move_prob: 1 2 3 4
        期望 CDF: [0.1, 0.3, 0.6, 1.0]
        """
        inp_file = self.create_minimal_system(tmp_path, mc_move_prob_line="1 2 3 4", mcsteps=2000)
        result, params_json = _run_with_params(inp_file, tmp_path)

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"
        assert params_json.exists(), "dump-params output missing"

        params = _load_params(params_json)
        actual_cdf = _fragment_cdf(params, "water")
        expected_cdf = _expected_cdf([1.0, 2.0, 3.0, 4.0])
        assert actual_cdf == pytest.approx(expected_cdf, abs=1e-6)

    def test_mc_move_prob_broadcasting(self, tmp_path):
        """
        验证单元素向量广播到所有fragments

        即使只有一个fragment，广播逻辑也应该工作
        """
        inp_file = self.create_minimal_system(tmp_path, mc_move_prob_line="2 3 4 1", mcsteps=2000)
        result, params_json = _run_with_params(inp_file, tmp_path)

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"
        assert params_json.exists(), "dump-params output missing"

        params = _load_params(params_json)
        actual_cdf = _fragment_cdf(params, "water")
        expected_cdf = _expected_cdf([2.0, 3.0, 4.0, 1.0])
        assert actual_cdf == pytest.approx(expected_cdf, abs=1e-6)

    def test_fragment_specific_overrides_global(self, tmp_path):
        """
        验证per-fragment显式设置优先于全局mc_move_prob

        设置全局 mc_move_prob: 1 1 1 1 (应该被忽略)
        设置per-fragment: 5 4 3 2
        期望 CDF 来自 per-fragment 设置
        """
        # Create system with both global and per-fragment settings
        pdb_file = tmp_path / "test.pdb"
        pdb_content = """CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1
ATOM      1  O   WAT     1      10.000  10.000  10.000  1.00  0.00
ATOM      2  H1  WAT     1      10.757  10.586  10.000  1.00  0.00
ATOM      3  H2  WAT     1       9.243  10.586  10.000  1.00  0.00
END
"""
        pdb_file.write_text(pdb_content)

        top_file = tmp_path / "test.top"
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

        atp_file = tmp_path / "atomtypes.atp"
        atp_file.write_text("O   15.9994\nH    1.008\n")

        ff_file = tmp_path / "ffnonbonded.itp"
        ff_file.write_text("""[ atomtypes ]
O   8   15.9994  -0.834  A  3.15061e-01  6.36386e-01
H   1   1.008     0.417  A  0.00000e+00  0.00000e+00
""")

        inp_file = tmp_path / "test.inp"
        # Global setting should be overridden by per-fragment
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
mcsteps: 2000
nprint: 200
eqsteps: 0

mc_move_prob: 1 1 1 1
attempt_prob_ins: 5.0
attempt_prob_del: 4.0
attempt_prob_trn: 3.0
attempt_prob_rot: 2.0

seed: 42

op_top: {tmp_path}/output.top
op_pdb: {tmp_path}/output.pdb
"""
        inp_file.write_text(inp_content)

        result, params_json = _run_with_params(inp_file, tmp_path)

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"
        assert params_json.exists(), "dump-params output missing"

        params = _load_params(params_json)
        actual_cdf = _fragment_cdf(params, "water")
        expected_cdf = _expected_cdf([5.0, 4.0, 3.0, 2.0])
        assert actual_cdf == pytest.approx(expected_cdf, abs=1e-6)

    def test_invalid_mc_move_prob_fallback(self, tmp_path):
        """
        验证错误的mc_move_prob参数数量时的回退行为

        提供少于4个参数，应该警告并使用默认值
        """
        inp_file = self.create_minimal_system(tmp_path, mc_move_prob_line="1 2", mcsteps=2000)  # Only 2 values
        result, params_json = _run_with_params(inp_file, tmp_path)

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"
        assert params_json.exists(), "dump-params output missing"

        params = _load_params(params_json)
        actual_cdf = _fragment_cdf(params, "water")
        expected_cdf = _expected_cdf([0.25, 0.25, 0.25, 0.25])
        assert actual_cdf == pytest.approx(expected_cdf, abs=1e-6)

    def test_requested_move_distribution_matches_cdf(self, tmp_path):
        """
        Verify requestedMove sampling distribution matches configured move CDF.

        This checks the actual move selection BEFORE fallback adjustments.
        """
        inp_file = self.create_minimal_system(
            tmp_path,
            mc_move_prob_line="4 3 2 1",
            mcsteps=6000,
            fragmuex_value=5.0,
        )
        inp_file.write_text(
            inp_file.read_text().replace(
                "box_size: 20.0 20.0 20.0",
                "box_size: 50.0 50.0 50.0",
            )
        )
        accept_log = tmp_path / "accept.jsonl"
        params_json = tmp_path / "params.json"

        result = subprocess.run(
            [
                str(GCMC_CPU_PATH),
                "--inp",
                str(inp_file),
                "--seed",
                "4242",
                "--prefix",
                str(tmp_path / "out"),
                "--dump-accept",
                str(accept_log),
                "--dump-params",
                str(params_json),
            ],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(tmp_path),
        )

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"
        assert accept_log.exists(), "Acceptance log not generated"
        assert params_json.exists(), "dump-params output missing"

        params = _load_params(params_json)
        cdf = _fragment_cdf(params, "water")
        expected = {
            "insertion": cdf[0],
            "deletion": cdf[1] - cdf[0],
            "translation": cdf[2] - cdf[1],
            "rotation": 1.0 - cdf[2],
        }

        records = _load_accept_records(accept_log)
        counts = _requested_move_counts(records)
        fractions = _requested_move_fractions(counts)

        for move in MOVE_TYPES:
            assert counts[move] > 0, f"No requestedMove records for {move}"
            assert abs(fractions[move] - expected[move]) < 0.05, (
                f"{move} fraction {fractions[move]:.3f} differs from expected {expected[move]:.3f}"
            )


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
