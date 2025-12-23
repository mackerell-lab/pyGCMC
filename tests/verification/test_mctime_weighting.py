"""
Test mctime weighting mechanism.

Verifies P1: mctime采样统计
- Fragment selection weights are normalized in dump-params
- JSONL insertion records cover all fragments (avoid silent single-fragment sampling)
"""

import pytest
import subprocess
import json
from pathlib import Path
from typing import List, Dict, Any
from collections import defaultdict


# Path to gcmc_cpu executable
GCMC_CPU_PATH = Path(__file__).parent.parent.parent / "build" / "bin" / "gcmc_cpu"


# ============================================================================
# JSONL Utilities (inline to avoid cross-directory conftest imports)
# ============================================================================

def read_jsonl(filepath: Path) -> List[Dict[str, Any]]:
    """Read JSONL file and return list of records."""
    records = []
    if not filepath.exists():
        return records
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError as e:
                print(f"Warning: Failed to parse line: {line[:50]}... Error: {e}")
    return records


def filter_by_move(records: List[Dict[str, Any]], move_type: str) -> List[Dict[str, Any]]:
    """Filter records by move type."""
    return [r for r in records if r.get('move') == move_type]


def count_by_species(records: List[Dict[str, Any]]) -> Dict[str, int]:
    """Count total attempts per species."""
    counts = defaultdict(int)
    for record in records:
        species = str(record.get('species', 'unknown')).lower()
        counts[species] += 1
    return dict(counts)

def load_params(filepath: Path) -> Dict[str, Any]:
    """Read JSON params file."""
    return json.loads(filepath.read_text())


def fragment_index(params: Dict[str, Any], name: str) -> int:
    names = params.get("fragment", {}).get("names", [])
    target = name.strip().upper()
    for idx, frag_name in enumerate(names):
        if str(frag_name).strip().upper() == target:
            return idx
    if not names:
        return 0
    raise AssertionError(f"Fragment name {name} not found in params: {names}")


def selection_prob(params: Dict[str, Any], name: str) -> float:
    idx = fragment_index(params, name)
    probs = params.get("fragment", {}).get("selection_prob", [])
    if idx >= len(probs):
        raise AssertionError(f"selection_prob missing index {idx} for fragment {name}")
    return float(probs[idx])


class TestMctimeWeighting:
    """Test mctime-based fragment selection weighting"""

    @staticmethod
    def create_two_fragment_system(
        tmpdir,
        mctime_weights,
        mcsteps=5000,
        seed=42,
        *,
        eqsteps: int = 1000,
        nprint: int = 2000,
        box_size: tuple[float, float, float] = (30.0, 30.0, 30.0),
        cutoff: float = 14.0,
        mc_move_prob: tuple[float, float, float, float] = (0.5, 0.5, 0.0, 0.0),
        fragmuex: tuple[float, float] = (-5.60, -6.50),
        fragconc: tuple[float, float] = (55.0, 24.0),
    ):
        """Create system with two fragments (water and methanol) with specified mctime weights"""

        # Initial PDB with one water molecule
        pdb_file = tmpdir / "test.pdb"
        pdb_content = """CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
ATOM      1  O   WAT     1      15.000  15.000  15.000  1.00  0.00
ATOM      2  H1  WAT     1      15.757  15.586  15.000  1.00  0.00
ATOM      3  H2  WAT     1      14.243  15.586  15.000  1.00  0.00
END
"""
        pdb_file.write_text(pdb_content)

        # TOP file with water and methanol
        top_file = tmpdir / "test.top"
        top_content = """[ defaults ]
1 2 yes 0.5 0.8333

[ atomtypes ]
O   8   15.9994  -0.834  A  3.15061e-01  6.36386e-01
H   1   1.008     0.417  A  0.00000e+00  0.00000e+00
C   6   12.011    0.145  A  3.50000e-01  2.76144e-01
OM  8   15.9994  -0.683  A  3.07000e-01  7.11280e-01

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
1  C   1  MEOH  C   1   0.145  12.011
2  OM  1  MEOH  O   2  -0.683  15.9994
3  H   1  MEOH  H1  3   0.040   1.008
4  H   1  MEOH  H2  4   0.040   1.008
5  H   1  MEOH  H3  5   0.040   1.008
6  H   1  MEOH  HO  6   0.418   1.008

[ bonds ]
1  2  1  0.1430  267776.0
1  3  1  0.1090  284512.0
1  4  1  0.1090  284512.0
1  5  1  0.1090  284512.0
2  6  1  0.0945  462750.4

[ angles ]
2  1  3  1  109.5  292.88
2  1  4  1  109.5  292.88
2  1  5  1  109.5  292.88
3  1  4  1  109.5  276.14
3  1  5  1  109.5  276.14
4  1  5  1  109.5  276.14
1  2  6  1  108.5  460.24

[ system ]
Water-Methanol System

[ molecules ]
WAT  1
"""
        top_file.write_text(top_content)

        # Atomtypes
        atp_file = tmpdir / "atomtypes.atp"
        atp_file.write_text("O   15.9994\\nH    1.008\\nC   12.011\\nOM  15.9994\\n")

        # Force field
        ff_file = tmpdir / "ffnonbonded.itp"
        ff_file.write_text("""[ atomtypes ]
O   8   15.9994  -0.834  A  3.15061e-01  6.36386e-01
H   1   1.008     0.417  A  0.00000e+00  0.00000e+00
C   6   12.011    0.145  A  3.50000e-01  2.76144e-01
OM  8   15.9994  -0.683  A  3.07000e-01  7.11280e-01
""")

        # INP file with two fragments and mctime weights
        mctime_str = " ".join(map(str, mctime_weights))
        mc_move_prob_str = " ".join(map(str, mc_move_prob))

        inp_file = tmpdir / "test.inp"
        inp_content = f"""par:{ff_file}
atomtypes:{atp_file}
top:{top_file}
pdb:{pdb_file}
protitp:{top_file}

fragname: water
fragconc: {fragconc[0]}
fragmuex: {fragmuex[0]}

fragname: methanol
fragconc: {fragconc[1]}
fragmuex: {fragmuex[1]}

mctime: {mctime_str}

box_size: {box_size[0]} {box_size[1]} {box_size[2]}
cutoff: {cutoff}
temperature: 300
mcsteps: {mcsteps}
nprint: {nprint}
eqsteps: {eqsteps}

moves_per_step: 1
mc_move_prob: {mc_move_prob_str}

seed: {seed}
"""
        inp_file.write_text(inp_content)

        return inp_file

    @staticmethod
    def parse_fragment_selection_counts(result_dir):
        """
        Extract fragment selection counts from gcmc_final.txt.

        Fragment selection count = insertAttempts ONLY
        (Deletion doesn't use mctime weighting - it selects from active fragments)
        """
        import re
        counts = {}

        # Look for gcmc_final.txt
        result_file = result_dir / "gcmc_final.txt"

        if not result_file.exists():
            return counts

        with open(result_file, 'r') as f:
            content = f.read()

        # Parse per-fragment statistics
        # Format:
        #   water:
        #     Final count: 30
        #     Density: ...
        #     Insert attempts: 37
        #     Delete attempts: 50138
        fragments = re.findall(r'^\s{2}(\w+):\s*$.*?Insert attempts:\s*(\d+)',
                              content, re.MULTILINE | re.DOTALL)

        for frag_name, insert_att in fragments:
            frag_name = frag_name.lower()
            counts[frag_name] = int(insert_att)

        return counts

    def test_mctime_weights_output(self, tmp_path):
        """
        P1验收测试：Fragment weights输出验证

        验证：
        - Fragment weights被正确归一化并反映在 dump-params 中
        - mctime: [5, 1] 应该产生 water=0.8333, methanol=0.1667
        """
        print(f"\\n=== mctime Weights Output Test ===")

        mctime_weights = [5.0, 1.0]
        expected_weights = {
            'water': 5.0 / 6.0,      # 0.8333...
            'methanol': 1.0 / 6.0    # 0.1667...
        }

        print(f"mctime weights: {mctime_weights}")
        print(f"Expected normalized: water={expected_weights['water']:.4f}, methanol={expected_weights['methanol']:.4f}")

        inp_file = self.create_two_fragment_system(
            tmp_path,
            mctime_weights=mctime_weights,
            mcsteps=2000,  # Short run but enough samples for weighting
            seed=42
        )

        params_json = tmp_path / "params.json"
        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--seed", "42", "--dump-params", str(params_json)],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(tmp_path)
        )

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"
        assert params_json.exists(), "dump-params output missing"

        params = load_params(params_json)
        observed = {
            "water": selection_prob(params, "water"),
            "methanol": selection_prob(params, "methanol"),
        }

        print(f"\\nObserved weights: {observed}")

        water_error = abs(observed["water"] - expected_weights["water"])
        methanol_error = abs(observed["methanol"] - expected_weights["methanol"])

        print(f"water: expected={expected_weights['water']:.4f}, actual={observed['water']:.4f}, error={water_error:.6f}")
        print(f"methanol: expected={expected_weights['methanol']:.4f}, actual={observed['methanol']:.4f}, error={methanol_error:.6f}")

        assert water_error < 1e-6, f"water weight error {water_error:.6f} exceeds 1e-6"
        assert methanol_error < 1e-6, f"methanol weight error {methanol_error:.6f} exceeds 1e-6"

        print(f"\\n✅ Fragment weights correctly normalized and output")

    # Historical note:
    # - test_mctime_basic_weighting (removed 2024-10-20):
    #   Replaced by test_mctime_sampling_distribution_jsonl (uses JSONL logs, more robust)
    # - test_mctime_chi_square (removed 2024-10-20):
    #   Replaced by test_mctime_sampling_distribution_jsonl (uses JSONL logs, more robust)
    # Archived tests available in: tmp/test_mctime_weighting_deprecated.py.archived

    def test_mctime_sampling_distribution_jsonl(self, tmp_path):
        """
        P1验收测试：mctime采样分布验证（使用JSONL acceptance log）

        验证：
        - 使用 --dump-accept 导出JSONL acceptance log
        - 使用 --dump-params 获取归一化的 fragment selection 概率
        - JSONL 至少覆盖两个片段的 insertion 记录
        """
        print(f"\n=== mctime Sampling Distribution Test (JSONL) ===")

        mctime_weights = [5.0, 1.0]  # water:methanol = 5:1
        total_weight = sum(mctime_weights)
        expected_probs = [w / total_weight for w in mctime_weights]

        print(f"mctime weights: {mctime_weights}")
        print(f"Expected probabilities: water={expected_probs[0]:.3f}, methanol={expected_probs[1]:.3f}")

        inp_file = self.create_two_fragment_system(
            tmp_path,
            mctime_weights=mctime_weights,
            # This test verifies fragment *selection* weighting (mctime), not equilibrium.
            # Keep the run short to avoid expensive growth to bulk-like densities.
            mcsteps=2000,
            seed=42,
            eqsteps=0,
            nprint=2000,
            box_size=(20.0, 20.0, 20.0),
            cutoff=8.0,
            mc_move_prob=(1.0, 0.0, 0.0, 0.0),
            fragmuex=(-50.0, -50.0),
        )

        accept_log = tmp_path / "accept.jsonl"
        params_json = tmp_path / "params.json"
        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file),
             "--dump-accept", str(accept_log), "--dump-params", str(params_json), "--seed", "42"],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(tmp_path)
        )

        if result.returncode != 0:
            print(f"\n❌ Simulation failed with return code {result.returncode}")
            print(f"STDERR:\n{result.stderr}")
        assert result.returncode == 0, f"Simulation failed: {result.stderr}"

        assert accept_log.exists(), f"Acceptance log not generated: {accept_log}"
        assert params_json.exists(), "dump-params output missing"

        params = load_params(params_json)
        observed_probs = [
            selection_prob(params, "water"),
            selection_prob(params, "methanol"),
        ]
        print(f"Observed selection probs: {observed_probs}")
        assert observed_probs[0] == pytest.approx(expected_probs[0], abs=1e-6)
        assert observed_probs[1] == pytest.approx(expected_probs[1], abs=1e-6)

        records = read_jsonl(accept_log)
        insertion_records = filter_by_move(records, 'insertion')
        counts = count_by_species(insertion_records)
        total_ins = sum(counts.values())
        assert total_ins >= 200, f"Too few insertion attempts recorded: {total_ins}"
        assert counts.get('water', 0) > 0, "No water insertion attempts recorded"
        assert counts.get('methanol', 0) > 0, "No methanol insertion attempts recorded"

        # Distribution check (avoid overfitting exact RNG stream):
        # observed insertion attempts should be consistent with the normalized mctime weights.
        frac_methanol = counts["methanol"] / total_ins
        assert frac_methanol == pytest.approx(expected_probs[1], abs=0.05)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
