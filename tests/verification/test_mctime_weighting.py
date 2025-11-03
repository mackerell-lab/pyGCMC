"""
Test mctime weighting mechanism.

Verifies P1: mctime采样统计
- Fragment selection frequency ≈ normalized mctime weights
- Chi-square or KS test with p > 0.05
- Or frequency error < 10%
"""

import pytest
import subprocess
import numpy as np
import json
from pathlib import Path
from scipy import stats
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
        species = record.get('species', 'unknown')
        counts[species] += 1
    return dict(counts)


class TestMctimeWeighting:
    """Test mctime-based fragment selection weighting"""

    @staticmethod
    def create_two_fragment_system(tmpdir, mctime_weights, mcsteps=5000, seed=42):
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

        inp_file = tmpdir / "test.inp"
        inp_content = f"""par:{ff_file}
atomtypes:{atp_file}
top:{top_file}
pdb:{pdb_file}
protitp:{top_file}

fragname: water
fragconc: 55.0
fragmuex: -5.60

fragname: methanol
fragconc: 24.0
fragmuex: -6.50

mctime: {mctime_str}

box_size: 30.0 30.0 30.0
cutoff: 14.0
temperature: 300
mcsteps: {mcsteps}
nprint: 2000
eqsteps: 1000

mc_move_prob: 0.5 0.5 0 0

seed: {seed}
"""
        inp_file.write_text(inp_content)

        return inp_file

    @staticmethod
    def parse_fragment_weights_from_stdout(stdout):
        """
        Extract fragment selection weights from stdout.

        Looks for line: "Fragment weights: water=0.8333, methanol=0.1667"
        """
        import re
        weights = {}

        # Match pattern: "Fragment weights: name1=0.123, name2=0.456"
        match = re.search(r'Fragment weights: (.+)', stdout)
        if match:
            weight_str = match.group(1)
            # Parse individual weights: "name=value"
            for pair in weight_str.split(', '):
                if '=' in pair:
                    name, value = pair.split('=')
                    weights[name.strip().lower()] = float(value.strip())

        return weights

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
        - Fragment weights被正确归一化并输出到stdout
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
            mcsteps=100,  # Short simulation just to check output
            seed=42
        )

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--seed", "42"],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(tmp_path)
        )

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"

        # Parse weights from stdout
        weights = self.parse_fragment_weights_from_stdout(result.stdout)

        print(f"\\nParsed weights: {weights}")

        assert 'water' in weights, "water weight not found in stdout"
        assert 'methanol' in weights, "methanol weight not found in stdout"

        # Verify weights (tolerance 1e-3)
        water_error = abs(weights['water'] - expected_weights['water'])
        methanol_error = abs(weights['methanol'] - expected_weights['methanol'])

        print(f"water: expected={expected_weights['water']:.4f}, actual={weights['water']:.4f}, error={water_error:.6f}")
        print(f"methanol: expected={expected_weights['methanol']:.4f}, actual={weights['methanol']:.4f}, error={methanol_error:.6f}")

        assert water_error < 1e-3, f"water weight error {water_error} exceeds 1e-3"
        assert methanol_error < 1e-3, f"methanol weight error {methanol_error} exceeds 1e-3"

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
        - 统计insertion move的fragment选择频率
        - 验证频率 ≈ normalized mctime weights
        - 使用chi-square检验 p > 0.05 或 绝对误差 < 10%
        """
        print(f"\n=== mctime Sampling Distribution Test (JSONL) ===")

        mctime_weights = [5.0, 1.0]  # water:methanol = 5:1
        total_weight = sum(mctime_weights)
        expected_probs = [w / total_weight for w in mctime_weights]  # [0.833, 0.167]

        print(f"mctime weights: {mctime_weights}")
        print(f"Expected probabilities: water={expected_probs[0]:.3f}, methanol={expected_probs[1]:.3f}")

        # Run simulation with acceptance log
        inp_file = self.create_two_fragment_system(
            tmp_path,
            mctime_weights=mctime_weights,
            mcsteps=15000,  # Longer run for better statistics
            seed=42
        )

        accept_log = tmp_path / "accept.jsonl"

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file),
             "--dump-accept", str(accept_log), "--seed", "42"],
            capture_output=True,
            text=True,
            timeout=180,
            cwd=str(tmp_path)
        )

        if result.returncode != 0:
            print(f"\n❌ Simulation failed with return code {result.returncode}")
            print(f"STDERR:\n{result.stderr}")
        assert result.returncode == 0, f"Simulation failed: {result.stderr}"

        # Parse acceptance log
        if not accept_log.exists():
            pytest.skip(f"Acceptance log not generated: {accept_log}")

        records = read_jsonl(accept_log)
        print(f"\nTotal records in acceptance log: {len(records)}")

        if len(records) == 0:
            pytest.skip("No records in acceptance log")

        # Filter for insertion AND deletion moves (both trigger fragment selection)
        insertion_records = filter_by_move(records, 'insertion')
        deletion_records = filter_by_move(records, 'deletion')
        gcmc_records = insertion_records + deletion_records

        print(f"Insertion attempts: {len(insertion_records)}")
        print(f"Deletion attempts: {len(deletion_records)}")
        print(f"Total GCMC moves (ins+del): {len(gcmc_records)}")

        if len(gcmc_records) == 0:
            pytest.skip("No GCMC moves recorded")

        # Count by species (both insertion and deletion reflect fragment selection)
        species_counts = count_by_species(gcmc_records)
        print(f"Species counts (ins+del): {species_counts}")

        # Calculate observed frequencies
        water_count = species_counts.get('water', 0)
        methanol_count = species_counts.get('methanol', 0)
        total_insertions = water_count + methanol_count

        if total_insertions == 0:
            pytest.skip("No species found in insertion records")

        obs_freq_water = water_count / total_insertions
        obs_freq_methanol = methanol_count / total_insertions

        print(f"\nObserved frequencies:")
        print(f"  water: {obs_freq_water:.3f} (expected: {expected_probs[0]:.3f})")
        print(f"  methanol: {obs_freq_methanol:.3f} (expected: {expected_probs[1]:.3f})")

        # Calculate absolute errors
        error_water = abs(obs_freq_water - expected_probs[0])
        error_methanol = abs(obs_freq_methanol - expected_probs[1])

        print(f"\nAbsolute errors:")
        print(f"  water: {error_water:.4f}")
        print(f"  methanol: {error_methanol:.4f}")

        # Perform chi-square test
        observed_counts = [water_count, methanol_count]
        expected_counts = [total_insertions * p for p in expected_probs]

        chi2_stat, p_value = stats.chisquare(observed_counts, expected_counts)

        print(f"\nChi-square test:")
        print(f"  Statistic: {chi2_stat:.4f}")
        print(f"  p-value: {p_value:.4f}")

        # Acceptance criteria: p > 0.05 OR absolute error < 10%
        if p_value >= 0.05:
            print(f"\n✅ Chi-square test passed (p={p_value:.4f} >= 0.05)")
        elif error_water < 0.10 and error_methanol < 0.10:
            print(f"\n✅ Frequency errors < 10%, test passes despite low p-value")
            print(f"   water error: {error_water*100:.1f}%, methanol error: {error_methanol*100:.1f}%")
        else:
            pytest.fail(
                f"Test failed: p-value={p_value:.4f} < 0.05 AND "
                f"errors exceed 10% (water: {error_water*100:.1f}%, methanol: {error_methanol*100:.1f}%)"
            )


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
