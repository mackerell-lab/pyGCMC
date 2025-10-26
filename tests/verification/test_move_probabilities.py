"""
Test mc_move_prob parsing and CDF generation.

Verifies P0.2: mc_move_prob生效验证
- Global mc_move_prob is parsed correctly
- Single-element vector is broadcasted to all fragments
- CDF is correctly normalized
"""

import pytest
import subprocess
import re
from pathlib import Path

# Path to gcmc_cpu executable
GCMC_CPU_PATH = Path(__file__).parent.parent.parent / "build" / "bin" / "gcmc_cpu"


class TestMoveProbabilities:
    """Test move probability parsing and CDF generation"""

    @staticmethod
    def create_minimal_system(tmpdir, mc_move_prob_line=None, fragment_specific=False):
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
fragmuex: -5.60

box_size: 20.0 20.0 20.0
cutoff: 10.0
temperature: 300
mcsteps: 10
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

    @staticmethod
    def extract_cdf_from_stdout(stdout):
        """Extract CDF values from stdout"""
        # Look for lines like: "  water cdf: [0.1, 0.3, 0.6, 1]"
        pattern = r'(\w+)\s+cdf:\s*\[([0-9.]+),\s*([0-9.]+),\s*([0-9.]+),\s*([0-9.]+)\]'
        matches = re.findall(pattern, stdout)

        cdfs = {}
        for match in matches:
            frag_name = match[0]
            cdf = [float(match[1]), float(match[2]), float(match[3]), float(match[4])]
            cdfs[frag_name] = cdf

        return cdfs

    def test_global_mc_move_prob_cdf(self, tmp_path):
        """
        P0.2验收测试：全局mc_move_prob设置，验证CDF正确归一化

        设置 mc_move_prob: 1 2 3 4
        期望 CDF: [0.1, 0.3, 0.6, 1.0]
        """
        inp_file = self.create_minimal_system(tmp_path, mc_move_prob_line="1 2 3 4")

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--seed", "42"],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(tmp_path)
        )

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"

        # Extract CDF from stdout
        cdfs = self.extract_cdf_from_stdout(result.stdout)

        # Should have one fragment (water)
        assert "water" in cdfs, f"Missing 'water' CDF in stdout. Got: {list(cdfs.keys())}"

        water_cdf = cdfs["water"]

        # Expected CDF: [1/10, 3/10, 6/10, 10/10] = [0.1, 0.3, 0.6, 1.0]
        expected_cdf = [0.1, 0.3, 0.6, 1.0]

        print(f"\n=== CDF Verification ===")
        print(f"Input mc_move_prob: 1 2 3 4")
        print(f"Expected CDF: {expected_cdf}")
        print(f"Actual CDF:   {water_cdf}")

        # Check each CDF value with tolerance
        for i, (expected, actual) in enumerate(zip(expected_cdf, water_cdf)):
            assert abs(expected - actual) < 0.001, \
                f"CDF[{i}]: expected={expected:.3f}, actual={actual:.3f}"
            print(f"✅ CDF[{i}]: {actual:.3f} (expected {expected:.3f})")

        print(f"\n✅ Global mc_move_prob CDF verification passed!")

    def test_mc_move_prob_broadcasting(self, tmp_path):
        """
        验证单元素向量广播到所有fragments

        即使只有一个fragment，广播逻辑也应该工作
        """
        inp_file = self.create_minimal_system(tmp_path, mc_move_prob_line="2 3 4 1")

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--seed", "42"],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(tmp_path)
        )

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"

        # Check that broadcast message appears in stdout
        assert "Broadcasted mc_move_prob" in result.stdout or "Parsed mc_move_prob" in result.stdout, \
            "Expected broadcast or parse message in stdout"

        # Extract CDF
        cdfs = self.extract_cdf_from_stdout(result.stdout)
        assert "water" in cdfs

        water_cdf = cdfs["water"]

        # Expected CDF: [2/10, 5/10, 9/10, 10/10] = [0.2, 0.5, 0.9, 1.0]
        expected_cdf = [0.2, 0.5, 0.9, 1.0]

        print(f"\n=== Broadcasting Verification ===")
        print(f"Input mc_move_prob: 2 3 4 1")
        print(f"Expected CDF: {expected_cdf}")
        print(f"Actual CDF:   {water_cdf}")

        for i, (expected, actual) in enumerate(zip(expected_cdf, water_cdf)):
            assert abs(expected - actual) < 0.001, \
                f"CDF[{i}]: expected={expected:.3f}, actual={actual:.3f}"

        print(f"✅ Broadcasting verification passed!")

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
mcsteps: 10
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

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--seed", "42"],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(tmp_path)
        )

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"

        # Extract CDF
        cdfs = self.extract_cdf_from_stdout(result.stdout)
        assert "water" in cdfs

        water_cdf = cdfs["water"]

        # Expected CDF from per-fragment [5,4,3,2]: [5/14, 9/14, 12/14, 14/14]
        expected_cdf = [5/14, 9/14, 12/14, 1.0]

        print(f"\n=== Priority Verification ===")
        print(f"Global mc_move_prob: 1 1 1 1 (should be ignored)")
        print(f"Per-fragment: 5 4 3 2")
        print(f"Expected CDF: {expected_cdf}")
        print(f"Actual CDF:   {water_cdf}")

        for i, (expected, actual) in enumerate(zip(expected_cdf, water_cdf)):
            assert abs(expected - actual) < 0.001, \
                f"CDF[{i}]: expected={expected:.3f}, actual={actual:.3f}"

        print(f"✅ Per-fragment priority verification passed!")

    def test_invalid_mc_move_prob_fallback(self, tmp_path):
        """
        验证错误的mc_move_prob参数数量时的回退行为

        提供少于4个参数，应该警告并使用默认值
        """
        inp_file = self.create_minimal_system(tmp_path, mc_move_prob_line="1 2")  # Only 2 values

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--seed", "42"],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(tmp_path)
        )

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"

        # Should see warning about requiring 4 values
        assert "WARNING" in result.stdout or "WARNING" in result.stderr, \
            "Expected warning about invalid mc_move_prob"

        # Extract CDF - should use defaults
        cdfs = self.extract_cdf_from_stdout(result.stdout)
        assert "water" in cdfs

        water_cdf = cdfs["water"]

        # Expected default CDF: [0.25, 0.5, 0.75, 1.0]
        expected_cdf = [0.25, 0.5, 0.75, 1.0]

        print(f"\n=== Fallback Verification ===")
        print(f"Invalid mc_move_prob: 1 2 (only 2 values)")
        print(f"Expected default CDF: {expected_cdf}")
        print(f"Actual CDF:   {water_cdf}")

        for i, (expected, actual) in enumerate(zip(expected_cdf, water_cdf)):
            assert abs(expected - actual) < 0.001, \
                f"CDF[{i}]: expected={expected:.3f}, actual={actual:.3f}"

        print(f"✅ Fallback to defaults verification passed!")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
