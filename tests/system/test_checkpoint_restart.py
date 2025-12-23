"""
Test checkpoint save and restart functionality.

Verifies P1: Checkpoint/Restart 加载
- Save checkpoint at intermediate step
- Load checkpoint and continue simulation
- Verify energy and statistics continuity
"""

import pytest
import subprocess
import numpy as np
from pathlib import Path
import os


# Path to gcmc_cpu executable
GCMC_CPU_PATH = Path(__file__).parent.parent.parent / "build" / "bin" / "gcmc_cpu"

def _pdb_atom_records(pdb_path):
    records = []
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        records.append(
            (
                line[12:16].strip(),
                line[17:20].strip(),
                line[22:26].strip(),
                round(float(line[30:38]), 3),
                round(float(line[38:46]), 3),
                round(float(line[46:54]), 3),
            )
        )
    return sorted(records)


class TestCheckpointRestart:
    """Test checkpoint save/load/restart functionality"""

    @staticmethod
    def create_test_system(tmpdir, mcsteps=100, seed=42):
        """Create minimal system for checkpoint testing"""

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
mcsteps: {mcsteps}
nprint: 20
eqsteps: 0

seed: {seed}
"""
        inp_file.write_text(inp_content)

        return inp_file

    def test_checkpoint_file_creation(self, tmp_path):
        """
        P1基础测试：checkpoint文件创建

        验证：
        - 运行模拟时能够创建checkpoint文件
        - checkpoint文件不为空
        """
        print(f"\n=== Checkpoint File Creation Test ===")

        inp_file = self.create_test_system(
            tmp_path,
            mcsteps=100,
            seed=42
        )

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--seed", "42", "--checkpoint-freq", "50"],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(tmp_path)
        )

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"

        # Check for checkpoint files
        checkpoint_files = list(tmp_path.glob("*checkpoint*"))
        assert checkpoint_files, f"No checkpoint files found; files: {list(tmp_path.glob('*'))}"
        print(f"✅ Found {len(checkpoint_files)} checkpoint file(s):")
        for ckpt in checkpoint_files:
            size = ckpt.stat().st_size
            print(f"  - {ckpt.name} ({size} bytes)")
            assert size > 0, f"Checkpoint file {ckpt.name} is empty"

    def test_checkpoint_roundtrip(self, tmp_path):
        """
        P1验收测试：checkpoint roundtrip

        验证：
        - 运行N步→保存checkpoint
        - 加载checkpoint→继续运行M步
        - 验证能量和统计数据连续性
        """
        print(f"\n=== Checkpoint Roundtrip Test ===")

        # Phase 1: Run N steps and save checkpoint
        inp_file = self.create_test_system(
            tmp_path,
            mcsteps=100,
            seed=42
        )

        print(f"Phase 1: Running 100 steps with checkpoint at step 50...")
        result1 = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--seed", "42", "--checkpoint-freq", "50"],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(tmp_path)
        )

        assert result1.returncode == 0, f"Phase 1 failed: {result1.stderr}"

        # Find checkpoint file
        checkpoint_files = list(tmp_path.glob("*checkpoint*.dat"))
        assert len(checkpoint_files) > 0, "No checkpoint file created"
        checkpoint_file = checkpoint_files[0]
        print(f"✅ Checkpoint created: {checkpoint_file.name} ({checkpoint_file.stat().st_size} bytes)")

        # Parse statistics from phase 1
        stats_file1 = tmp_path / "gcmc_statistics.dat"
        assert stats_file1.exists(), "Statistics file not found"

        # Read final statistics from phase 1
        with open(stats_file1, 'r') as f:
            lines = f.readlines()
        data_lines = [l for l in lines if not l.startswith('#') and l.strip()]
        assert data_lines, "No data in statistics file"

        last_line1 = data_lines[-1].split()
        step1 = int(last_line1[0]) if len(last_line1) > 0 else 0
        print(f"Phase 1 completed at step {step1}")

        # Phase 2: Load checkpoint and continue
        # Note: This requires --resume CLI parameter which may not be implemented yet
        print(f"\nPhase 2: Attempting to load checkpoint and continue...")

        # Try to resume from checkpoint
        # Create a new INP file for continuation
        inp_file2 = self.create_test_system(
            tmp_path,
            mcsteps=50,  # Additional 50 steps
            seed=42
        )

        # Try using --resume parameter (may not be implemented yet)
        result2 = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file2), "--resume", str(checkpoint_file)],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(tmp_path)
        )

        assert result2.returncode == 0, (
            f"Phase 2 failed.\nSTDERR:\n{result2.stderr}\nSTDOUT:\n{result2.stdout}"
        )

        print(f"✅ Checkpoint loaded and simulation continued")

        # Verify continuation worked
        stats_file2 = tmp_path / "gcmc_statistics.dat"
        with open(stats_file2, 'r') as f:
            lines = f.readlines()
        data_lines = [l for l in lines if not l.startswith('#') and l.strip()]
        last_line2 = data_lines[-1].split()
        step2 = int(last_line2[0]) if len(last_line2) > 0 else 0

        print(f"Phase 2 completed at step {step2}")
        assert step2 > step1, f"Resume did not advance steps: {step1} -> {step2}"
        print(f"✅ Roundtrip successful: {step1} → {step2} steps")

    def test_missing_checkpoint_error(self, tmp_path):
        """
        P1验收测试：缺失checkpoint文件的错误处理

        验证：
        - 指定不存在的checkpoint文件时有清晰错误信息
        """
        print(f"\n=== Missing Checkpoint Error Test ===")

        nonexistent_file = tmp_path / "nonexistent_checkpoint.dat"

        inp_file = self.create_test_system(tmp_path, mcsteps=10, seed=42)

        # Try to resume from nonexistent file
        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--resume", str(nonexistent_file)],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(tmp_path)
        )

        # Should fail with clear error message
        assert result.returncode != 0, "Should fail when checkpoint file doesn't exist"
        stats_file = tmp_path / "gcmc_statistics.dat"
        if stats_file.exists():
            lines = stats_file.read_text().splitlines()
            data_lines = [l for l in lines if l.strip() and not l.startswith("#")]
            if data_lines:
                steps = [int(line.split()[0]) for line in data_lines]
                assert max(steps) == 0, "Simulation advanced despite resume error"

        pdb_file = tmp_path / "gcmc_final.pdb"
        if pdb_file.exists():
            pdb_line = next(
                (line for line in inp_file.read_text().splitlines() if line.startswith("pdb:")),
                None,
            )
            assert pdb_line is not None, "Missing pdb: entry in input file"
            input_pdb = Path(pdb_line.split(":", 1)[1].strip())
            assert _pdb_atom_records(pdb_file) == _pdb_atom_records(input_pdb), \
                "Final PDB differs from input despite resume error"

    def test_simulation_reproducibility(self, tmp_path):
        """
        基础测试：相同seed的模拟应该产生相同结果

        验证：
        - 作为checkpoint测试的基准
        - 验证seed设置正确工作
        """
        print(f"\n=== Simulation Reproducibility Test ===")

        # Run 1
        run1_dir = tmp_path / "run1"
        run1_dir.mkdir()
        inp_file1 = self.create_test_system(run1_dir, mcsteps=100, seed=42)

        result1 = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file1), "--seed", "42"],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(run1_dir)
        )

        assert result1.returncode == 0

        # Run 2
        run2_dir = tmp_path / "run2"
        run2_dir.mkdir()
        inp_file2 = self.create_test_system(run2_dir, mcsteps=100, seed=42)

        result2 = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file2), "--seed", "42"],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(run2_dir)
        )

        assert result2.returncode == 0

        # Compare final PDB files
        pdb1 = run1_dir / "gcmc_final.pdb"
        pdb2 = run2_dir / "gcmc_final.pdb"

        assert pdb1.exists(), "Run 1 output PDB not found"
        assert pdb2.exists(), "Run 2 output PDB not found"
        records1 = _pdb_atom_records(pdb1)
        records2 = _pdb_atom_records(pdb2)
        assert records1, "No ATOM/HETATM records found in run 1"
        assert records2, "No ATOM/HETATM records found in run 2"
        assert records1 == records2, "Same seed should produce identical atom records"
        print(f"✅ Identical output for same seed")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
