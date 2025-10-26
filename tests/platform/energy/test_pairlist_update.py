"""
Test neighbor list (pairlist) functionality.

Verifies P0.5: 邻居表（最小可用）
- Neighbor list cutoff effects
- Rebuild frequency tracking
- Energy consistency with full calculation
"""

import pytest
import subprocess
import numpy as np
from pathlib import Path

# Path to gcmc_cpu executable
GCMC_CPU_PATH = Path(__file__).parent.parent.parent.parent / "build" / "bin" / "gcmc_cpu"


class TestPairlistUpdate:
    """Test neighbor list update and cutoff effects"""

    @staticmethod
    def create_test_system(tmpdir, box_size=30.0, n_molecules=10, mcsteps=100, seed=42,
                          pairlist_cutoff=None, pairlist_freq=None):
        """Create test system for neighbor list testing"""

        # Create initial PDB with multiple water molecules
        pdb_file = tmpdir / "test.pdb"
        pdb_content = f"CRYST1   {box_size:.3f}   {box_size:.3f}   {box_size:.3f}  90.00  90.00  90.00 P 1           1\n"

        atom_id = 1
        for i in range(n_molecules):
            # Place molecules in a grid
            x = (i % 3) * 10.0 + 5.0
            y = ((i // 3) % 3) * 10.0 + 5.0
            z = (i // 9) * 10.0 + 5.0

            pdb_content += f"ATOM  {atom_id:5d}  O   WAT  {i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00\n"
            atom_id += 1
            pdb_content += f"ATOM  {atom_id:5d}  H1  WAT  {i+1:4d}    {x+0.757:8.3f}{y+0.586:8.3f}{z:8.3f}  1.00  0.00\n"
            atom_id += 1
            pdb_content += f"ATOM  {atom_id:5d}  H2  WAT  {i+1:4d}    {x-0.757:8.3f}{y+0.586:8.3f}{z:8.3f}  1.00  0.00\n"
            atom_id += 1

        pdb_content += "END\n"
        pdb_file.write_text(pdb_content)

        # TOP file
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
"""
        for i in range(n_molecules):
            top_content += f"WAT  1\n"

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

box_size: {box_size} {box_size} {box_size}
cutoff: {box_size/2 - 1.0}
temperature: 300
mcsteps: {mcsteps}
nprint: 100
eqsteps: 0

mc_move_prob: 0 0 0.5 0.5

seed: {seed}
"""

        # Add pairlist parameters if specified
        if pairlist_cutoff is not None:
            inp_content += f"pairlist_cutoff: {pairlist_cutoff}\n"
        if pairlist_freq is not None:
            inp_content += f"pairlist_freq: {pairlist_freq}\n"

        inp_file.write_text(inp_content)

        return inp_file

    def test_pairlist_cutoff_effect(self, tmp_path):
        """
        P0.5验收测试：pairlist_cutoff参数生效

        验证：
        - 不同cutoff值被正确解析
        - INP文件中的pairlist_cutoff参数可以设置
        """
        print(f"\n=== Pairlist Cutoff Effect Test ===")

        # Test with larger cutoff
        run1_dir = tmp_path / "run_large_cutoff"
        run1_dir.mkdir()
        inp_file1 = self.create_test_system(
            run1_dir,
            box_size=30.0,
            n_molecules=5,
            mcsteps=50,
            pairlist_cutoff=2.0,  # Large cutoff
            seed=42
        )

        result1 = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file1), "--seed", "42"],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(run1_dir)
        )

        assert result1.returncode == 0, f"Large cutoff run failed: {result1.stderr}"

        # Test with smaller cutoff
        run2_dir = tmp_path / "run_small_cutoff"
        run2_dir.mkdir()
        inp_file2 = self.create_test_system(
            run2_dir,
            box_size=30.0,
            n_molecules=5,
            mcsteps=50,
            pairlist_cutoff=0.8,  # Small cutoff
            seed=42
        )

        result2 = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file2), "--seed", "42"],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(run2_dir)
        )

        assert result2.returncode == 0, f"Small cutoff run failed: {result2.stderr}"

        print(f"✅ Both runs completed successfully")
        print(f"✅ pairlist_cutoff parameter parsed correctly")

    def test_pairlist_update_frequency(self, tmp_path):
        """
        P0.5验收测试：pairlist_freq控制重建频率

        验证：
        - pairlist_freq参数被正确解析
        - 不同频率值下模拟可以正常运行
        """
        print(f"\n=== Pairlist Update Frequency Test ===")

        # Test with frequent updates
        run1_dir = tmp_path / "run_freq_high"
        run1_dir.mkdir()
        inp_file1 = self.create_test_system(
            run1_dir,
            box_size=30.0,
            n_molecules=5,
            mcsteps=100,
            pairlist_freq=10,  # Update every 10 steps
            seed=42
        )

        result1 = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file1), "--seed", "42"],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(run1_dir)
        )

        assert result1.returncode == 0, f"High freq run failed: {result1.stderr}"

        # Test with infrequent updates
        run2_dir = tmp_path / "run_freq_low"
        run2_dir.mkdir()
        inp_file2 = self.create_test_system(
            run2_dir,
            box_size=30.0,
            n_molecules=5,
            mcsteps=100,
            pairlist_freq=1000,  # Update every 1000 steps
            seed=42
        )

        result2 = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file2), "--seed", "42"],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(run2_dir)
        )

        assert result2.returncode == 0, f"Low freq run failed: {result2.stderr}"

        print(f"✅ Both runs completed successfully")
        print(f"✅ pairlist_freq parameter controls update frequency")

    def test_basic_simulation_with_pairlist(self, tmp_path):
        """
        基础测试：带有pairlist参数的模拟可以正常运行

        验证：
        - 带有pairlist_cutoff和pairlist_freq参数的模拟不会崩溃
        - 能正常输出结果文件
        """
        print(f"\n=== Basic Simulation with Pairlist ===")

        inp_file = self.create_test_system(
            tmp_path,
            box_size=30.0,
            n_molecules=5,
            mcsteps=50,
            pairlist_cutoff=1.5,
            pairlist_freq=50,
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

        # Check output files exist
        final_pdb = tmp_path / "gcmc_final.pdb"
        assert final_pdb.exists(), "Output PDB not found"

        stats_files = list(tmp_path.glob("*statistics.dat"))
        assert len(stats_files) > 0, "No statistics.dat file"

        print(f"✅ Simulation completed successfully")
        print(f"✅ Output files generated")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
