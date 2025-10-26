"""
Reproducibility tests for gcmc_cpu with fixed seed.

Verifies that same seed produces identical results across multiple runs.
Part of P0.1: Unified RNG reproducibility.
"""

import pytest
import subprocess
import numpy as np
from pathlib import Path

# Path to gcmc_cpu executable
GCMC_CPU_PATH = Path(__file__).parent.parent.parent / "build" / "bin" / "gcmc_cpu"


class TestSeedReproducibility:
    """Test RNG reproducibility with fixed seed"""

    @staticmethod
    def create_test_system(tmpdir, seed, mcsteps=1000):
        """Create minimal water system with fixed seed"""

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

        # INP file with fixed seed
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
nprint: 200
eqsteps: 10

mc_move_prob: 0.25 0.25 0.25 0.25
seed: {seed}

op_top: {tmpdir}/output.top
op_pdb: {tmpdir}/output.pdb
"""
        inp_file.write_text(inp_content)

        return inp_file

    @staticmethod
    def parse_statistics_dat(dat_file):
        """Parse statistics.dat file and return final statistics"""
        content = dat_file.read_text()
        lines = content.strip().split('\n')
        data_lines = [l for l in lines if not l.startswith('#') and l.strip()]

        if not data_lines:
            return None

        # Parse last line
        last_line = data_lines[-1].split()

        return {
            'step': int(last_line[0]),
            'energy': float(last_line[1]),
            'n_total': int(last_line[2]),
            'accept_rate': float(last_line[3]),
            'ins_att': int(last_line[4]),
            'ins_acc': int(last_line[5]),
            'del_att': int(last_line[6]),
            'del_acc': int(last_line[7]),
            'trn_att': int(last_line[8]),
            'trn_acc': int(last_line[9]),
            'rot_att': int(last_line[10]),
            'rot_acc': int(last_line[11]),
        }

    def test_seed_reproducibility_exact(self, tmp_path):
        """
        P0.1 验收测试：相同seed两次运行，所有统计指标完全相同。

        Verifies:
        - Insert/Delete/Translate/Rotate attempts identical
        - Insert/Delete/Translate/Rotate accepted identical
        - Final molecule count identical
        - Final energy identical
        """
        seed = 42
        mcsteps = 1000

        # Run 1
        run1_dir = tmp_path / "run1"
        run1_dir.mkdir()
        inp_file1 = self.create_test_system(run1_dir, seed, mcsteps)

        result1 = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file1), "--seed", str(seed)],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(run1_dir)
        )

        assert result1.returncode == 0, f"Run 1 failed: {result1.stderr}"

        # Get statistics from Run 1
        dat_files1 = list(run1_dir.glob("*statistics.dat"))
        assert len(dat_files1) > 0, "Run 1: No statistics DAT file"
        stats1 = self.parse_statistics_dat(dat_files1[0])
        assert stats1 is not None, "Run 1: Failed to parse statistics"

        # Run 2 (same seed, clean directory)
        run2_dir = tmp_path / "run2"
        run2_dir.mkdir()
        inp_file2 = self.create_test_system(run2_dir, seed, mcsteps)

        result2 = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file2), "--seed", str(seed)],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(run2_dir)
        )

        assert result2.returncode == 0, f"Run 2 failed: {result2.stderr}"

        # Get statistics from Run 2
        dat_files2 = list(run2_dir.glob("*statistics.dat"))
        assert len(dat_files2) > 0, "Run 2: No statistics DAT file"
        stats2 = self.parse_statistics_dat(dat_files2[0])
        assert stats2 is not None, "Run 2: Failed to parse statistics"

        # Strict equality checks - all values must be EXACTLY identical
        print(f"\n=== Reproducibility Check (seed={seed}) ===")
        print(f"Run 1: {stats1}")
        print(f"Run 2: {stats2}")

        # Check all fields
        for key in ['ins_att', 'ins_acc', 'del_att', 'del_acc',
                    'trn_att', 'trn_acc', 'rot_att', 'rot_acc',
                    'n_total']:
            assert stats1[key] == stats2[key], \
                f"❌ {key}: Run1={stats1[key]} != Run2={stats2[key]}"
            print(f"✅ {key}: {stats1[key]} (exact match)")

        # Energy should match (float comparison with tiny tolerance for numerical precision)
        energy_diff = abs(stats1['energy'] - stats2['energy'])
        assert energy_diff < 1e-10, \
            f"❌ energy: Run1={stats1['energy']:.10f} != Run2={stats2['energy']:.10f} (diff={energy_diff})"
        print(f"✅ energy: {stats1['energy']:.6f} (diff={energy_diff:.2e})")

        print(f"\n✅ Perfect reproducibility achieved with seed={seed}!")

    def test_different_seeds_produce_different_results(self, tmp_path):
        """
        Sanity check: different seeds should produce different results.

        This verifies that seed is actually being used.
        """
        mcsteps = 500

        # Run with seed=42
        run1_dir = tmp_path / "run_seed42"
        run1_dir.mkdir()
        inp_file1 = self.create_test_system(run1_dir, seed=42, mcsteps=mcsteps)

        result1 = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file1), "--seed", "42"],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(run1_dir)
        )

        assert result1.returncode == 0
        dat_files1 = list(run1_dir.glob("*statistics.dat"))
        stats1 = self.parse_statistics_dat(dat_files1[0])

        # Run with seed=12345
        run2_dir = tmp_path / "run_seed12345"
        run2_dir.mkdir()
        inp_file2 = self.create_test_system(run2_dir, seed=12345, mcsteps=mcsteps)

        result2 = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file2), "--seed", "12345"],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(run2_dir)
        )

        assert result2.returncode == 0
        dat_files2 = list(run2_dir.glob("*statistics.dat"))
        stats2 = self.parse_statistics_dat(dat_files2[0])

        # At least one field should be different (very high probability)
        fields_different = [
            stats1[key] != stats2[key]
            for key in ['ins_att', 'del_att', 'trn_att', 'rot_att', 'n_total']
        ]

        assert any(fields_different), \
            "Different seeds produced identical results - RNG not working!"

        print(f"\n✅ Different seeds produce different results (as expected)")
        print(f"Seed 42:    InsAtt={stats1['ins_att']}, N={stats1['n_total']}")
        print(f"Seed 12345: InsAtt={stats2['ins_att']}, N={stats2['n_total']}")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
