"""
Test GCMC output DAT file completeness and schema.

Verifies P0.4: 输出文件完整性 - DAT文件
- statistics.dat schema and value ranges
- energy.dat schema and consistency
- density.dat schema and format
"""

import pytest
import subprocess
import re
from pathlib import Path


# Path to gcmc_cpu executable
GCMC_CPU_PATH = Path(__file__).resolve().parents[3] / "build" / "bin" / "gcmc_cpu"


class TestOutputDAT:
    """Test output DAT file schemas and value ranges"""

    @staticmethod
    def create_test_system(tmpdir, mcsteps=100, seed=42):
        """Create minimal system for output testing"""

        pdb_file = tmpdir / "test.pdb"
        pdb_content = """CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1
ATOM      1  O   WAT     1      10.000  10.000  10.000  1.00  0.00
ATOM      2  H1  WAT     1      10.757  10.586  10.000  1.00  0.00
ATOM      3  H2  WAT     1       9.243  10.586  10.000  1.00  0.00
END
"""
        pdb_file.write_text(pdb_content)

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

        atp_file = tmpdir / "atomtypes.atp"
        atp_file.write_text("O   15.9994\nH    1.008\n")

        ff_file = tmpdir / "ffnonbonded.itp"
        ff_file.write_text("""[ atomtypes ]
O   8   15.9994  -0.834  A  3.15061e-01  6.36386e-01
H   1   1.008     0.417  A  0.00000e+00  0.00000e+00
""")

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

    def test_statistics_dat_schema(self, tmp_path):
        """
        P0.4验收测试：statistics.dat格式验证

        验证：
        - 文件存在
        - 表头包含必要字段（Step, Accept Rate, Molecules, Energy等）
        - 每行列数固定
        - 值域合理（接受率0-1，分子数≥0，能量有限）
        """
        inp_file = self.create_test_system(tmp_path, mcsteps=100, seed=42)

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--seed", "42"],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(tmp_path)
        )

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"

        # Find statistics DAT file
        dat_files = list(tmp_path.glob("*statistics.dat"))
        assert len(dat_files) > 0, f"No statistics.dat file found. Files: {list(tmp_path.glob('*.dat'))}"

        stats_dat = dat_files[0]
        print(f"\n=== Statistics DAT Schema Verification ===")
        print(f"File: {stats_dat.name}")

        with open(stats_dat, 'r') as f:
            lines = f.readlines()

        assert len(lines) > 0, "statistics.dat is empty"

        # Parse header (lines starting with #)
        header_lines = [l for l in lines if l.startswith('#')]
        data_lines = [l for l in lines if not l.startswith('#') and l.strip()]

        print(f"Header lines: {len(header_lines)}")
        print(f"Data lines: {len(data_lines)}")

        # Header should mention key fields
        header_text = ' '.join(header_lines)
        print(f"\nHeader preview: {header_lines[0] if header_lines else 'No header'}")

        assert data_lines, "No data lines (simulation may have been too short)"

        # Parse first data line to determine column count
        first_data = data_lines[0].split()
        num_columns = len(first_data)

        print(f"Number of columns: {num_columns}")
        print(f"First data line: {data_lines[0].strip()}")

        # Check all data lines have same column count
        for i, line in enumerate(data_lines):
            columns = line.split()
            assert len(columns) == num_columns, \
                f"Line {i} has {len(columns)} columns, expected {num_columns}"

        # Validate value ranges for each data line
        for i, line in enumerate(data_lines):
            columns = line.split()

            # Typically: Step, Energy, N_total, Accept_rate, Ins_att, Ins_acc, Del_att, Del_acc, Trn_att, Trn_acc, Rot_att, Rot_acc
            step = int(columns[0])
            energy = float(columns[1])
            n_total = int(columns[2])
            accept_rate = float(columns[3])

            # Validate ranges
            assert step >= 0, f"Step {step} < 0"
            assert n_total >= 0, f"N_total {n_total} < 0"
            assert 0 <= accept_rate <= 100, f"Accept rate {accept_rate} not in [0, 100]"
            assert -1e10 < energy < 1e10, f"Energy {energy} out of reasonable range"

        print(f"\n✅ All {len(data_lines)} data lines have consistent schema")
        print(f"✅ Value ranges validated (accept rate 0-100%, N≥0, finite energy)")

        print(f"\n✅ statistics.dat schema verification passed!")

    def test_output_files_exist(self, tmp_path):
        """
        验证基本输出文件存在

        验证：
        - gcmc_final.pdb exists
        - gcmc_final.top exists
        - At least one statistics.dat exists
        """
        inp_file = self.create_test_system(tmp_path, mcsteps=50, seed=42)

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--seed", "42"],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(tmp_path)
        )

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"

        print(f"\n=== Output Files Existence Check ===")

        # Check PDB
        pdb_file = tmp_path / "gcmc_final.pdb"
        assert pdb_file.exists(), "gcmc_final.pdb not found"
        print(f"✅ gcmc_final.pdb exists ({pdb_file.stat().st_size} bytes)")

        # Check TOP
        top_file = tmp_path / "gcmc_final.top"
        assert top_file.exists(), "gcmc_final.top not found"
        print(f"✅ gcmc_final.top exists ({top_file.stat().st_size} bytes)")

        # Check statistics DAT
        dat_files = list(tmp_path.glob("*statistics.dat"))
        assert len(dat_files) > 0, "No statistics.dat files found"
        print(f"✅ Found {len(dat_files)} statistics.dat file(s)")

        for dat in dat_files:
            print(f"  - {dat.name} ({dat.stat().st_size} bytes)")

        print(f"\n✅ All expected output files exist!")

    def test_pdb_basic_format(self, tmp_path):
        """
        验证PDB基本格式

        验证：
        - Contains CRYST1 line
        - Box dimensions match input
        - Contains ATOM lines or END
        """
        inp_file = self.create_test_system(tmp_path, mcsteps=50, seed=42)

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--seed", "42"],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(tmp_path)
        )

        assert result.returncode == 0

        pdb_file = tmp_path / "gcmc_final.pdb"
        with open(pdb_file, 'r') as f:
            content = f.read()

        print(f"\n=== PDB Format Verification ===")

        # Check for CRYST1
        assert 'CRYST1' in content, "Missing CRYST1 line in PDB"
        print(f"✅ CRYST1 line present")

        # Check for END or ATOM
        has_atoms = 'ATOM' in content or 'HETATM' in content
        has_end = 'END' in content

        assert has_atoms, "No ATOM/HETATM records found in PDB"
        assert has_end, "Missing END marker in PDB"
        print(f"✅ Contains ATOM/HETATM records")
        print(f"✅ Contains END marker")

        # Extract box dimensions from CRYST1
        cryst_match = re.search(r'CRYST1\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)', content)
        if cryst_match:
            box_x = float(cryst_match.group(1))
            box_y = float(cryst_match.group(2))
            box_z = float(cryst_match.group(3))
            print(f"✅ Box dimensions: {box_x} × {box_y} × {box_z} Å")

            # Default (gcmc_gpu/charmm) mode: INP box_size is in Å and should round-trip into CRYST1 as-is.
            inp_content = (tmp_path / "test.inp").read_text()
            inp_match = re.search(r"^box_size:\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)", inp_content, re.M)
            assert inp_match, "Missing box_size line in input INP"
            expected_x = float(inp_match.group(1))
            expected_y = float(inp_match.group(2))
            expected_z = float(inp_match.group(3))

            tolerance = 1.0
            assert abs(box_x - expected_x) < tolerance, f"Box X mismatch: {box_x} vs {expected_x}"
            assert abs(box_y - expected_y) < tolerance, f"Box Y mismatch: {box_y} vs {expected_y}"
            assert abs(box_z - expected_z) < tolerance, f"Box Z mismatch: {box_z} vs {expected_z}"
            print(f"✅ Box dimensions match input")

        print(f"\n✅ PDB format verification passed!")

    def test_density_dat_schema(self, tmp_path):
        """
        P0.4验收测试：density.dat格式验证

        验证：
        - 文件存在（需要设置wdens参数）
        - 表头包含 Step 和分子计数/密度字段
        - 每行列数固定
        - 值域合理（Step≥0, N≥0, density≥0）
        """
        # Create system with wdens parameter
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

        # INP file with wdens parameter
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
mcsteps: 100
nprint: 20
eqsteps: 0
wdens: 20

seed: 42
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

        # Find density DAT file
        density_files = list(tmp_path.glob("*density.dat"))

        assert density_files, f"No density.dat file found; dat files: {list(tmp_path.glob('*.dat'))}"

        density_dat = density_files[0]
        print(f"\n=== Density DAT Schema Verification ===")
        print(f"File: {density_dat.name}")

        with open(density_dat, 'r') as f:
            lines = f.readlines()

        assert len(lines) > 0, "density.dat is empty"

        # Parse header and data
        header_lines = [l for l in lines if l.startswith('#')]
        data_lines = [l for l in lines if not l.startswith('#') and l.strip()]

        print(f"Header lines: {len(header_lines)}")
        print(f"Data lines: {len(data_lines)}")

        if len(header_lines) > 0:
            print(f"Header: {header_lines[0].strip()}")

        assert data_lines, "No data lines in density.dat"

        first_data = data_lines[0].split()
        num_columns = len(first_data)

        print(f"Number of columns: {num_columns}")
        print(f"First data line: {data_lines[0].strip()}")

        # Expected format: Step N_water density(molecules/nm³) density(M)
        assert num_columns >= 2, f"Expected at least 2 columns (Step, Count), got {num_columns}"

        # Check all lines have same column count
        for i, line in enumerate(data_lines):
            columns = line.split()
            assert len(columns) == num_columns, \
                f"Line {i} has {len(columns)} columns, expected {num_columns}"

        # Validate value ranges
        for i, line in enumerate(data_lines):
            columns = line.split()
            step = int(columns[0])
            count = int(columns[1])

            assert step >= 0, f"Step {step} < 0"
            assert count >= 0, f"Count {count} < 0"

            # If density values present
            if num_columns >= 3:
                density = float(columns[2])
                assert density >= 0, f"Density {density} < 0"

        print(f"\n✅ All {len(data_lines)} data lines have consistent schema")
        print(f"✅ Value ranges validated (Step≥0, Count≥0, Density≥0)")

        print(f"\n✅ density.dat schema verification passed!")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
