"""
Functional validation tests for gcmc_cpu
测试 gcmc_cpu 的基本功能是否正确运行

Tests:
1. mc_move_prob parameter parsing and functionality
2. Statistics DAT file output
3. Basic GCMC simulation runs
4. INP parameter handling
"""

import pytest
import numpy as np
import subprocess
import tempfile
import re
from pathlib import Path

# Path to gcmc_cpu executable
GCMC_CPU_PATH = Path(__file__).parent.parent.parent / "build" / "bin" / "gcmc_cpu"


class TestMcMoveProb:
    """Test mc_move_prob parameter functionality"""

    @staticmethod
    def create_minimal_system(tmpdir, mc_move_prob_values=None):
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
mcsteps: 100
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
        inp_file = self.create_minimal_system(tmp_path, [0.25, 0.25, 0.25, 0.25])

        # Run gcmc_cpu in tmp_path directory
        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file)],
            capture_output=True,
            text=True,
            timeout=30,
            cwd=str(tmp_path)  # Run in tmp_path to ensure files are created there
        )

        # Check that it ran successfully
        assert result.returncode == 0, f"gcmc_cpu failed: {result.stderr}"

        # Check that mc_move_prob was parsed
        assert "[INP] Parsed mc_move_prob:" in result.stdout
        assert "0.25 (ins)" in result.stdout
        assert "0.25 (del)" in result.stdout
        assert "0.25 (trn)" in result.stdout
        assert "0.25 (rot)" in result.stdout

        # Check that probabilities were applied
        assert "Move probability cdf" in result.stdout or "move probabilities" in result.stdout
        # Should show the CDF values

    def test_mc_move_prob_biased_insertion(self, tmp_path):
        """Test mc_move_prob with bias toward insertion (0.7, 0.1, 0.1, 0.1)"""
        inp_file = self.create_minimal_system(tmp_path, [0.7, 0.1, 0.1, 0.1])

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file)],
            capture_output=True,
            text=True,
            timeout=30,
            cwd=str(tmp_path)
        )

        assert result.returncode == 0, f"gcmc_cpu failed: {result.stderr}"

        # Check that biased probabilities were parsed
        assert "[INP] Parsed mc_move_prob:" in result.stdout
        assert "0.7 (ins)" in result.stdout
        assert "0.1 (del)" in result.stdout

        # Check that move probabilities reflect the bias
        # Should show something like Ins=0.7, Del=0.1, Trn=0.1, Rot=0.1
        prob_match = re.search(r'move probabilities:.*Ins=([\d.]+)', result.stdout)
        if prob_match:
            ins_prob = float(prob_match.group(1))
            # Should be close to 0.7
            assert abs(ins_prob - 0.7) < 0.01, f"Insert probability {ins_prob} not close to 0.7"

    def test_mc_move_prob_before_fragname(self, tmp_path):
        """Test that mc_move_prob works even when placed before fragname (time-order independent)"""

        # Create files
        pdb_file, top_file, atp_file, ff_file = self._create_basic_files(tmp_path)

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
mcsteps: 50
nprint: 20
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

        assert result.returncode == 0, f"gcmc_cpu failed: {result.stderr}"

        # Check that broadcasting happened
        assert "[INP] Broadcasted mc_move_prob to 2 fragments" in result.stdout

        # Both fragments should have probability cdf output
        # Look for CDF output showing both fragments
        assert "cdf:" in result.stdout, "Should show probability CDF for fragments"

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

        # Should have run expected number of steps
        assert "Step" in result.stdout
        assert "mcsteps: 100" in result.stdout or "100" in result.stdout

        # Should report statistics
        assert "Acceptance" in result.stdout or "acceptance" in result.stdout
        assert "Energy" in result.stdout

        # Output files should exist (default prefix is "gcmc")
        top_files = list(tmp_path.glob("*.top"))
        pdb_files = list(tmp_path.glob("*.pdb"))
        assert len(top_files) > 0, f"Output TOP file not created. Files: {list(tmp_path.iterdir())}"
        assert len(pdb_files) > 0, f"Output PDB file not created. Files: {list(tmp_path.iterdir())}"

    def test_multiple_fragments(self, tmp_path):
        """Test GCMC with multiple fragment types"""

        # Create basic files
        pdb_file = tmp_path / "test.pdb"
        pdb_file.write_text("""CRYST1   25.000   25.000   25.000  90.00  90.00  90.00 P 1           1
ATOM      1  O   WAT     1      12.000  12.000  12.000  1.00  0.00
ATOM      2  H1  WAT     1      12.757  12.586  12.000  1.00  0.00
ATOM      3  H2  WAT     1      11.243  12.586  12.000  1.00  0.00
END
""")

        top_file = tmp_path / "test.top"
        top_file.write_text("""[ defaults ]
1 2 yes 0.5 0.8333
[ atomtypes ]
O   8   15.9994  -0.834  A  3.15061e-01  6.36386e-01
H   1   1.008     0.417  A  0.00000e+00  0.00000e+00
C   6   12.011    0.000  A  3.39967e-01  4.57730e-01
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

        atp_file = tmp_path / "atomtypes.atp"
        atp_file.write_text("O   15.9994\nH    1.008\nC   12.011\n")

        ff_file = tmp_path / "ffnonbonded.itp"
        ff_file.write_text("""[ atomtypes ]
O   8   15.9994  -0.834  A  3.15061e-01  6.36386e-01
H   1   1.008     0.417  A  0.00000e+00  0.00000e+00
C   6   12.011    0.000  A  3.39967e-01  4.57730e-01
""")

        # INP with two fragments
        inp_file = tmp_path / "test.inp"
        inp_content = f"""par:{ff_file}
atomtypes:{atp_file}
top:{top_file}
pdb:{pdb_file}
protitp:{top_file}

fragname: water methanol
fragconc: 55.0 1.0
fragmuex: -5.60 -3.0
mctime: 10 1

box_size: 25.0 25.0 25.0
cutoff: 10.0
temperature: 300
mcsteps: 50
nprint: 25
eqsteps: 5

mc_move_prob: 0.3 0.3 0.2 0.2

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

        # Should complete
        assert result.returncode == 0, f"Multi-fragment simulation failed: {result.stderr}"

        # Should show both fragments in output
        assert "water" in result.stdout.lower()
        # Note: methanol might not appear if not inserted, but simulation should still run

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

        # Should report energy
        energy_match = re.search(r'Energy:\s*([-\d.]+)', result.stdout)
        assert energy_match, "No energy value found in output"

        energy = float(energy_match.group(1))
        # Energy can be any value, just check it's numeric
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
fragmuex: -5.60

mc_move_prob: 0.5 0.5

box_size: 20.0 20.0 20.0
cutoff: 10.0
temperature: 300
mcsteps: 10
nprint: 10
eqsteps: 0

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

        # Should show warning about wrong number of values
        assert "[WARNING] mc_move_prob requires 4 values" in result.stdout or \
               "[WARNING] mc_move_prob requires 4 values" in result.stderr


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v", "-s"])
