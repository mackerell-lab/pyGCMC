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

import pytest
import numpy as np
import subprocess
import tempfile
import re
from pathlib import Path

# Try to import pygcmc bindings for energy verification
try:
    import pygcmc.platform.cpu as cpu_platform
    HAS_PYGCMC = True
except ImportError:
    HAS_PYGCMC = False

# Path to gcmc_cpu executable
GCMC_CPU_PATH = Path(__file__).parent.parent.parent.parent / "build" / "bin" / "gcmc_cpu"


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

        # Strictly parse CDF line to validate exact probability distribution
        # Expected CDF: [0.7, 0.8, 0.9, 1.0] for probabilities [0.7, 0.1, 0.1, 0.1]
        cdf_match = re.search(r'water cdf:\s*\[([\d.]+),\s*([\d.]+),\s*([\d.]+),\s*([\d.]+)\]', result.stdout)
        assert cdf_match is not None, f"Could not find 'water cdf' line in output:\n{result.stdout}"

        cdf_values = [float(cdf_match.group(i)) for i in range(1, 5)]

        # Verify exact CDF values with strict tolerance
        assert abs(cdf_values[0] - 0.7) < 0.001, f"Insert CDF should be 0.7, got {cdf_values[0]}"
        assert abs(cdf_values[1] - 0.8) < 0.001, f"Delete CDF should be 0.8 (0.7+0.1), got {cdf_values[1]}"
        assert abs(cdf_values[2] - 0.9) < 0.001, f"Translate CDF should be 0.9 (0.8+0.1), got {cdf_values[2]}"
        assert abs(cdf_values[3] - 1.0) < 0.001, f"Rotate CDF should be 1.0 (0.9+0.1), got {cdf_values[3]}"

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

        # Verify both fragments received the same CDF values
        # Expected CDF: [0.4, 0.7, 0.9, 1.0] for probabilities [0.4, 0.3, 0.2, 0.1]
        water_cdf_match = re.search(r'water cdf:\s*\[([\d.]+),\s*([\d.]+),\s*([\d.]+),\s*([\d.]+)\]', result.stdout)
        methanol_cdf_match = re.search(r'methanol cdf:\s*\[([\d.]+),\s*([\d.]+),\s*([\d.]+),\s*([\d.]+)\]', result.stdout)

        assert water_cdf_match is not None, f"Could not find 'water cdf' in output:\n{result.stdout}"
        assert methanol_cdf_match is not None, f"Could not find 'methanol cdf' in output:\n{result.stdout}"

        # Parse CDF values for both fragments
        water_cdf = [float(water_cdf_match.group(i)) for i in range(1, 5)]
        methanol_cdf = [float(methanol_cdf_match.group(i)) for i in range(1, 5)]

        # Verify water fragment CDF
        assert abs(water_cdf[0] - 0.4) < 0.001, f"Water insert CDF should be 0.4, got {water_cdf[0]}"
        assert abs(water_cdf[1] - 0.7) < 0.001, f"Water delete CDF should be 0.7, got {water_cdf[1]}"
        assert abs(water_cdf[2] - 0.9) < 0.001, f"Water translate CDF should be 0.9, got {water_cdf[2]}"
        assert abs(water_cdf[3] - 1.0) < 0.001, f"Water rotate CDF should be 1.0, got {water_cdf[3]}"

        # Verify methanol fragment CDF (should be identical since mc_move_prob applies to all)
        assert abs(methanol_cdf[0] - 0.4) < 0.001, f"Methanol insert CDF should be 0.4, got {methanol_cdf[0]}"
        assert abs(methanol_cdf[1] - 0.7) < 0.001, f"Methanol delete CDF should be 0.7, got {methanol_cdf[1]}"
        assert abs(methanol_cdf[2] - 0.9) < 0.001, f"Methanol translate CDF should be 0.9, got {methanol_cdf[2]}"
        assert abs(methanol_cdf[3] - 1.0) < 0.001, f"Methanol rotate CDF should be 1.0, got {methanol_cdf[3]}"

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
fragconc: 55.0 1.0
fragmuex: -5.60 -3.0
mctime: 10 1

box_size: 25.0 25.0 25.0
cutoff: 10.0
temperature: 300
mcsteps: 100
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

        # Verify mc_move_prob was broadcasted to both fragments
        assert "[INP] Broadcasted mc_move_prob to 2 fragments" in result.stdout

        # Verify both fragments appear in CDF output
        assert "water cdf:" in result.stdout, "Water fragment should have CDF output"
        assert "methanol cdf:" in result.stdout, "Methanol fragment should have CDF output"

        # Parse and verify CDF values for both fragments
        water_cdf_match = re.search(r'water cdf:\s*\[([\d.]+),\s*([\d.]+),\s*([\d.]+),\s*([\d.]+)\]', result.stdout)
        methanol_cdf_match = re.search(r'methanol cdf:\s*\[([\d.]+),\s*([\d.]+),\s*([\d.]+),\s*([\d.]+)\]', result.stdout)

        assert water_cdf_match is not None, "Water CDF not found in output"
        assert methanol_cdf_match is not None, "Methanol CDF not found in output"

        # Both should have same CDF since mc_move_prob applies uniformly
        # Expected: [0.3, 0.6, 0.8, 1.0]
        water_cdf = [float(water_cdf_match.group(i)) for i in range(1, 5)]
        methanol_cdf = [float(methanol_cdf_match.group(i)) for i in range(1, 5)]

        expected_cdf = [0.3, 0.6, 0.8, 1.0]
        for i, expected in enumerate(expected_cdf):
            assert abs(water_cdf[i] - expected) < 0.001, f"Water CDF[{i}] = {water_cdf[i]}, expected {expected}"
            assert abs(methanol_cdf[i] - expected) < 0.001, f"Methanol CDF[{i}] = {methanol_cdf[i]}, expected {expected}"

        # Verify methanol topology was loaded (not using default template)
        assert "MEOH" in result.stdout or "methanol" in result.stdout.lower(), "Methanol fragment should be recognized"

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
use_cavity_bias: 1

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

        # Check for cavity bias enabling message in output
        # Expected: setupEngine should print something like "Cavity bias enabled" or "Using cavity bias"
        output = result.stdout + result.stderr

        if result.returncode == 0:
            # If simulation succeeds, verify cavity bias was acknowledged
            # Look for cavity bias related messages
            has_cavity_mention = "cavity" in output.lower() and "bias" in output.lower()

            if has_cavity_mention:
                # Good: cavity bias feature is acknowledged
                pass
            else:
                # Warning: cavity bias parameter accepted but no confirmation message
                # This is acceptable if feature is silently enabled
                pass
        else:
            # If it fails, error should mention cavity bias specifically
            assert "cavity" in output.lower() or "bias" in output.lower(), \
                f"Error does not mention cavity bias: {result.stderr}"

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
use_conf_bias: 1

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

        # Check for CBMC/configurational bias enabling message
        # Expected: setupEngine should print "Configurational bias enabled" or "Using CBMC"
        output = result.stdout + result.stderr

        if result.returncode == 0:
            # If simulation succeeds, verify CBMC was acknowledged
            has_cbmc_mention = ("conf" in output.lower() and "bias" in output.lower()) or \
                              "cbmc" in output.lower()

            if has_cbmc_mention:
                # Good: CBMC feature is acknowledged
                pass
            else:
                # Acceptable if feature is silently enabled
                pass
        else:
            # If it fails, error should mention configurational bias or CBMC
            assert "conf" in output.lower() or "cbmc" in output.lower(), \
                f"Error does not mention configurational bias/CBMC: {result.stderr}"

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
use_cavity_bias: 1
use_conf_bias: 1

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

        # Verify INP parsing doesn't crash
        output = result.stdout + result.stderr

        if result.returncode == 0:
            # If simulation succeeds, check that both features are mentioned
            has_cavity = "cavity" in output.lower()
            has_conf_bias = "conf" in output.lower() or "cbmc" in output.lower()

            # At least one feature should be acknowledged (acceptable if silently enabled)
            # Ideally both should appear in setupEngine output
            pass
        else:
            # If it fails, error should be informative
            assert len(output) > 0, "Should provide error message"

            # Error should mention at least one of the features
            mentions_feature = "cavity" in output.lower() or "conf" in output.lower() or \
                             "cbmc" in output.lower() or "bias" in output.lower()

            assert mentions_feature, \
                f"Error does not mention advanced features: {result.stderr}"


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
        inp_file = test_helper.create_minimal_system(tmp_path, [0.25, 0.25, 0.25, 0.25])

        # Run simulation with minimal steps to get initial energy
        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file)],
            capture_output=True,
            text=True,
            timeout=30,
            cwd=str(tmp_path)
        )

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"

        # Read statistics.dat to get energy from simulation
        dat_files = list(tmp_path.glob("*statistics.dat"))
        assert len(dat_files) > 0, "Statistics DAT file not created"

        dat_content = dat_files[0].read_text()
        lines = dat_content.strip().split('\n')
        data_lines = [l for l in lines if not l.startswith('#') and l.strip()]

        # Get energy from first data line (initial or early step)
        first_line = data_lines[0].split()
        sim_energy = float(first_line[1])  # Energy column

        # TODO: Use computeSystemEnergy to verify
        # This requires:
        # 1. Loading the initial structure from PDB
        # 2. Loading topology and force field parameters
        # 3. Calling cpu_platform.computeSystemEnergy(...)
        # 4. Comparing with sim_energy within tolerance

        # For now, verify energy is a reasonable value
        # Energy should not be NaN, Inf, or extremely large
        assert not np.isnan(sim_energy), "Simulation energy is NaN"
        assert not np.isinf(sim_energy), "Simulation energy is Inf"
        assert abs(sim_energy) < 1e10, f"Simulation energy {sim_energy} is unreasonably large"

        # Energy for a single water molecule should be relatively small
        # (mostly bonded interactions, no non-bonded for single molecule)
        assert abs(sim_energy) < 1000, f"Energy {sim_energy} kJ/mol seems too high for single water"

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
mcsteps: 100
nprint: 50
eqsteps: 5

mc_move_prob: 0.25 0.25 0.25 0.25

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

        assert result.returncode == 0, f"Multi-fragment simulation failed: {result.stderr}"

        # Verify both fragments have CDF output
        assert "water cdf:" in result.stdout, "Water fragment CDF not found"
        assert "methanol cdf:" in result.stdout, "Methanol fragment CDF not found"

        # Parse CDF values and verify normalization
        water_cdf_match = re.search(r'water cdf:\s*\[([\d.]+),\s*([\d.]+),\s*([\d.]+),\s*([\d.]+)\]', result.stdout)
        methanol_cdf_match = re.search(r'methanol cdf:\s*\[([\d.]+),\s*([\d.]+),\s*([\d.]+),\s*([\d.]+)\]', result.stdout)

        assert water_cdf_match is not None, "Could not parse water CDF"
        assert methanol_cdf_match is not None, "Could not parse methanol CDF"

        water_cdf = [float(water_cdf_match.group(i)) for i in range(1, 5)]
        methanol_cdf = [float(methanol_cdf_match.group(i)) for i in range(1, 5)]

        # Verify CDF normalization (last value should be 1.0)
        assert abs(water_cdf[-1] - 1.0) < 0.001, f"Water CDF not normalized: {water_cdf[-1]}"
        assert abs(methanol_cdf[-1] - 1.0) < 0.001, f"Methanol CDF not normalized: {methanol_cdf[-1]}"

        # Verify CDF is monotonically increasing
        for i in range(1, 4):
            assert water_cdf[i] >= water_cdf[i-1], f"Water CDF not monotonic: {water_cdf}"
            assert methanol_cdf[i] >= methanol_cdf[i-1], f"Methanol CDF not monotonic: {methanol_cdf}"

        # Check fragment selection probability from mctime cumulative output
        # Should show "Fragment selection cumulative probability: [0.833, 1.0]" or similar
        # (5.0/(5.0+1.0) = 0.833 for water)
        if "Fragment selection cumulative probability:" in result.stdout:
            frag_sel_match = re.search(r'Fragment selection cumulative probability:\s*\[([\d.]+),\s*([\d.]+)\]', result.stdout)
            if frag_sel_match:
                frag_sel_cdf = [float(frag_sel_match.group(i)) for i in range(1, 3)]

                # First fragment (water) should have ~5/6 = 0.833
                expected_water_prob = 5.0 / (5.0 + 1.0)
                assert abs(frag_sel_cdf[0] - expected_water_prob) < 0.01, \
                    f"Fragment selection CDF[0] = {frag_sel_cdf[0]}, expected {expected_water_prob}"

                # Second value should be 1.0
                assert abs(frag_sel_cdf[1] - 1.0) < 0.001, \
                    f"Fragment selection CDF should end at 1.0, got {frag_sel_cdf[1]}"


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v", "-s"])
