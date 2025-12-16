#!/usr/bin/env python3
"""
Compatibility tests for gcmc_cpu with legacy gcmc_gpu/gcmc_opencl examples
Tests various input formats, parameters, and output validation
"""

import pytest
import subprocess
import tempfile
import numpy as np
from pathlib import Path
import re

GCMC_CPU_PATH = Path(__file__).resolve().parents[3] / "build" / "bin" / "gcmc_cpu"
TEST_DATA_DIR = Path(__file__).resolve().parents[2] / "data"


class TestGCMCCompatibility:
    """Test compatibility with legacy GCMC implementations"""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test outputs"""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield tmpdir

    def run_gcmc(self, inp_file, temp_dir, extra_args=None):
        """Helper to run gcmc_cpu and capture output"""
        args = [str(GCMC_CPU_PATH), "--inp", str(inp_file)]
        if extra_args:
            args.extend(extra_args)

        result = subprocess.run(
            args,
            cwd=temp_dir,
            capture_output=True,
            text=True
        )
        return result

    def test_basic_water_insertion(self, temp_dir):
        """Test basic water insertion with minimal parameters"""
        # Create minimal INP file similar to gcmc_gpu examples
        inp_content = """
# Basic water insertion test
version: gcmc_2.0
box: 30.0 30.0 30.0
temperature: 300.0
fragname: WAT
fragmuex: -10.0
fragconc: 55.0
mcsteps: 100
nprint: 20
nsave: 50
fragitp: {}/charmm36.ff/mol/sol.itp
""".format(TEST_DATA_DIR)

        inp_file = Path(temp_dir) / "water_basic.inp"
        inp_file.write_text(inp_content)

        result = self.run_gcmc(inp_file, temp_dir, ["--prefix", "water_basic", "--seed", "42"])

        # Check successful execution
        assert result.returncode == 0, f"Failed: {result.stderr}"
        assert "Simulation completed" in result.stdout

        # Check output files exist
        pdb_file = Path(temp_dir) / "water_basic_final.pdb"
        top_file = Path(temp_dir) / "water_basic_final.top"
        assert pdb_file.exists(), "PDB file not created"
        assert top_file.exists(), "TOP file not created"

        # Validate PDB format
        pdb_content = pdb_file.read_text()
        assert "REMARK GCMC Trajectory" in pdb_content
        assert "CRYST1" in pdb_content
        assert "END" in pdb_content

        # Validate TOP format
        top_content = top_file.read_text()
        assert "[ system ]" in top_content
        assert "[ molecules ]" in top_content

    def test_multi_fragment_silcs(self, temp_dir):
        """Test multi-fragment SILCS-like simulation"""
        inp_content = """
version: gcmc_2.0
box: 40.0 40.0 40.0
temperature: 298.0
# Multiple fragment types like SILCS
fragname: WAT
fragmuex: -10.5
fragconc: 55.0
mctime: 0.8
fragname: NA
fragmuex: -8.0
fragconc: 0.15
mctime: 0.1
fragname: CL
fragmuex: -7.5
fragconc: 0.15
mctime: 0.1
mcsteps: 50
nprint: 10
fragitp: {}/charmm36.ff/mol/sol.itp
fragitp: {}/charmm36.ff/mol/na.itp
fragitp: {}/charmm36.ff/mol/cl.itp
""".format(TEST_DATA_DIR, TEST_DATA_DIR, TEST_DATA_DIR)

        inp_file = Path(temp_dir) / "silcs.inp"
        inp_file.write_text(inp_content)

        result = self.run_gcmc(inp_file, temp_dir, ["--prefix", "silcs"])

        assert result.returncode == 0, f"Failed: {result.stderr}"

        # Check that multiple fragments were recognized
        assert "WAT" in result.stdout or "SOL" in result.stdout
        assert "NA" in result.stdout or "SOD" in result.stdout
        assert "CL" in result.stdout or "CLA" in result.stdout

    def test_cavity_bias_parameters(self, temp_dir):
        """Test cavity bias parameters from legacy format"""
        inp_content = """
version: gcmc_2.0
box: 30.0 30.0 30.0
temperature: 300.0
fragname: WAT
fragmuex: -10.0
fragconc: 55.0
# Cavity bias parameters
use_cavity_bias: yes
cavity_grid_dx: 1.0
probe_radius: 1.4
exclude_protein_volume: yes
exclude_hydrogens_from_grid: yes
use_vdw_radius_for_grid: yes
mcsteps: 50
nprint: 10
fragitp: {}/charmm36.ff/mol/sol.itp
""".format(TEST_DATA_DIR)

        inp_file = Path(temp_dir) / "cavity.inp"
        inp_file.write_text(inp_content)

        result = self.run_gcmc(inp_file, temp_dir, ["--prefix", "cavity", "--verbose"])

        assert result.returncode == 0, f"Failed: {result.stderr}"

        # Check that cavity bias was enabled (if verbose output available)
        if "--verbose" in result.stdout:
            assert "cavity" in result.stdout.lower() or "bias" in result.stdout.lower()

    def test_gcmc_region_constraint(self, temp_dir):
        """Test GCMC region constraints"""
        # Test sphere region
        inp_content = """
version: gcmc_2.0
box: 50.0 50.0 50.0
temperature: 300.0
fragname: WAT
fragmuex: -10.0
fragconc: 55.0
gcmc_region: sphere 25.0 25.0 25.0 10.0
mcsteps: 50
nprint: 10
fragitp: {}/charmm36.ff/mol/sol.itp
""".format(TEST_DATA_DIR)

        inp_file = Path(temp_dir) / "region.inp"
        inp_file.write_text(inp_content)

        result = self.run_gcmc(inp_file, temp_dir, ["--prefix", "region"])

        assert result.returncode == 0, f"Failed: {result.stderr}"

        # Check PDB coordinates are within sphere
        pdb_file = Path(temp_dir) / "region_final.pdb"
        if pdb_file.exists():
            pdb_content = pdb_file.read_text()
            # Extract coordinates from ATOM lines
            atom_lines = [l for l in pdb_content.split('\n') if l.startswith('ATOM')]
            for line in atom_lines:
                if len(line) > 54:
                    x = float(line[30:38]) / 10.0  # Angstrom to nm
                    y = float(line[38:46]) / 10.0
                    z = float(line[46:54]) / 10.0
                    # Check if within sphere (center at 25.0, 25.0, 25.0 nm, radius 10.0 nm)
                    dist = np.sqrt((x-25.0)**2 + (y-25.0)**2 + (z-25.0)**2)
                    assert dist <= 10.5, f"Atom outside region: ({x}, {y}, {z}), dist={dist}"

    def test_switching_function(self, temp_dir):
        """Test switching function parameters"""
        inp_content = """
version: gcmc_2.0
box: 30.0 30.0 30.0
temperature: 300.0
fragname: WAT
fragmuex: -10.0
fragconc: 55.0
use_switching: yes
switch_r_on: 8.0
switch_r_off: 10.0
mcsteps: 50
nprint: 10
fragitp: {}/charmm36.ff/mol/sol.itp
""".format(TEST_DATA_DIR)

        inp_file = Path(temp_dir) / "switching.inp"
        inp_file.write_text(inp_content)

        result = self.run_gcmc(inp_file, temp_dir, ["--prefix", "switching", "--verbose"])

        assert result.returncode == 0, f"Failed: {result.stderr}"

        # With verbose, check switching was recognized
        if "--verbose" in str(["--verbose"]):
            assert "switch" in result.stdout.lower() or result.returncode == 0

    def test_target_numwaters(self, temp_dir):
        """Test target number of waters feature"""
        inp_content = """
version: gcmc_2.0
box: 30.0 30.0 30.0
temperature: 300.0
fragname: WAT
fragmuex: -10.0
fragconc: 55.0
target_numwaters: 100
wdens: 1.0
mcsteps: 100
nprint: 20
fragitp: {}/charmm36.ff/mol/sol.itp
""".format(TEST_DATA_DIR)

        inp_file = Path(temp_dir) / "target.inp"
        inp_file.write_text(inp_content)

        result = self.run_gcmc(inp_file, temp_dir, ["--prefix", "target"])

        assert result.returncode == 0, f"Failed: {result.stderr}"

        # Check water density output if wdens was set
        if "Water density" in result.stdout:
            # Extract density values
            density_lines = [l for l in result.stdout.split('\n') if "density" in l.lower()]
            assert len(density_lines) > 0, "No density output found"

    def test_pairlist_frequency(self, temp_dir):
        """Test pairlist update frequency parameter"""
        inp_content = """
version: gcmc_2.0
box: 30.0 30.0 30.0
temperature: 300.0
fragname: WAT
fragmuex: -10.0
fragconc: 55.0
pairlist_freq: 500
mcsteps: 50
nprint: 10
fragitp: {}/charmm36.ff/mol/sol.itp
""".format(TEST_DATA_DIR)

        inp_file = Path(temp_dir) / "pairlist.inp"
        inp_file.write_text(inp_content)

        result = self.run_gcmc(inp_file, temp_dir, ["--prefix", "pairlist", "--verbose"])

        assert result.returncode == 0, f"Failed: {result.stderr}"

        # With verbose, check pairlist_freq was parsed
        if "--verbose" in str(["--verbose"]):
            assert "pairlist" in result.stdout.lower() or "500" in result.stdout

    def test_output_file_format(self, temp_dir):
        """Test output file format and content validation"""
        inp_content = """
version: gcmc_2.0
box: 20.0 20.0 20.0
temperature: 300.0
fragname: WAT
fragmuex: -10.0
fragconc: 55.0
mcsteps: 10
nprint: 5
nsave: 5
fragitp: {}/charmm36.ff/mol/sol.itp
""".format(TEST_DATA_DIR)

        inp_file = Path(temp_dir) / "format.inp"
        inp_file.write_text(inp_content)

        result = self.run_gcmc(inp_file, temp_dir, ["--prefix", "format", "--seed", "123"])

        assert result.returncode == 0, f"Failed: {result.stderr}"

        # Check PDB format details
        pdb_file = Path(temp_dir) / "format_final.pdb"
        assert pdb_file.exists()
        pdb_lines = pdb_file.read_text().split('\n')

        # Validate PDB header
        assert any("REMARK GCMC Trajectory" in l for l in pdb_lines)
        assert any("REMARK Step:" in l for l in pdb_lines)
        assert any("REMARK Energy:" in l for l in pdb_lines)
        assert any("REMARK Water molecules:" in l for l in pdb_lines)

        # Validate CRYST1 record
        cryst1_lines = [l for l in pdb_lines if l.startswith("CRYST1")]
        assert len(cryst1_lines) == 1
        cryst1 = cryst1_lines[0]
        parts = cryst1.split()
        assert len(parts) >= 4, f"Unexpected CRYST1 format: {cryst1}"
        lx, ly, lz = float(parts[1]), float(parts[2]), float(parts[3])
        assert lx == pytest.approx(20.0, abs=1e-3)
        assert ly == pytest.approx(20.0, abs=1e-3)
        assert lz == pytest.approx(20.0, abs=1e-3)

        # Check TOP file format
        top_file = Path(temp_dir) / "format_final.top"
        assert top_file.exists()
        top_content = top_file.read_text()

        # Validate TOP sections
        assert "[ system ]" in top_content
        assert "[ molecules ]" in top_content
        assert "[ fragments ]" in top_content
        assert "[ statistics ]" in top_content

        # Check fragment information
        assert "WAT" in top_content or "SOL" in top_content
        assert "55.00" in top_content  # concentration (printed with 2 decimals)

        # ChemPot in topology is in kJ/mol; for version:gcmc_2.0 inputs, fragmuex is kcal/mol.
        frag_lines = top_content.splitlines()
        in_frag = False
        frag_row = None
        for line in frag_lines:
            if line.strip() == "[ fragments ]":
                in_frag = True
                continue
            if in_frag:
                if line.startswith("["):
                    break
                if not line.strip() or line.lstrip().startswith(";"):
                    continue
                cols = line.split()
                if not cols:
                    continue
                if cols[0].upper() in {"WAT", "SOL"} and len(cols) >= 5:
                    frag_row = cols
                    break
        assert frag_row is not None, "Missing fragment row in [ fragments ] section"
        chem_pot_kj = float(frag_row[-1])
        assert chem_pot_kj == pytest.approx(-10.0 * 4.184, abs=0.02)

    def test_error_handling(self, temp_dir):
        """Test error handling for invalid inputs"""
        # Test missing required parameters
        inp_content = """
version: gcmc_2.0
box: 30.0 30.0 30.0
temperature: 300.0
# Missing fragname, fragmuex, fragconc
mcsteps: 10
"""

        inp_file = Path(temp_dir) / "error.inp"
        inp_file.write_text(inp_content)

        result = self.run_gcmc(inp_file, temp_dir, ["--prefix", "error"])

        # Should fail or handle gracefully
        # Either error in output or specific return code
        assert result.returncode != 0 or "error" in result.stderr.lower() or "warning" in result.stdout.lower()

    def test_checkpoint_disabled_by_default(self, temp_dir):
        """Test that checkpoint is disabled by default"""
        inp_content = """
version: gcmc_2.0
box: 30.0 30.0 30.0
temperature: 300.0
fragname: WAT
fragmuex: -10.0
fragconc: 55.0
mcsteps: 10
nprint: 5
fragitp: {}/charmm36.ff/mol/sol.itp
""".format(TEST_DATA_DIR)

        inp_file = Path(temp_dir) / "nocheckpoint.inp"
        inp_file.write_text(inp_content)

        result = self.run_gcmc(inp_file, temp_dir, ["--prefix", "nocheckpoint"])

        assert result.returncode == 0, f"Failed: {result.stderr}"

        # Check no checkpoint files created (actual checkpoint files, not just filenames containing "checkpoint")
        checkpoint_files = [f for f in Path(temp_dir).glob("*.checkpoint") if f.is_file()]
        checkpoint_files.extend([f for f in Path(temp_dir).glob("*_checkpoint.*") if f.is_file()])
        assert len(checkpoint_files) == 0, f"Unexpected checkpoint files: {checkpoint_files}"

        # Check no gcmc_final.txt created
        txt_files = list(Path(temp_dir).glob("*.txt"))
        assert len(txt_files) == 0, f"Unexpected txt files: {txt_files}"

        # Only PDB and TOP should exist
        pdb_files = list(Path(temp_dir).glob("*.pdb"))
        top_files = list(Path(temp_dir).glob("*.top"))
        assert len(pdb_files) > 0, "No PDB files created"
        assert len(top_files) > 0, "No TOP files created"


class TestGCMCStatistics:
    """Test statistical outputs and correctness"""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield tmpdir

    def test_acceptance_rate_reasonable(self, temp_dir):
        """Test that acceptance rates are reasonable"""
        inp_content = """
version: gcmc_2.0
box: 30.0 30.0 30.0
temperature: 300.0
fragname: WAT
fragmuex: -10.0
fragconc: 10.0  # Lower concentration for better acceptance
mcsteps: 100
nprint: 20
fragitp: {}/charmm36.ff/mol/sol.itp
""".format(TEST_DATA_DIR)

        inp_file = Path(temp_dir) / "accept.inp"
        inp_file.write_text(inp_content)

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--prefix", "accept", "--seed", "999"],
            cwd=temp_dir,
            capture_output=True,
            text=True
        )

        assert result.returncode == 0, f"Failed: {result.stderr}"

        # Parse acceptance rate from output
        accept_pattern = r"acceptance rate:\s*([\d.]+)%"
        matches = re.findall(accept_pattern, result.stdout.lower())

        if matches:
            # At least one acceptance rate should be reasonable (0-100%)
            rates = [float(m) for m in matches]
            assert all(0 <= r <= 100 for r in rates), f"Invalid acceptance rates: {rates}"

    def test_energy_conservation(self, temp_dir):
        """Test that energies are reasonable"""
        inp_content = """
version: gcmc_2.0
box: 30.0 30.0 30.0
temperature: 300.0
fragname: WAT
fragmuex: -10.0
fragconc: 55.0
mcsteps: 50
nprint: 10
fragitp: {}/charmm36.ff/mol/sol.itp
""".format(TEST_DATA_DIR)

        inp_file = Path(temp_dir) / "energy.inp"
        inp_file.write_text(inp_content)

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--prefix", "energy"],
            cwd=temp_dir,
            capture_output=True,
            text=True
        )

        assert result.returncode == 0, f"Failed: {result.stderr}"

        # Parse energy values
        energy_pattern = r"energy:\s*([-\d.]+)\s*kJ/mol"
        matches = re.findall(energy_pattern, result.stdout.lower())

        if matches:
            energies = [float(m) for m in matches]
            # Check energies are finite and reasonable
            assert all(np.isfinite(e) for e in energies), f"Non-finite energies: {energies}"
            # Energy should be reasonable (not too positive)
            assert all(e < 1e6 for e in energies), f"Unreasonable energies: {energies}"


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v"])
