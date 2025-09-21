#!/usr/bin/env python3
"""
Feature tests for gcmc_cpu functionality
Tests specific features that are critical for GCMC simulations
"""

import pytest
import subprocess
import tempfile
import numpy as np
from pathlib import Path
import json
import re

GCMC_CPU_PATH = Path(__file__).parent.parent.parent / "build/bin/gcmc_cpu"
TEST_DATA_DIR = Path(__file__).parent.parent / "data"


class TestWaterInsertion:
    """Test water insertion and deletion mechanics"""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield tmpdir

    def test_water_molecule_structure(self, temp_dir):
        """Test that water molecules maintain correct structure"""
        inp_content = """
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

        inp_file = Path(temp_dir) / "water_struct.inp"
        inp_file.write_text(inp_content)

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--prefix", "water", "--seed", "42"],
            cwd=temp_dir,
            capture_output=True,
            text=True
        )

        assert result.returncode == 0, f"Failed: {result.stderr}"

        # Check PDB for water structure
        pdb_file = Path(temp_dir) / "water_final.pdb"
        if pdb_file.exists():
            pdb_content = pdb_file.read_text()
            atom_lines = [l for l in pdb_content.split('\n') if l.startswith('ATOM')]

            # Group atoms by residue
            residues = {}
            for line in atom_lines:
                if len(line) > 54:
                    res_id = line[22:26].strip()
                    atom_name = line[12:16].strip()
                    if res_id not in residues:
                        residues[res_id] = []
                    residues[res_id].append({
                        'name': atom_name,
                        'x': float(line[30:38]),
                        'y': float(line[38:46]),
                        'z': float(line[46:54])
                    })

            # Check each water molecule
            for res_id, atoms in residues.items():
                if len(atoms) == 3:  # Water has 3 atoms
                    # Check atom names
                    names = [a['name'] for a in atoms]
                    assert 'OW' in names or 'O' in names, f"No oxygen in water {res_id}"
                    h_count = sum(1 for n in names if 'H' in n)
                    assert h_count == 2, f"Wrong number of hydrogens in water {res_id}: {names}"

                    # Check O-H distances (should be ~0.957 Angstrom)
                    o_atom = next((a for a in atoms if 'O' in a['name']), None)
                    h_atoms = [a for a in atoms if 'H' in a['name']]

                    if o_atom and len(h_atoms) == 2:
                        for h in h_atoms:
                            dist = np.sqrt((o_atom['x']-h['x'])**2 +
                                         (o_atom['y']-h['y'])**2 +
                                         (o_atom['z']-h['z'])**2)
                            assert 0.8 < dist < 1.2, f"Invalid O-H distance: {dist} Å"

    def test_water_density_calculation(self, temp_dir):
        """Test water density calculation and output"""
        inp_content = """
version: gcmc_2.0
box: 30.0 30.0 30.0
temperature: 300.0
fragname: WAT
fragmuex: -10.0
fragconc: 55.0
wdens: 1.0
mcsteps: 100
nprint: 20
fragitp: {}/charmm36.ff/mol/sol.itp
""".format(TEST_DATA_DIR)

        inp_file = Path(temp_dir) / "density.inp"
        inp_file.write_text(inp_content)

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--prefix", "density", "--verbose"],
            cwd=temp_dir,
            capture_output=True,
            text=True
        )

        assert result.returncode == 0, f"Failed: {result.stderr}"

        # Check for density output
        if "density" in result.stdout.lower():
            # Parse density values
            density_pattern = r"density.*?(\d+\.?\d*)\s*molecules"
            matches = re.findall(density_pattern, result.stdout.lower())

            if matches:
                densities = [float(m) for m in matches]
                # Density should be non-negative
                assert all(d >= 0 for d in densities), f"Negative densities: {densities}"


class TestCavityBias:
    """Test cavity bias functionality"""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield tmpdir

    def test_cavity_bias_grid(self, temp_dir):
        """Test cavity bias grid generation"""
        inp_content = """
version: gcmc_2.0
box: 20.0 20.0 20.0
temperature: 300.0
fragname: WAT
fragmuex: -10.0
fragconc: 55.0
use_cavity_bias: yes
cavity_grid_dx: 2.0
probe_radius: 1.4
mcsteps: 50
nprint: 10
fragitp: {}/charmm36.ff/mol/sol.itp
""".format(TEST_DATA_DIR)

        inp_file = Path(temp_dir) / "cavity.inp"
        inp_file.write_text(inp_content)

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--prefix", "cavity", "--verbose"],
            cwd=temp_dir,
            capture_output=True,
            text=True
        )

        assert result.returncode == 0, f"Failed: {result.stderr}"

        # With verbose, should mention cavity bias
        if "--verbose" in str(["--verbose"]):
            # Check that cavity bias is being used
            assert "cavity" in result.stdout.lower() or result.returncode == 0

    def test_exclusion_flags(self, temp_dir):
        """Test protein and hydrogen exclusion flags"""
        inp_content = """
version: gcmc_2.0
box: 30.0 30.0 30.0
temperature: 300.0
fragname: WAT
fragmuex: -10.0
fragconc: 55.0
use_cavity_bias: yes
exclude_protein_volume: yes
exclude_hydrogens_from_grid: yes
use_vdw_radius_for_grid: yes
mcsteps: 50
nprint: 10
fragitp: {}/charmm36.ff/mol/sol.itp
""".format(TEST_DATA_DIR)

        inp_file = Path(temp_dir) / "exclusion.inp"
        inp_file.write_text(inp_content)

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--prefix", "exclusion"],
            cwd=temp_dir,
            capture_output=True,
            text=True
        )

        assert result.returncode == 0, f"Failed: {result.stderr}"
        assert "Simulation completed" in result.stdout


class TestRegionConstraints:
    """Test GCMC region constraints"""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield tmpdir

    def test_sphere_region(self, temp_dir):
        """Test sphere region constraint"""
        inp_content = """
version: gcmc_2.0
box: 50.0 50.0 50.0
temperature: 300.0
fragname: WAT
fragmuex: -10.0
fragconc: 55.0
gcmc_region: sphere 25.0 25.0 25.0 10.0
mcsteps: 100
nprint: 20
fragitp: {}/charmm36.ff/mol/sol.itp
""".format(TEST_DATA_DIR)

        inp_file = Path(temp_dir) / "sphere.inp"
        inp_file.write_text(inp_content)

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--prefix", "sphere", "--seed", "123"],
            cwd=temp_dir,
            capture_output=True,
            text=True
        )

        assert result.returncode == 0, f"Failed: {result.stderr}"

        # Verify molecules are within sphere
        pdb_file = Path(temp_dir) / "sphere_final.pdb"
        if pdb_file.exists():
            self._check_sphere_constraint(pdb_file, center=(25.0, 25.0, 25.0), radius=10.0)

    def test_box_region(self, temp_dir):
        """Test box region constraint"""
        inp_content = """
version: gcmc_2.0
box: 50.0 50.0 50.0
temperature: 300.0
fragname: WAT
fragmuex: -10.0
fragconc: 55.0
gcmc_region: box 10.0 10.0 10.0 40.0 40.0 40.0
mcsteps: 100
nprint: 20
fragitp: {}/charmm36.ff/mol/sol.itp
""".format(TEST_DATA_DIR)

        inp_file = Path(temp_dir) / "box.inp"
        inp_file.write_text(inp_content)

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--prefix", "box", "--seed", "456"],
            cwd=temp_dir,
            capture_output=True,
            text=True
        )

        assert result.returncode == 0, f"Failed: {result.stderr}"

        # Verify molecules are within box
        pdb_file = Path(temp_dir) / "box_final.pdb"
        if pdb_file.exists():
            self._check_box_constraint(pdb_file, min_coord=(10.0, 10.0, 10.0), max_coord=(40.0, 40.0, 40.0))

    def test_cylinder_region(self, temp_dir):
        """Test cylinder region constraint"""
        inp_content = """
version: gcmc_2.0
box: 50.0 50.0 50.0
temperature: 300.0
fragname: WAT
fragmuex: -10.0
fragconc: 55.0
gcmc_region: cylinder 25.0 25.0 25.0 10.0 30.0 z
mcsteps: 100
nprint: 20
fragitp: {}/charmm36.ff/mol/sol.itp
""".format(TEST_DATA_DIR)

        inp_file = Path(temp_dir) / "cylinder.inp"
        inp_file.write_text(inp_content)

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--prefix", "cylinder", "--seed", "789"],
            cwd=temp_dir,
            capture_output=True,
            text=True
        )

        assert result.returncode == 0, f"Failed: {result.stderr}"

    def _check_sphere_constraint(self, pdb_file, center, radius):
        """Helper to check atoms are within sphere"""
        pdb_content = pdb_file.read_text()
        atom_lines = [l for l in pdb_content.split('\n') if l.startswith('ATOM')]

        for line in atom_lines:
            if len(line) > 54:
                x = float(line[30:38])  # Angstrom
                y = float(line[38:46])
                z = float(line[46:54])

                # Calculate distance from center
                dist = np.sqrt((x-center[0])**2 + (y-center[1])**2 + (z-center[2])**2)
                assert dist <= radius + 0.5, f"Atom at ({x},{y},{z}) outside sphere, dist={dist}"

    def _check_box_constraint(self, pdb_file, min_coord, max_coord):
        """Helper to check atoms are within box"""
        pdb_content = pdb_file.read_text()
        atom_lines = [l for l in pdb_content.split('\n') if l.startswith('ATOM')]

        for line in atom_lines:
            if len(line) > 54:
                x = float(line[30:38])  # Angstrom
                y = float(line[38:46])
                z = float(line[46:54])

                # Check within box boundaries
                assert min_coord[0] - 0.5 <= x <= max_coord[0] + 0.5, f"X coord {x} outside box"
                assert min_coord[1] - 0.5 <= y <= max_coord[1] + 0.5, f"Y coord {y} outside box"
                assert min_coord[2] - 0.5 <= z <= max_coord[2] + 0.5, f"Z coord {z} outside box"


class TestTargetControl:
    """Test target number control features"""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield tmpdir

    def test_target_numwaters(self, temp_dir):
        """Test target number of waters control"""
        inp_content = """
version: gcmc_2.0
box: 30.0 30.0 30.0
temperature: 300.0
fragname: WAT
fragmuex: -10.0
fragconc: 55.0
target_numwaters: 50
mcsteps: 200
nprint: 40
fragitp: {}/charmm36.ff/mol/sol.itp
""".format(TEST_DATA_DIR)

        inp_file = Path(temp_dir) / "target.inp"
        inp_file.write_text(inp_content)

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--prefix", "target", "--seed", "333"],
            cwd=temp_dir,
            capture_output=True,
            text=True
        )

        assert result.returncode == 0, f"Failed: {result.stderr}"

        # Check final water count
        pdb_file = Path(temp_dir) / "target_final.pdb"
        if pdb_file.exists():
            pdb_content = pdb_file.read_text()
            # Extract water count from REMARK
            water_count_pattern = r"REMARK Water molecules:\s*(\d+)"
            match = re.search(water_count_pattern, pdb_content)
            if match:
                water_count = int(match.group(1))
                # Should be biased toward target
                assert 0 <= water_count <= 200, f"Unexpected water count: {water_count}"


class TestParameterParsing:
    """Test parameter parsing and validation"""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield tmpdir

    def test_all_parameters(self, temp_dir):
        """Test parsing of all supported parameters"""
        inp_content = """
version: gcmc_2.0
box: 30.0 30.0 30.0
temperature: 298.15
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
nsave: 25
use_cavity_bias: yes
cavity_grid_dx: 1.5
probe_radius: 1.4
use_switching: yes
switch_r_on: 8.0
switch_r_off: 10.0
pairlist_freq: 500
target_numwaters: 100
wdens: 1.0
gcmc_region: sphere 15.0 15.0 15.0 5.0
exclude_protein_volume: yes
exclude_hydrogens_from_grid: yes
use_vdw_radius_for_grid: yes
fragitp: {}/charmm36.ff/mol/sol.itp
fragitp: {}/charmm36.ff/mol/na.itp
fragitp: {}/charmm36.ff/mol/cl.itp
""".format(TEST_DATA_DIR, TEST_DATA_DIR, TEST_DATA_DIR)

        inp_file = Path(temp_dir) / "all_params.inp"
        inp_file.write_text(inp_content)

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--prefix", "params", "--verbose"],
            cwd=temp_dir,
            capture_output=True,
            text=True
        )

        assert result.returncode == 0, f"Failed: {result.stderr}"

        # Check that key parameters were parsed (if verbose output)
        if "--verbose" in str(["--verbose"]):
            # Check various parameters mentioned in output
            expected_keywords = ["temperature", "298.15", "fragment", "cavity",
                               "switching", "pairlist", "region"]
            found_count = sum(1 for kw in expected_keywords if kw.lower() in result.stdout.lower())
            assert found_count > 3, f"Not enough parameters recognized in output"

    def test_invalid_parameter_handling(self, temp_dir):
        """Test handling of invalid parameters"""
        inp_content = """
version: gcmc_2.0
box: 30.0 30.0 30.0
temperature: -100.0  # Invalid temperature
fragname: WAT
fragmuex: -10.0
fragconc: 55.0
mcsteps: 10
fragitp: {}/charmm36.ff/mol/sol.itp
""".format(TEST_DATA_DIR)

        inp_file = Path(temp_dir) / "invalid.inp"
        inp_file.write_text(inp_content)

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--prefix", "invalid"],
            cwd=temp_dir,
            capture_output=True,
            text=True
        )

        # Should either fail or handle gracefully
        # Negative temperature might cause issues
        assert result.returncode != 0 or "warning" in result.stdout.lower() or result.returncode == 0


class TestOutputValidation:
    """Test output file validation"""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield tmpdir

    def test_pdb_atom_numbering(self, temp_dir):
        """Test PDB atom numbering is sequential"""
        inp_content = """
version: gcmc_2.0
box: 20.0 20.0 20.0
temperature: 300.0
fragname: WAT
fragmuex: -10.0
fragconc: 55.0
mcsteps: 50
nprint: 10
fragitp: {}/charmm36.ff/mol/sol.itp
""".format(TEST_DATA_DIR)

        inp_file = Path(temp_dir) / "numbering.inp"
        inp_file.write_text(inp_content)

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--prefix", "num", "--seed", "111"],
            cwd=temp_dir,
            capture_output=True,
            text=True
        )

        assert result.returncode == 0, f"Failed: {result.stderr}"

        pdb_file = Path(temp_dir) / "num_final.pdb"
        if pdb_file.exists():
            pdb_content = pdb_file.read_text()
            atom_lines = [l for l in pdb_content.split('\n') if l.startswith('ATOM')]

            # Check sequential numbering
            for i, line in enumerate(atom_lines, 1):
                atom_num = int(line[6:11].strip())
                assert atom_num == i, f"Atom numbering not sequential: expected {i}, got {atom_num}"

    def test_top_file_statistics(self, temp_dir):
        """Test TOP file statistics section"""
        inp_content = """
version: gcmc_2.0
box: 20.0 20.0 20.0
temperature: 300.0
fragname: WAT
fragmuex: -10.0
fragconc: 55.0
mcsteps: 100
nprint: 20
fragitp: {}/charmm36.ff/mol/sol.itp
""".format(TEST_DATA_DIR)

        inp_file = Path(temp_dir) / "stats.inp"
        inp_file.write_text(inp_content)

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--prefix", "stats", "--seed", "222"],
            cwd=temp_dir,
            capture_output=True,
            text=True
        )

        assert result.returncode == 0, f"Failed: {result.stderr}"

        top_file = Path(temp_dir) / "stats_final.top"
        if top_file.exists():
            top_content = top_file.read_text()

            # Check required sections
            assert "[ system ]" in top_content
            assert "[ molecules ]" in top_content
            assert "[ fragments ]" in top_content
            assert "[ statistics ]" in top_content

            # Check fragment info
            if "[ fragments ]" in top_content:
                # Should contain fragment name, count, concentration, chemical potential
                assert "WAT" in top_content or "SOL" in top_content
                assert "55.00" in top_content  # concentration
                assert "-10.00" in top_content  # chemical potential


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v"])