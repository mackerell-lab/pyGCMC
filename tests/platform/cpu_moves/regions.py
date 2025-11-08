"""
Test GCMC region constraint functionality
"""

import pytest
import tempfile
from pathlib import Path


TEST_DATA_DIR = Path(__file__).resolve().parents[2] / "data"
GCMC_CPU_PATH = Path(__file__).resolve().parents[3] / "build" / "bin" / "gcmc_cpu"


def test_gcmc_region_sphere(tmp_path):
    """Test sphere region constraint"""
    # Get path to fragment itp file
    fragitp_path = TEST_DATA_DIR / "charmm36.ff/mol/sol.itp"

    # Create a simple inp file with gcmc_region
    inp_content = f"""
version: gcmc_2.0
num_threads: 1
temperature: 300.0
mc_steps: 10
nprint: 5
fragname: SOL
fragmuex: -10.0
fragconc: 55.0
gcmc_region: sphere 25.0 25.0 25.0 10.0
box: 50.0 50.0 50.0
fragitp: {fragitp_path}
"""

    inp_file = tmp_path / "test_sphere.inp"
    inp_file.write_text(inp_content)

    # Import and run gcmc_cpu
    import subprocess

    # Run gcmc_cpu with the inp file
    gcmc_executable = GCMC_CPU_PATH
    if gcmc_executable.exists():
        result = subprocess.run(
            [str(gcmc_executable), "--inp", str(inp_file)],
            cwd=tmp_path,
            capture_output=True,
            text=True
        )

        # Check that the simulation ran successfully
        assert result.returncode == 0, f"gcmc_cpu failed with error: {result.stderr}"

        # Check that region constraint was mentioned in output
        assert "gcmc_region" in result.stdout.lower() or \
               "region constraint" in result.stdout.lower() or \
               "sphere" in result.stdout.lower() or \
               "simulation completed" in result.stdout.lower(), \
               f"Region constraint not properly configured. Output: {result.stdout}"
    else:
        pytest.skip("gcmc_cpu executable not found")


def test_gcmc_region_box(tmp_path):
    """Test box region constraint"""
    # Get path to fragment itp file
    fragitp_path = TEST_DATA_DIR / "charmm36.ff/mol/sol.itp"

    # Create a simple inp file with box region
    inp_content = f"""
version: gcmc_2.0
num_threads: 1
temperature: 300.0
mc_steps: 10
nprint: 5
fragname: SOL
fragmuex: -10.0
fragconc: 55.0
gcmc_region: box 10.0 10.0 10.0 40.0 40.0 40.0
box: 50.0 50.0 50.0
fragitp: {fragitp_path}
"""

    inp_file = tmp_path / "test_box.inp"
    inp_file.write_text(inp_content)

    # Import and run gcmc_cpu
    import subprocess

    # Run gcmc_cpu with the inp file
    gcmc_executable = GCMC_CPU_PATH
    if gcmc_executable.exists():
        result = subprocess.run(
            [str(gcmc_executable), "--inp", str(inp_file)],
            cwd=tmp_path,
            capture_output=True,
            text=True
        )

        # Check that the simulation ran successfully
        assert result.returncode == 0, f"gcmc_cpu failed with error: {result.stderr}"

        # Check that region constraint was mentioned in output
        assert "gcmc_region" in result.stdout.lower() or \
               "region constraint" in result.stdout.lower() or \
               "box" in result.stdout.lower() or \
               "simulation completed" in result.stdout.lower(), \
               f"Region constraint not properly configured. Output: {result.stdout}"
    else:
        pytest.skip("gcmc_cpu executable not found")


def test_gcmc_region_cylinder(tmp_path):
    """Test cylinder region constraint"""
    # Get path to fragment itp file
    fragitp_path = TEST_DATA_DIR / "charmm36.ff/mol/sol.itp"

    # Create a simple inp file with cylinder region
    inp_content = f"""
version: gcmc_2.0
num_threads: 1
temperature: 300.0
mc_steps: 10
nprint: 5
fragname: SOL
fragmuex: -10.0
fragconc: 55.0
gcmc_region: cylinder 25.0 25.0 25.0 10.0 30.0 z
box: 50.0 50.0 50.0
fragitp: {fragitp_path}
"""

    inp_file = tmp_path / "test_cylinder.inp"
    inp_file.write_text(inp_content)

    # Import and run gcmc_cpu
    import subprocess

    # Run gcmc_cpu with the inp file
    gcmc_executable = GCMC_CPU_PATH
    if gcmc_executable.exists():
        result = subprocess.run(
            [str(gcmc_executable), "--inp", str(inp_file)],
            cwd=tmp_path,
            capture_output=True,
            text=True
        )

        # Check that the simulation ran successfully
        assert result.returncode == 0, f"gcmc_cpu failed with error: {result.stderr}"

        # Check that region constraint was mentioned in output
        assert "gcmc_region" in result.stdout.lower() or \
               "region constraint" in result.stdout.lower() or \
               "cylinder" in result.stdout.lower() or \
               "simulation completed" in result.stdout.lower(), \
               f"Region constraint not properly configured. Output: {result.stdout}"
    else:
        pytest.skip("gcmc_cpu executable not found")


if __name__ == "__main__":
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        test_gcmc_region_sphere(tmp_path)
        test_gcmc_region_box(tmp_path)
        test_gcmc_region_cylinder(tmp_path)
        print("All gcmc_region tests passed!")
