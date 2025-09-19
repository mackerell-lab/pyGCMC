"""
Shared fixtures for CLI tests
"""

import pytest
import subprocess
import tempfile
import shutil
from pathlib import Path


def get_gcmc_cpu_path():
    """Get path to gcmc_cpu executable"""
    # Check if we're in build directory
    build_dir = Path(__file__).parent.parent.parent.parent / "build"
    gcmc_cpu = build_dir / "bin" / "gcmc_cpu"
    
    if gcmc_cpu.exists():
        return str(gcmc_cpu)
    
    # Try relative to current directory
    gcmc_cpu = Path("build/bin/gcmc_cpu")
    if gcmc_cpu.exists():
        return str(gcmc_cpu)
    
    # Try in PATH
    result = subprocess.run(["which", "gcmc_cpu"], capture_output=True, text=True)
    if result.returncode == 0:
        return result.stdout.strip()
    
    raise FileNotFoundError("gcmc_cpu executable not found")


@pytest.fixture
def gcmc_cpu():
    """Fixture to get gcmc_cpu path"""
    return get_gcmc_cpu_path()


@pytest.fixture
def test_data_dir():
    """Get test data directory"""
    return Path(__file__).parent.parent.parent / "data"


@pytest.fixture
def temp_dir():
    """Create temporary directory for test outputs"""
    temp = tempfile.mkdtemp(prefix="test_gcmc_cpu_")
    yield temp
    shutil.rmtree(temp, ignore_errors=True)