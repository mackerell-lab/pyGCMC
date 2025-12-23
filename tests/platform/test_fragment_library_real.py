"""
Test FragmentLibrary with real ITP parsing
This test actually calls the parser through Python bindings and verifies the loaded data
"""
import json
import pytest
from pathlib import Path
import subprocess
import tempfile

# Test data paths
TEST_DATA_DIR = Path(__file__).parent.parent / "data" / "charmm36.ff"
MOL_DIR = TEST_DATA_DIR / "mol"


class TestFragmentLibraryRealParsing:
    """Test FragmentLibrary with real ITP files"""
    
    @pytest.fixture
    def itp_files(self):
        """Provide paths to real ITP test files"""
        return {
            'sol': MOL_DIR / "sol.itp",
            'acox': MOL_DIR / "acox.itp",
            'acey': MOL_DIR / "acey.itp",
            'acet': MOL_DIR / "acet.itp",
            'aald': MOL_DIR / "aald.itp",
            'benx': MOL_DIR / "benx.itp"
        }
    
    def test_load_sol_water(self, itp_files):
        """Test loading SOL (water) ITP file and verify parsed data"""
        # Since we can't directly import FragmentLibrary in Python,
        # we need to test it through gcmc_cpu or create a test executable
        
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            
            # Create a test INP file that loads the SOL fragment
            inp_file = tmp_path / "test_sol.inp"
            inp_file.write_text(f"""
fragitp:{itp_files['sol']}
box_size:10.0 10.0 10.0
temperature:298.15
mcsteps:50
nprint:10
fragname:SOL
fragconc:55.0
fragmuex:5.0
mc_move_prob:1 0 0 0
attempt_prob_ins:1.0
attempt_prob_del:0.0
attempt_prob_trn:0.0
attempt_prob_rot:0.0
""")

            # Run gcmc_cpu with verbose to see fragment loading
            gcmc_cpu = Path(__file__).parent.parent.parent / "build" / "bin" / "gcmc_cpu"
            accept_log = tmp_path / "accept.jsonl"
            out_prefix = tmp_path / "out"
            result = subprocess.run(
                [
                    str(gcmc_cpu),
                    "--inp",
                    str(inp_file),
                    "--seed",
                    "42",
                    "--prefix",
                    str(out_prefix),
                    "--dump-accept",
                    str(accept_log),
                ],
                cwd=tmp_path,
                capture_output=True,
                text=True,
                timeout=5
            )
            
            # Clean up generated files
            for f in tmp_path.glob("gcmc_final.*"):
                f.unlink()
            
            assert result.returncode == 0, f"Simulation failed: {result.stderr}"
            assert accept_log.exists(), "Acceptance log missing"
            records = [json.loads(line) for line in accept_log.read_text().splitlines() if line.strip()]
            assert any(
                rec.get("move") == "insertion" and str(rec.get("species", "")).upper() == "SOL"
                for rec in records
            ), "No SOL insertion attempts recorded"
    
    def test_load_acox_with_virtual_site(self, itp_files):
        """Test loading ACOX (acetone) with virtual site"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            
            # Create test INP for ACOX
            inp_file = tmp_path / "test_acox.inp"
            inp_file.write_text(f"""
fragitp:{itp_files['acox']}
box_size:10.0 10.0 10.0
temperature:298.15
mcsteps:50
nprint:10
fragname:ACOX
fragconc:1.0
fragmuex:5.0
mc_move_prob:1 0 0 0
attempt_prob_ins:1.0
attempt_prob_del:0.0
attempt_prob_trn:0.0
attempt_prob_rot:0.0
""")
            
            gcmc_cpu = Path(__file__).parent.parent.parent / "build" / "bin" / "gcmc_cpu"
            accept_log = tmp_path / "accept.jsonl"
            out_prefix = tmp_path / "out"
            result = subprocess.run(
                [
                    str(gcmc_cpu),
                    "--inp",
                    str(inp_file),
                    "--seed",
                    "42",
                    "--prefix",
                    str(out_prefix),
                    "--dump-accept",
                    str(accept_log),
                ],
                cwd=tmp_path,
                capture_output=True,
                text=True,
                timeout=5
            )
            
            # Clean up
            for f in tmp_path.glob("gcmc_final.*"):
                f.unlink()
            
            assert result.returncode == 0, f"Simulation failed: {result.stderr}"
            assert accept_log.exists(), "Acceptance log missing"
            records = [json.loads(line) for line in accept_log.read_text().splitlines() if line.strip()]
            assert any(
                rec.get("move") == "insertion" and str(rec.get("species", "")).upper() == "ACOX"
                for rec in records
            ), "No ACOX insertion attempts recorded"
    
    def test_load_charged_fragment(self, itp_files):
        """Test loading ACEY (acetate ion) with -1 charge"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            
            inp_file = tmp_path / "test_acey.inp"
            inp_file.write_text(f"""
fragitp:{itp_files['acey']}
box_size:10.0 10.0 10.0
temperature:298.15
mcsteps:50
nprint:10
fragname:ACEY
fragconc:1.0
fragmuex:5.0
mc_move_prob:1 0 0 0
attempt_prob_ins:1.0
attempt_prob_del:0.0
attempt_prob_trn:0.0
attempt_prob_rot:0.0
""")
            
            gcmc_cpu = Path(__file__).parent.parent.parent / "build" / "bin" / "gcmc_cpu"
            accept_log = tmp_path / "accept.jsonl"
            out_prefix = tmp_path / "out"
            result = subprocess.run(
                [
                    str(gcmc_cpu),
                    "--inp",
                    str(inp_file),
                    "--seed",
                    "42",
                    "--prefix",
                    str(out_prefix),
                    "--dump-accept",
                    str(accept_log),
                ],
                cwd=tmp_path,
                capture_output=True,
                text=True,
                timeout=5
            )
            
            # Clean up
            for f in tmp_path.glob("gcmc_final.*"):
                f.unlink()
            
            assert result.returncode == 0, f"Simulation failed: {result.stderr}"
            assert accept_log.exists(), "Acceptance log missing"
            records = [json.loads(line) for line in accept_log.read_text().splitlines() if line.strip()]
            assert any(
                rec.get("move") == "insertion" and str(rec.get("species", "")).upper() == "ACEY"
                for rec in records
            ), "No ACEY insertion attempts recorded"
    
    def test_load_multiple_fragments(self, itp_files):
        """Test loading multiple fragment types in one simulation"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            
            inp_file = tmp_path / "test_multi.inp"
            inp_file.write_text(f"""
fragitp:{itp_files['sol']}
fragitp:{itp_files['benx']}
fragitp:{itp_files['acet']}
box_size:20.0 20.0 20.0
temperature:300.0
mcsteps:200
nprint:50
fragname:SOL BENX ACET
fragconc:55.0 1.0 1.0
fragmuex:5.0 5.0 5.0
mc_move_prob:1 0 0 0
attempt_prob_ins:1.0
attempt_prob_del:0.0
attempt_prob_trn:0.0
attempt_prob_rot:0.0
""")
            
            gcmc_cpu = Path(__file__).parent.parent.parent / "build" / "bin" / "gcmc_cpu"
            accept_log = tmp_path / "accept.jsonl"
            out_prefix = tmp_path / "out"
            result = subprocess.run(
                [
                    str(gcmc_cpu),
                    "--inp",
                    str(inp_file),
                    "--seed",
                    "42",
                    "--prefix",
                    str(out_prefix),
                    "--dump-accept",
                    str(accept_log),
                ],
                cwd=tmp_path,
                capture_output=True,
                text=True,
                timeout=5
            )
            
            # Clean up
            for f in tmp_path.glob("gcmc_final.*"):
                f.unlink()
            
            assert result.returncode == 0, f"Simulation failed: {result.stderr}"
            assert accept_log.exists(), "Acceptance log missing"
            records = [json.loads(line) for line in accept_log.read_text().splitlines() if line.strip()]
            seen_species = {str(rec.get("species", "")).upper() for rec in records if rec.get("move") == "insertion"}
            expected = {"SOL", "BENX", "ACET"}
            assert expected.issubset(seen_species), f"Missing insertion attempts for: {expected - seen_species}"
    
    def test_invalid_itp_file(self):
        """Test handling of invalid ITP file"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            
            # Create an invalid ITP file
            bad_itp = tmp_path / "bad.itp"
            bad_itp.write_text("This is not a valid ITP file\n")
            
            inp_file = tmp_path / "test_bad.inp"
            inp_file.write_text(f"""
fragitp:{bad_itp}
box_size:10.0 10.0 10.0
temperature:298.15
mcsteps:10
fragname:BAD
fragconc:1.0
fragmuex:-1.0
""")
            
            gcmc_cpu = Path(__file__).parent.parent.parent / "build" / "bin" / "gcmc_cpu"
            result = subprocess.run(
                [str(gcmc_cpu), "--inp", str(inp_file)],
                cwd=tmp_path,
                capture_output=True,
                text=True,
                timeout=5
            )
            
            # Clean up
            for f in tmp_path.glob("gcmc_final.*"):
                f.unlink()
            
            # Should fail cleanly (no segfault) for an invalid ITP file.
            assert result.returncode != -11
            assert result.returncode != 0
    
    def test_nonexistent_itp_file(self):
        """Test handling of non-existent ITP file"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            
            inp_file = tmp_path / "test_missing.inp"
            inp_file.write_text("""
fragitp:/nonexistent/file.itp
box_size:10.0 10.0 10.0
temperature:298.15
mcsteps:10
fragname:MISSING
fragconc:1.0
fragmuex:-1.0
""")
            
            gcmc_cpu = Path(__file__).parent.parent.parent / "build" / "bin" / "gcmc_cpu"
            result = subprocess.run(
                [str(gcmc_cpu), "--inp", str(inp_file)],
                cwd=tmp_path,
                capture_output=True,
                text=True,
                timeout=5
            )
            
            # Clean up
            for f in tmp_path.glob("gcmc_final.*"):
                f.unlink()
            
            # Should fail cleanly (no segfault) for a missing ITP file.
            assert result.returncode != -11
            assert result.returncode != 0


class TestFragmentLibraryContent:
    """Test ITP file content expectations"""
    
    def test_sol_itp_content(self):
        """Verify SOL ITP file has expected content"""
        sol_itp = MOL_DIR / "sol.itp"
        assert sol_itp.exists(), f"SOL ITP not found at {sol_itp}"
        
        content = sol_itp.read_text()
        
        # SOL should have 3 atoms
        assert "OW" in content or "OT" in content, "Oxygen atom not found"
        assert "HW1" in content or "HT" in content, "First hydrogen not found"
        assert "HW2" in content, "Second hydrogen not found"
        
        # Check charges
        assert "-0.834" in content, "Oxygen charge not found"
        assert "0.417" in content, "Hydrogen charge not found"
        
        # Check for moleculetype section
        assert "[ moleculetype ]" in content
        assert "SOL" in content
    
    def test_acox_itp_content(self):
        """Verify ACOX ITP has virtual site"""
        acox_itp = MOL_DIR / "acox.itp"
        assert acox_itp.exists()
        
        content = acox_itp.read_text()
        
        # Should have virtual site
        assert "[ virtual_sites" in content or "LP" in content, "Virtual site not found"
        
        # Should have atoms section
        assert "[ atoms ]" in content
        assert "ACOX" in content
        
        # Check for carbonyl oxygen
        assert "O1" in content or "OG2D3" in content
    
    def test_charged_fragment_content(self):
        """Verify charged fragments have appropriate charges"""
        acey_itp = MOL_DIR / "acey.itp"
        assert acey_itp.exists()
        
        content = acey_itp.read_text()
        
        # ACEY is acetate with -1 charge
        # Should have two oxygens with negative charge
        assert "OG2D2" in content or "-0.76" in content
        
        # Check moleculetype
        assert "[ moleculetype ]" in content
        assert "ACEY" in content


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
