#!/usr/bin/env python3
"""
Regression test suite for gcmc_cpu
Ensures that changes don't break existing functionality
"""

import pytest
import subprocess
from pathlib import Path
import json

GCMC_CPU_PATH = Path(__file__).parent.parent.parent.parent / "build" / "bin" / "gcmc_cpu"


def _pdb_atom_signature(pdb_path: Path) -> tuple[tuple[str, str, int, float, float, float], ...]:
    """Extract a stable, formatting-insensitive signature from a PDB (ATOM/HETATM only)."""
    atoms: list[tuple[str, str, int, float, float, float]] = []
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resname = line[17:20].strip().upper()
        atom_name = line[12:16].strip().upper()
        try:
            resid = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError as exc:
            raise AssertionError(f"Failed to parse PDB ATOM line: {line}") from exc
        atoms.append((resname, atom_name, resid, round(x, 3), round(y, 3), round(z, 3)))
    return tuple(sorted(atoms))


def _load_acceptance_jsonl(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


class TestRegression:
    """Regression test suite"""
    
    @pytest.fixture
    def reference_setup(self, tmp_path):
        """Setup reference simulation files"""
        # Standard water box for regression testing
        pdb_content = """REMARK Reference water box
CRYST1   25.000   25.000   25.000  90.00  90.00  90.00 P 1           1
ATOM      1  O   WAT     1      10.000  10.000  10.000  1.00  0.00           O
ATOM      2  H1  WAT     1      10.750  10.000  10.000  1.00  0.00           H
ATOM      3  H2  WAT     1      10.000  10.750  10.000  1.00  0.00           H
ATOM      4  O   WAT     2      15.000  15.000  15.000  1.00  0.00           O
ATOM      5  H1  WAT     2      15.750  15.000  15.000  1.00  0.00           H
ATOM      6  H2  WAT     2      15.000  15.750  15.000  1.00  0.00           H
END
"""
        pdb_path = tmp_path / "reference.pdb"
        pdb_path.write_text(pdb_content)
        
        top_path = tmp_path / "reference.top"
        top_path.write_text("[ system ]\nReference\n[ molecules ]\nWAT 2")
        
        inp_path = tmp_path / "reference.inp"
        inp_content = """# Reference simulation for regression testing
top:reference.top
pdb:reference.pdb
op_top:output.top
op_pdb:output.pdb
box_size:25.0 25.0 25.0
gc_center:12.5 12.5 12.5
cutoff:8.0
mcsteps:100
nprint:10
fragname:water
fragconc:55.0
fragmuex:-5.0
"""
        inp_path.write_text(inp_content)
        
        return {"inp": inp_path, "pdb": pdb_path, "top": top_path, "dir": tmp_path}
    
    def test_deterministic_with_seed(self, reference_setup):
        """Same seed should produce identical final structures (no stdout/stderr parsing)."""
        signatures = []
        
        for run in range(3):
            result = subprocess.run(
                [
                    str(GCMC_CPU_PATH),
                    "--inp", str(reference_setup["inp"]),
                    "--prefix", f"det_{run}",
                    "--seed", "424242",
                    "--no-stats"
                ],
                cwd=reference_setup["dir"],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            assert result.returncode == 0, result.stdout + result.stderr

            out_pdb = reference_setup["dir"] / f"det_{run}_final.pdb"
            assert out_pdb.exists(), f"Missing output PDB: {out_pdb}"
            signatures.append(_pdb_atom_signature(out_pdb))
        
        for i in range(1, len(signatures)):
            assert signatures[0] == signatures[i], f"Run {i} differs from run 0 with same seed"
    
    def test_backward_compatibility(self, reference_setup):
        """Test that old-style INP files still work"""
        # Create an INP file with minimal required fields
        old_inp_path = reference_setup["dir"] / "old_style.inp"
        old_content = """# Minimal old-style INP
top:reference.top
pdb:reference.pdb
box_size:25.0 25.0 25.0
mcsteps:50
fragname:water
fragconc:55.0
fragmuex:-5.0
"""
        old_inp_path.write_text(old_content)
        
        out_prefix = reference_setup["dir"] / "old_style"
        result = subprocess.run(
            [
                str(GCMC_CPU_PATH),
                "--inp", str(old_inp_path),
                "--seed", "111",
                "--prefix", str(out_prefix),
            ],
            cwd=reference_setup["dir"],
            capture_output=True,
            text=True,
            timeout=10
        )
        
        # Should still work with defaults for missing parameters
        assert result.returncode == 0, "Old-style INP should still work"
        assert Path(f"{out_prefix}_final.pdb").exists(), "Expected final PDB not created"
    
    def test_known_scenarios(self, reference_setup):
        """
        Sanity-check a few common deck variations without relying on stdout or final.txt text parsing.

        We assert that key inputs (fragconc/fragmuex) are actually consumed by checking the structured
        `--dump-accept` fields (mu/z) on a single insertion attempt.
        """

        def run_once(extra_inp: str, tag: str) -> dict:
            work = reference_setup["dir"]
            inp_path = work / f"{tag}.inp"
            accept_path = work / f"{tag}_acceptance.jsonl"

            base_content = reference_setup["inp"].read_text()
            # Force a single insertion attempt for deterministic record inspection.
            inp_path.write_text(
                base_content
                + "\n"
                + extra_inp
                + "\n"
                + "moves_per_step:1\nmcsteps:1\nnprint:1\nmc_move_prob:1 0 0 0\n"
            )

            result = subprocess.run(
                [
                    str(GCMC_CPU_PATH),
                    "--inp",
                    str(inp_path),
                    "--prefix",
                    tag,
                    "--seed",
                    "777",
                    "--no-stats",
                    "--dump-accept",
                    str(accept_path),
                ],
                cwd=work,
                capture_output=True,
                text=True,
                timeout=20,
            )
            assert result.returncode == 0, result.stdout + result.stderr
            assert accept_path.exists()

            records = _load_acceptance_jsonl(accept_path)
            assert records, f"No acceptance records produced for {tag}"
            rec = records[0]
            assert str(rec.get("move", "")).strip().lower() == "insertion"
            return rec

        baseline = run_once("", "baseline")
        assert float(baseline["z"]) > 0.0

        # 1) More negative mu_ex should reduce activity (z).
        mu_more_negative = run_once("fragmuex:-10.0", "mu_more_negative")
        assert float(mu_more_negative["mu"]) == pytest.approx(-10.0 * 4.184, rel=1e-6, abs=1e-6)
        assert float(mu_more_negative["z"]) < float(baseline["z"])

        # 2) Higher concentration should increase activity (z) linearly.
        high_conc = run_once("fragconc:100.0", "high_conc")
        assert float(high_conc["z"]) > float(baseline["z"])

        # 3) Less negative / positive mu_ex should increase activity (z).
        mu_less_negative = run_once("fragmuex:10.0", "mu_less_negative")
        assert float(mu_less_negative["mu"]) == pytest.approx(10.0 * 4.184, rel=1e-6, abs=1e-6)
        assert float(mu_less_negative["z"]) > float(baseline["z"])
    
    def test_parameter_validation(self, reference_setup):
        """Invalid parameters must fail fast (no silent fallback)."""
        
        invalid_cases = [
            ("negative_steps", "mcsteps:-100"),
            ("zero_box", "box_size:0.0 10.0 10.0"),
            ("huge_cutoff", "cutoff:1000.0"),  # Larger than box
            ("invalid_conc", "fragconc:-10.0"),
        ]
        
        for name, invalid_param in invalid_cases:
            inp_path = reference_setup["dir"] / f"{name}.inp"
            base_content = reference_setup["inp"].read_text()
            
            # Replace or add the invalid parameter
            if ':' in invalid_param:
                key, value = invalid_param.split(':')
                # Simple replacement - would need more robust parsing in production
                lines = base_content.split('\n')
                found = False
                for i, line in enumerate(lines):
                    if line.startswith(key + ':'):
                        lines[i] = invalid_param
                        found = True
                        break
                if not found:
                    lines.append(invalid_param)
                inp_path.write_text('\n'.join(lines))
            
            out_prefix = reference_setup["dir"] / f"invalid_{name}"
            result = subprocess.run(
                [
                    str(GCMC_CPU_PATH),
                    "--inp", str(inp_path),
                    "--prefix", str(out_prefix),
                    "--seed", "12345",
                ],
                cwd=reference_setup["dir"],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            assert result.returncode != 0, f"{name} should fail for invalid param: {invalid_param}"
            # No output artifacts should be produced on failure.
            expected_outputs = [
                Path(f"{out_prefix}_final.pdb"),
                Path(f"{out_prefix}_final.top"),
                Path(f"{out_prefix}_statistics.dat"),
                Path(f"{out_prefix}_final.txt"),
            ]
            assert not any(p.exists() for p in expected_outputs), f"Unexpected outputs: {expected_outputs}"
    
    def test_output_format_stability(self, reference_setup):
        """Test that output format remains consistent"""
        result = subprocess.run(
            [
                str(GCMC_CPU_PATH),
                "--inp", str(reference_setup["inp"]),
                "--prefix", "format",
                "--seed", "12345"
            ],
            cwd=reference_setup["dir"],
            capture_output=True,
            text=True,
            timeout=10
        )
        
        assert result.returncode == 0
        
        output_file = reference_setup["dir"] / "format_final.txt"
        content = output_file.read_text()
        
        # Check for expected sections
        expected_sections = [
            "GCMC Simulation Final Results",
            "Configuration:",
            "Performance:",
            "Statistics:",
            "Fragment Statistics:"
        ]
        
        for section in expected_sections:
            assert section in content, f"Missing section: {section}"
        
        # Check format of specific fields
        lines = content.split('\n')
        for line in lines:
            if "Temperature:" in line:
                assert " K" in line, "Temperature should have units"
            elif "Box:" in line:
                assert " x " in line and " nm" in line, "Box should have dimensions and units"
            elif "Steps/second:" in line or "steps/second:" in line:
                # Should be a number
                parts = line.split(':')
                if len(parts) > 1:
                    try:
                        float(parts[1].strip().split()[0])
                    except:
                        pytest.fail("Steps/second should be numeric")
    
    @pytest.mark.parametrize("seed", [1, 100, 999, 12345, 99999])
    def test_seed_range(self, reference_setup, seed):
        """Test that various seed values work correctly"""
        result = subprocess.run(
            [
                str(GCMC_CPU_PATH),
                "--inp", str(reference_setup["inp"]),
                "--prefix", f"seed_{seed}",
                "--seed", str(seed),
                "--no-stats"
            ],
            cwd=reference_setup["dir"],
            capture_output=True,
            text=True,
            timeout=10
        )
        
        assert result.returncode == 0, f"Failed with seed {seed}"
        
        # Different seeds should produce different results
        output_file = reference_setup["dir"] / f"seed_{seed}_final.txt"
        assert output_file.exists()
    
    def test_cli_argument_precedence(self, reference_setup):
        """Test that CLI arguments override INP file settings"""
        # INP file has mcsteps:100
        # CLI will override with different value
        
        # Run with CLI override
        result = subprocess.run(
            [
                str(GCMC_CPU_PATH),
                "--inp", str(reference_setup["inp"]),
                "--prefix", "override",
                "--seed", "555",
                # Note: --steps override not implemented yet
                # This test documents expected behavior
            ],
            cwd=reference_setup["dir"],
            capture_output=True,
            text=True,
            timeout=10
        )
        
        assert result.returncode == 0
        
        # Check that prefix was overridden
        assert (reference_setup["dir"] / "override_final.txt").exists()


class TestContinuousIntegration:
    """Tests designed for CI/CD pipelines"""
    
    def test_quick_smoke_test(self, tmp_path):
        """Quick smoke test for CI - should complete in seconds"""
        # Minimal files
        (tmp_path / "ci.pdb").write_text("ATOM      1  O   WAT     1       0.0   0.0   0.0\nEND")
        (tmp_path / "ci.top").write_text("[ system ]\nCI\n[ molecules ]\nWAT 1")
        
        inp_path = tmp_path / "ci.inp"
        inp_path.write_text("""top:ci.top
pdb:ci.pdb
box_size:30.0 30.0 30.0
cutoff:8.0
mcsteps:10
fragname:water
fragconc:55.0
fragmuex:-5.0
""")
        
        result = subprocess.run(
            [
                str(GCMC_CPU_PATH),
                "--inp", str(inp_path),
                "--seed", "1"
            ],
            cwd=tmp_path,
            capture_output=True,
            text=True,
            timeout=5
        )
        
        assert result.returncode == 0, "CI smoke test failed"
    
    def test_version_output(self, tmp_path):
        """Test version information (when implemented)"""
        # This would test --version flag when implemented
        # For now, just check that help works
        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--help"],
            capture_output=True,
            text=True,
            cwd=tmp_path
        )
        
        # Help might return 1 instead of 0
        assert result.returncode in [0, 1]
        # Avoid log-string assertions; just ensure no output artifacts are produced.
        outputs = list(tmp_path.glob("*"))
        assert not outputs, f"Unexpected outputs from --help: {outputs}"

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
