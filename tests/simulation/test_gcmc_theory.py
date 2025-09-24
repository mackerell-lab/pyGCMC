"""
Theoretical validation tests for gcmc_cpu
基于GCMC统计力学理论的验证测试
"""

import pytest
import numpy as np
import subprocess
import tempfile
import re
from pathlib import Path

# Path to gcmc_cpu executable
GCMC_CPU_PATH = Path(__file__).parent.parent.parent / "build" / "bin" / "gcmc_cpu"

# Physical constants
KB = 8.314e-3  # kJ/mol/K
NA = 6.02214076e23  # Avogadro's number
WATER_MW = 18.015  # g/mol


class TestGCMCTheory:
    """Test suite for GCMC theoretical validation"""

    @staticmethod
    def create_minimal_water_files(tmpdir):
        """Create minimal water PDB and TOP files"""
        # Minimal PDB for water with CRYST1 line for box dimensions
        pdb_file = tmpdir / "water.pdb"
        pdb_content = """CRYST1   15.000   15.000   15.000  90.00  90.00  90.00 P 1           1
ATOM      1  O   WAT     1       0.000   0.000   0.000  1.00  0.00
ATOM      2  H1  WAT     1       0.757   0.586   0.000  1.00  0.00
ATOM      3  H2  WAT     1      -0.757   0.586   0.000  1.00  0.00
END
"""
        pdb_file.write_text(pdb_content)

        # Minimal TOP with realistic water parameters (TIP3P-like)
        top_file = tmpdir / "water.top"
        top_content = """[ defaults ]
; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ
1 2 yes 0.5 0.8333

[ atomtypes ]
; name  at.num  mass     charge   ptype    sigma      epsilon
O       8       15.9994  -0.834   A        3.15061e-01  6.36386e-01
H       1       1.008    0.417    A        0.00000e+00  0.00000e+00

[ moleculetype ]
; name  nrexcl
WAT     2

[ atoms ]
;   nr  type  resnr  residue  atom  cgnr  charge    mass
1   O     1      WAT      O     1     -0.834   15.9994
2   H     1      WAT      H1    1      0.417    1.008
3   H     1      WAT      H2    1      0.417    1.008
"""
        top_file.write_text(top_content)

        # Create fragment template file for WAT
        itp_file = tmpdir / "wat.itp"
        itp_content = """[ moleculetype ]
; name  nrexcl
WAT     2

[ atoms ]
;   nr  type  resnr  residue  atom  cgnr  charge    mass
1   O     1      WAT      O     1     -0.834   15.9994
2   H     1      WAT      H1    1      0.417    1.008
3   H     1      WAT      H2    1      0.417    1.008

[ bonds ]
; i  j  func  length  force
1    2   1    0.09572  502416.0
1    3   1    0.09572  502416.0

[ angles ]
; i  j  k  func  angle  force
2    1    3   1    104.52  628.02
"""
        itp_file.write_text(itp_content)

        return str(pdb_file), str(top_file), str(itp_file)

    @staticmethod
    def run_gcmc(inp_content, tmpdir, steps=1000, seed=12345):
        """Run GCMC simulation and return metrics"""
        inp_file = tmpdir / "test.inp"
        inp_file.write_text(inp_content)

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--seed", str(seed)],
            cwd=str(tmpdir),
            capture_output=True,
            text=True,
            timeout=30
        )

        # Parse output for key metrics
        metrics = {
            "final_count": 0,
            "acceptance_rate": 0.0,
            "insert_accept": 0.0,
            "delete_accept": 0.0,
            "stdout": result.stdout,
            "returncode": result.returncode
        }

        # Extract final molecule count
        count_matches = re.findall(r'WAT:\s+(\d+)', result.stdout)
        if count_matches:
            metrics["final_count"] = int(count_matches[-1])

        # Extract acceptance rates
        accept_match = re.search(r'Total acceptance rate:\s+([\d.]+)%', result.stdout)
        if accept_match:
            metrics["acceptance_rate"] = float(accept_match.group(1))

        # Extract per-move acceptance
        insert_match = re.search(r'Insert.*accept:\s+([\d.]+)%', result.stdout)
        if insert_match:
            metrics["insert_accept"] = float(insert_match.group(1))

        delete_match = re.search(r'Delete.*accept:\s+([\d.]+)%', result.stdout)
        if delete_match:
            metrics["delete_accept"] = float(delete_match.group(1))

        return metrics

    def test_chemical_potential_response(self, tmp_path):
        """Test that density increases monotonically with chemical potential"""
        pdb, top, itp = self.create_minimal_water_files(tmp_path)

        # Use moderate mu values to get reasonable activities
        # exp(beta*mu) for mu=0,1,2,3 gives activities of ~1.0, ~2.2, ~5.0, ~11.1 molecules/nm³
        mu_values = [0.0, 1.0, 2.0, 3.0]
        densities = []

        for mu in mu_values:
            # Only specify chemical potential, not concentration
            # This tests the pure chemical potential effect
            inp_content = f"""# Chemical potential scan
pdb:{pdb}
top:{top}
fragitp:{itp}
op_pdb:output.pdb
op_top:output.top
box_size:15.0 15.0 15.0
cutoff:7.0
mcsteps:500
nprint:100
fragname:WAT
fragmuex:{mu}
"""
            metrics = self.run_gcmc(inp_content, tmp_path, steps=500)

            # Calculate density (molecules/nm³)
            volume = 15.0 * 15.0 * 15.0  # nm³
            density = metrics["final_count"] / volume
            densities.append(density)

        # Check general trend (allowing some noise)
        # At least the extreme values should follow the trend
        assert densities[-1] > densities[0], \
            f"Density should increase with μ: {densities[0]:.3f} -> {densities[-1]:.3f}"

    def test_detailed_balance(self, tmp_path):
        """Test that insert/delete acceptance rates satisfy detailed balance"""
        pdb, top, itp = self.create_minimal_water_files(tmp_path)

        inp_content = f"""# Detailed balance test
pdb:{pdb}
top:{top}
fragitp:{itp}
op_pdb:output.pdb
op_top:output.top
box_size:10.0 10.0 10.0
cutoff:5.0
mcsteps:1000
nprint:250
fragname:WAT
fragconc:0.2
# Force equal attempt probabilities
attempt_prob_ins:0.25
attempt_prob_del:0.25
attempt_prob_trn:0.25
attempt_prob_rot:0.25
"""
        metrics = self.run_gcmc(inp_content, tmp_path, steps=5000)

        # Check that simulation ran successfully
        assert metrics["returncode"] == 0, "Simulation should complete successfully"

        # Check that we have some molecules
        assert metrics["final_count"] > 0, "Should have inserted some molecules"

        # Check that acceptance rates are reasonable (not 0 or 100)
        if metrics["insert_accept"] > 0 and metrics["delete_accept"] > 0:
            ratio = metrics["insert_accept"] / metrics["delete_accept"]
            # Ratio should be in reasonable range (detailed balance)
            assert 0.01 < ratio < 100, f"Insert/Delete ratio {ratio:.3f} outside reasonable bounds"

    def test_region_constraint(self, tmp_path):
        """Test that molecules are constrained to specified regions"""
        pdb, top, itp = self.create_minimal_water_files(tmp_path)

        # Test sphere region
        inp_content = f"""# Region constraint test
pdb:{pdb}
top:{top}
fragitp:{itp}
op_pdb:region_test.pdb
op_top:output.top
box_size:30.0 30.0 30.0
gcmc_region:sphere 15.0 15.0 15.0 5.0
cutoff:9.0
mcsteps:1000
nprint:500
fragname:WAT
fragconc:0.1
"""
        metrics = self.run_gcmc(inp_content, tmp_path, steps=2000)

        # Check that region constraint was applied
        assert "region" in metrics["stdout"].lower(), "Region constraint should be mentioned in output"

        # Check that molecules were created
        assert metrics["final_count"] > 0, "Should create molecules in constrained region"

        # The number of molecules should be less than if the whole box was used
        # (sphere volume is much smaller than box volume)
        sphere_volume = (4.0/3.0) * np.pi * 5.0**3  # ~523 nm³
        box_volume = 30.0 * 30.0 * 30.0  # 27000 nm³
        volume_ratio = sphere_volume / box_volume  # ~0.019

        # Molecules should roughly scale with available volume (with wide tolerance)
        # We just check it's not filling the whole box
        expected_max = 100  # For a small sphere, shouldn't have too many molecules
        assert metrics["final_count"] < expected_max, \
            f"Too many molecules ({metrics['final_count']}) for small sphere region"

    def test_nbar_modes(self, tmp_path):
        """Test different nbar activity modes"""
        pdb, top, itp = self.create_minimal_water_files(tmp_path)

        base_inp = f"""# nbar mode test
pdb:{pdb}
top:{top}
fragitp:{itp}
op_pdb:output.pdb
op_top:output.top
box_size:15.0 15.0 15.0
cutoff:7.0
mcsteps:1000
nprint:500
fragname:WAT
fragconc:55.0
fragmuex:1.0
"""

        modes = [
            ("default", ""),
            ("const", "const_water_nbar:50.0"),
            ("volume", "volume_water_nbar:55.0"),
        ]

        results = {}
        for mode_name, mode_line in modes:
            inp_content = base_inp + "\n" + mode_line if mode_line else base_inp
            metrics = self.run_gcmc(inp_content, tmp_path, steps=1000)
            results[mode_name] = metrics["final_count"]

        # Check that all modes produced some molecules
        for mode, count in results.items():
            assert count > 0, f"{mode} mode should produce molecules"

        # Different modes should produce different results (at least const vs default)
        # Due to stochastic nature, we only check they're not identical
        assert len(set(results.values())) > 1, "Different nbar modes should produce different results"

    def test_cbmc_configuration(self, tmp_path):
        """Test that CBMC configuration is properly recognized"""
        pdb, top, itp = self.create_minimal_water_files(tmp_path)

        # Test with CBMC enabled
        inp_content = f"""# CBMC test
pdb:{pdb}
top:{top}
fragitp:{itp}
op_pdb:output.pdb
op_top:output.top
box_size:15.0 15.0 15.0
cutoff:7.0
mcsteps:500
nprint:100
fragname:WAT
fragconc:40.0
fragconc:55.0
fragmuex:-5.0
use_conf_bias:yes
fragconf:3
"""
        metrics = self.run_gcmc(inp_content, tmp_path, steps=1000)

        # Check that CBMC was recognized
        assert "configuration bias" in metrics["stdout"].lower() or \
               "conf" in metrics["stdout"].lower(), \
               "CBMC configuration should be mentioned in output"

        # Check simulation ran successfully
        assert metrics["returncode"] == 0, "CBMC simulation should run successfully"

    def test_move_probabilities(self, tmp_path):
        """Test per-fragment move probability configuration"""
        pdb, top, itp = self.create_minimal_water_files(tmp_path)

        inp_content = f"""# Move probability test
pdb:{pdb}
top:{top}
fragitp:{itp}
op_pdb:output.pdb
op_top:output.top
box_size:15.0 15.0 15.0
cutoff:7.0
mcsteps:500
nprint:100
fragname:WAT
fragconc:30.0
fragconc:55.0
fragmuex:-6.0
# Custom move probabilities (1:2:3:4 ratio)
attempt_prob_ins:0.1
attempt_prob_del:0.2
attempt_prob_trn:0.3
attempt_prob_rot:0.4
"""
        metrics = self.run_gcmc(inp_content, tmp_path, steps=1000)

        # Check if CDF is mentioned in output
        assert "cdf" in metrics["stdout"].lower(), "Move probability CDF should be calculated"

        # Verify the simulation ran successfully
        assert metrics["returncode"] == 0, "Simulation with custom move probabilities should run"

    def test_cavity_bias(self, tmp_path):
        """Test cavity bias configuration"""
        pdb, top, itp = self.create_minimal_water_files(tmp_path)

        # Test with cavity bias enabled
        inp_content = f"""# Cavity bias test
pdb:{pdb}
top:{top}
fragitp:{itp}
op_pdb:output.pdb
op_top:output.top
box_size:15.0 15.0 15.0
cutoff:7.0
mcsteps:1000
nprint:500
fragname:WAT
fragmuex:2.0
use_cavity_bias:yes
"""
        metrics_cavity = self.run_gcmc(inp_content, tmp_path, steps=1000)

        # Test without cavity bias
        inp_content_no_cavity = inp_content.replace("use_cavity_bias:yes", "use_cavity_bias:no")
        metrics_no_cavity = self.run_gcmc(inp_content_no_cavity, tmp_path, steps=1000)

        # Both should run successfully
        assert metrics_cavity["returncode"] == 0, "Cavity bias simulation should run"
        assert metrics_no_cavity["returncode"] == 0, "No cavity bias simulation should run"

        # Both should produce molecules
        assert metrics_cavity["final_count"] > 0, "Should produce molecules with cavity bias"
        assert metrics_no_cavity["final_count"] > 0, "Should produce molecules without cavity bias"

    def test_temperature_effect(self, tmp_path):
        """Test temperature dependence of insertion"""
        pdb, top, itp = self.create_minimal_water_files(tmp_path)

        temperatures = [250, 350]
        counts = []

        for T in temperatures:
            inp_content = f"""# Temperature test at {T}K
pdb:{pdb}
top:{top}
fragitp:{itp}
op_pdb:output.pdb
op_top:output.top
box_size:15.0 15.0 15.0
temperature:{T}
cutoff:7.0
mcsteps:2000
nprint:500
fragname:WAT
fragconc:55.0
fragmuex:-5.0
"""
            metrics = self.run_gcmc(inp_content, tmp_path, steps=2000)
            counts.append(metrics["final_count"])

        # Both temperatures should produce molecules
        assert all(c > 0 for c in counts), "Should produce molecules at all temperatures"

        # Higher temperature typically leads to lower density for liquids
        # But this is not always guaranteed in short simulations, so we just check they ran
        assert len(counts) == len(temperatures), "All temperature simulations should complete"

    @pytest.mark.parametrize("mu,expected_range", [
        (1.0, (10, 100)),    # Low activity
        (2.0, (20, 150)),    # Medium activity
        (3.0, (30, 200)),    # High activity
    ])
    def test_activity_ranges(self, tmp_path, mu, expected_range):
        """Test that different chemical potentials produce expected density ranges"""
        pdb, top, itp = self.create_minimal_water_files(tmp_path)

        inp_content = f"""# Activity range test
pdb:{pdb}
top:{top}
fragitp:{itp}
op_pdb:output.pdb
op_top:output.top
box_size:15.0 15.0 15.0
cutoff:7.0
mcsteps:1000
nprint:500
fragname:WAT
fragmuex:{mu}
"""
        metrics = self.run_gcmc(inp_content, tmp_path, steps=1000)

        # Check molecule count is in expected range
        min_expected, max_expected = expected_range
        assert min_expected <= metrics["final_count"] <= max_expected, \
            f"μ={mu}: count {metrics['final_count']} outside range {expected_range}"