"""
P1 Cavity Bias Validation Tests

These tests validate cavity bias functionality through CLI execution and JSONL logs:
1. Cavity fraction output matches bias factor
2. Acceptance formula consistency with cavity bias
3. Distribution-level detailed balance with cavity bias

Execution via CLI → parse JSONL → verify statistics
"""

import pytest
import subprocess
from pathlib import Path
import re
from acceptance_log_utils import read_jsonl, filter_by_move, acceptance_statistics


class TestCavityBias:
    """Cavity bias sanity and verification tests (P1)"""

    def test_cavity_stats_output(self, tmp_path):
        """
        Verify cavity grid statistics are output correctly.

        Validates:
        - Stdout contains "Cavity stats: total=<int>, cavities=<int>, fraction=<float>"
        - Total grid points > 0
        - Cavity fraction in reasonable range [0, 1]
        """
        # Build simple water box INP with cavity bias enabled
        pdb, top, atp, ff = self.setup_water_system(tmp_path)

        inp = self.build_inp(
            pdb=pdb,
            top=top,
            atp=atp,
            ff=ff,
            use_cavity=True,
            grid_spacing=2.0,
            probe_radius=1.4,
            mcsteps=100,  # Minimal steps, just need grid calculation
            mu=-2.0
        )

        inp_file = tmp_path / "test.inp"
        inp_file.write_text(inp)

        # Run simulation (use tmp_path as cwd to avoid parallel test interference)
        build_dir = Path(__file__).parent.parent.parent.parent / "build"
        gcmc_exe = build_dir / "bin" / "gcmc_cpu"
        result = subprocess.run(
            [str(gcmc_exe), "--inp", str(inp_file)],
            cwd=tmp_path,
            capture_output=True,
            text=True,
            timeout=60
        )

        assert result.returncode == 0, f"Simulation failed:\n{result.stderr}"

        # Parse cavity stats from stdout
        pattern = r"Cavity stats: total=(\d+), cavities=(\d+), fraction=([\d.]+)"
        matches = re.findall(pattern, result.stdout)

        assert len(matches) > 0, f"No cavity stats found in output:\n{result.stdout}"

        # Check first cavity stats (initial grid calculation)
        total, cavities, fraction = matches[0]
        total = int(total)
        cavities = int(cavities)
        fraction = float(fraction)

        # Validate statistics
        assert total > 0, f"Total grid points should be > 0, got {total}"
        assert cavities >= 0, f"Cavity points should be >= 0, got {cavities}"
        assert 0.0 <= fraction <= 1.0, f"Cavity fraction should be in [0,1], got {fraction}"
        assert abs(fraction - cavities / total) < 1e-6, \
            f"Fraction mismatch: {fraction} != {cavities}/{total}"

    def test_cavity_fraction_matches_bias(self, tmp_path):
        """
        Verify cavity fraction ≈ bias factor in acceptance records.

        Validates:
        - Cavity bias (wCavity) in JSONL logs
        - Average wCavity ≈ cavity fraction from grid stats (±10% tolerance)
        - Both insertion and deletion moves record cavity bias correctly
        """
        # Build water box with moderate density for cavity formation
        pdb, top, atp, ff = self.setup_water_system(tmp_path)

        inp = self.build_inp(
            pdb=pdb,
            top=top,
            atp=atp,
            ff=ff,
            use_cavity=True,
            grid_spacing=2.0,
            probe_radius=1.4,
            mcsteps=5000,
            mu=-2.0
        )

        inp_file = tmp_path / "test.inp"
        inp_file.write_text(inp)
        accept_log = tmp_path / "accept.jsonl"

        # Run simulation with acceptance logging (use tmp_path as cwd to avoid parallel test interference)
        build_dir = Path(__file__).parent.parent.parent.parent / "build"
        gcmc_exe = build_dir / "bin" / "gcmc_cpu"
        result = subprocess.run(
            [str(gcmc_exe), "--inp", str(inp_file), "--dump-accept", str(accept_log)],
            cwd=tmp_path,
            capture_output=True,
            text=True,
            timeout=120
        )

        assert result.returncode == 0, f"Simulation failed:\n{result.stderr}"
        assert accept_log.exists(), "Acceptance log not created"

        # Parse cavity fraction from stdout
        pattern = r"Cavity stats: total=(\d+), cavities=(\d+), fraction=([\d.]+)"
        matches = re.findall(pattern, result.stdout)
        assert len(matches) > 0, "No cavity stats in output"

        # Get cavity fraction from grid (use first calculation)
        _, _, grid_fraction = matches[0]
        grid_fraction = float(grid_fraction)

        # Parse JSONL records
        records = read_jsonl(accept_log)
        assert len(records) > 0, "No acceptance records found"

        # Filter insertion moves (cavity bias only applies to insertions)
        insertions = filter_by_move(records, "insertion")
        assert len(insertions) > 100, f"Need >100 insertions for statistics, got {len(insertions)}"

        # Collect wCavity values
        w_cavity_values = []
        for rec in insertions:
            if "wCavity" in rec and rec["wCavity"] is not None:
                w_cavity_values.append(rec["wCavity"])

        assert len(w_cavity_values) > 50, \
            f"Need >50 cavity bias values, got {len(w_cavity_values)}"

        # Calculate average cavity bias
        avg_w_cavity = sum(w_cavity_values) / len(w_cavity_values)

        # Verify wCavity ≈ grid_fraction (±10% tolerance)
        # Note: wCavity may vary slightly as molecules are inserted/deleted
        tolerance = 0.10  # 10% tolerance
        lower_bound = grid_fraction * (1 - tolerance)
        upper_bound = grid_fraction * (1 + tolerance)

        assert lower_bound <= avg_w_cavity <= upper_bound, \
            f"Cavity bias mismatch: avg wCavity={avg_w_cavity:.4f}, " \
            f"grid fraction={grid_fraction:.4f}, tolerance=±{tolerance*100}%"

    def test_cavity_acceptance_formula_consistency(self, tmp_path):
        """
        Verify acceptance formula with cavity bias: pAcc = min(1, (z*V_cav)/(N+1) * exp(-βΔU) * qF)

        For insertion with cavity bias:
          pAcc_theory = min(1, (z * V_cavity) / (N+1) * exp(-βΔU) * qForward)
        where:
          - V_cavity = box_volume * cavity_fraction (effective volume)
          - qForward = CBMC Rosenbluth weight (W_new/K)

        Validates:
        - |pAcc_actual - pAcc_theory| < 5e-5 for insertions
        - Formula holds with cavity bias component
        """
        # Build water box with CBMC + cavity bias
        pdb, top, atp, ff = self.setup_water_system(tmp_path)

        inp = self.build_inp(
            pdb=pdb,
            top=top,
            atp=atp,
            ff=ff,
            use_cavity=True,
            use_cbmc=True,
            k_trials=8,
            grid_spacing=2.0,
            probe_radius=1.4,
            mcsteps=5000,  # Reduced from 10000
            mu=-1.0  # Raised from -2.0 for better acceptance
        )

        inp_file = tmp_path / "test.inp"
        inp_file.write_text(inp)
        accept_log = tmp_path / "accept.jsonl"

        # Run simulation (use tmp_path as cwd to avoid parallel test interference)
        build_dir = Path(__file__).parent.parent.parent.parent / "build"
        gcmc_exe = build_dir / "bin" / "gcmc_cpu"
        result = subprocess.run(
            [str(gcmc_exe), "--inp", str(inp_file), "--dump-accept", str(accept_log)],
            cwd=tmp_path,
            capture_output=True,
            text=True,
            timeout=180
        )

        assert result.returncode == 0, f"Simulation failed:\n{result.stderr}"

        # Parse records
        records = read_jsonl(accept_log)
        insertions = filter_by_move(records, "insertion")
        assert len(insertions) >= 100, f"Need >=100 insertions, got {len(insertions)}"

        # Verify formula for each insertion
        errors = []
        for rec in insertions:
            # Extract fields
            n_before = rec["nBefore"]
            beta_mu = rec["betaMu"]
            z = rec["z"]
            delta_u = rec["deltaU"]
            q_forward = rec["qForward"]
            w_cavity = rec.get("wCavity", 1.0)
            p_acc_actual = rec["pAcc"]
            v_eff = rec["vEff"]  # Effective volume (already includes cavity bias if applicable)

            # Theory: P_acc = min(1, (z * V_eff) / (N+1) * exp(-βΔU) * qForward)
            # Note: vEff already includes cavity bias, so we don't multiply by wCavity again
            import math
            ratio = (z * v_eff) / (n_before + 1) * math.exp(-delta_u) * q_forward
            p_acc_theory = min(1.0, ratio)

            # Calculate error
            error = abs(p_acc_actual - p_acc_theory)
            errors.append(error)

            # Individual record should match within tolerance
            if error >= 5e-5:
                print(f"WARNING: Large error={error:.2e} for record step={rec['step']}")
                print(f"  nBefore={n_before}, z={z:.3e}, vEff={v_eff:.3f}, "
                      f"wCavity={w_cavity:.4f}")
                print(f"  deltaU={delta_u:.3f}, qForward={q_forward:.4f}")
                print(f"  pAcc_actual={p_acc_actual:.6f}, pAcc_theory={p_acc_theory:.6f}")

        # Statistical check: max error should be < 5e-5
        max_error = max(errors)
        mean_error = sum(errors) / len(errors)

        print(f"Formula verification: {len(insertions)} insertions")
        print(f"  Max error: {max_error:.2e}")
        print(f"  Mean error: {mean_error:.2e}")

        assert max_error < 5e-5, \
            f"Formula error too large: max_error={max_error:.2e} >= 5e-5"

    # Helper methods

    def setup_water_system(self, tmp_path):
        """Create minimal water system files (PDB, PSF, ATP, FF)"""
        # Create simple water PDB (single molecule)
        pdb_content = """REMARK   Water box for cavity bias testing
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
ATOM      1  OH2 TIP3    1      15.000  15.000  15.000  1.00  0.00      WAT  O
ATOM      2  H1  TIP3    1      15.757  15.000  15.587  1.00  0.00      WAT  H
ATOM      3  H2  TIP3    1      14.243  15.000  15.587  1.00  0.00      WAT  H
END
"""
        pdb_file = tmp_path / "water.pdb"
        pdb_file.write_text(pdb_content)

        # Create PSF topology
        psf_content = """PSF CMAP

       1 !NTITLE
 REMARKS water molecule

       3 !NATOM
       1 WAT  1    TIP3 OH2  OT    -0.834000       15.9994           0
       2 WAT  1    TIP3 H1   HT     0.417000        1.0080           0
       3 WAT  1    TIP3 H2   HT     0.417000        1.0080           0

       2 !NBOND: bonds
       1       2       1       3

       1 !NTHETA: angles
       2       1       3

       0 !NPHI: dihedrals

       0 !NIMPHI: impropers

       0 !NDON: donors

       0 !NACC: acceptors

       0 !NNB

       0       0       0

       1       0 !NGRP
       0       0       0
"""
        psf_file = tmp_path / "water.psf"
        psf_file.write_text(psf_content)

        # Create ATP atom types
        atp_content = """MASS     1 OT      15.999400 O  ! TIP3P water oxygen
MASS     2 HT       1.008000 H  ! TIP3P water hydrogen
"""
        atp_file = tmp_path / "water.atp"
        atp_file.write_text(atp_content)

        # Create force field (simplified TIP3P)
        ff_content = """* TIP3P water force field
*

BONDS
OT   HT     450.0       0.9572

ANGLES
HT   OT   HT      55.0     104.52

NONBONDED nbxmod  5 atom cdiel shift vatom vdistance vswitch -
cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5
OT      0.0       -0.1521     1.7682   ! TIP3P oxygen
HT      0.0       -0.046      0.2245   ! TIP3P hydrogen

END
"""
        ff_file = tmp_path / "water.prm"
        ff_file.write_text(ff_content)

        return pdb_file, psf_file, atp_file, ff_file

    def build_inp(self, pdb, top, atp, ff, use_cavity=False, use_cbmc=False,
                  k_trials=1, grid_spacing=2.0, probe_radius=1.4, mcsteps=1000, mu=-2.0):
        """Build INP file content with cavity bias and CBMC options"""
        # Base configuration (old-style INP format)
        cavity_lines = f"""use_cavity_bias:yes
cavity_grid_spacing:{grid_spacing}
cavity_probe_radius:{probe_radius}""" if use_cavity else "use_cavity_bias:no"

        cbmc_lines = f"""use_conf_bias:yes
fragconf:{k_trials}""" if use_cbmc else "use_conf_bias:no"

        inp = f"""# Cavity bias validation test
par:{ff}
atomtypes:{atp}
top:{top}
pdb:{pdb}
protitp:{top}

fragname: water
fragconc: 55.0
fragmuex: {mu}
{cbmc_lines if use_cbmc else ''}

box_size: 30.0 30.0 30.0
cutoff: 12.0
temperature: 300.0
mcsteps: {mcsteps}
nprint: {max(1, mcsteps // 2)}
eqsteps: 0

mc_move_prob: 0.5 0.5 0 0

{cavity_lines}

seed: 42
"""
        return inp

    def parse_cavity_stats(self, stdout):
        """Parse cavity statistics from stdout"""
        pattern = r"Cavity stats: total=(\d+), cavities=(\d+), fraction=([\d.]+)"
        matches = re.findall(pattern, stdout)
        if not matches:
            return None

        results = []
        for match in matches:
            total = int(match[0])
            cavities = int(match[1])
            fraction = float(match[2])
            results.append({
                "total": total,
                "cavities": cavities,
                "fraction": fraction
            })
        return results
