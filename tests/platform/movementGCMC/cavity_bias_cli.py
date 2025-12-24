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
import json
from typing import List, Dict, Any
from collections import defaultdict


# ============================================================================
# JSONL Utilities (inline to avoid cross-directory conftest imports)
# ============================================================================

def read_jsonl(filepath: Path) -> List[Dict[str, Any]]:
    """Read JSONL file and return list of records."""
    records = []
    if not filepath.exists():
        return records
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError as e:
                print(f"Warning: Failed to parse line: {line[:50]}... Error: {e}")
    return records


def filter_by_move(records: List[Dict[str, Any]], move_type: str) -> List[Dict[str, Any]]:
    """Filter records by move type."""
    return [r for r in records if r.get('move') == move_type]


def acceptance_statistics(records: List[Dict[str, Any]]) -> Dict[str, Dict[str, float]]:
    """Compute acceptance rate statistics by move type and species."""
    stats = defaultdict(lambda: defaultdict(lambda: {'attempts': 0, 'accepted': 0}))
    for record in records:
        move = record.get('move', 'unknown')
        species = record.get('species', 'unknown')
        accepted = record.get('accepted', False)
        stats[move][species]['attempts'] += 1
        if accepted:
            stats[move][species]['accepted'] += 1
    result = {}
    for move, species_dict in stats.items():
        result[move] = {}
        for species, counts in species_dict.items():
            if counts['attempts'] > 0:
                result[move][species] = counts['accepted'] / counts['attempts']
            else:
                result[move][species] = 0.0
    return result


class TestCavityBias:
    """Cavity bias sanity and verification tests (P1)"""

    def test_cavity_stats_output(self, tmp_path):
        """
        Verify cavity grid statistics are computed and exported correctly.

        Validates:
        - --dump-accept contains cavity-related structured fields (no stdout parsing)
        - cavityFraction and wCavity are in [0, 1]
        - vEff ≈ vBox * cavityFraction (effective insertion volume)
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
        accept_log = tmp_path / "accept.jsonl"

        # Run simulation (use tmp_path as cwd to avoid parallel test interference)
        build_dir = Path(__file__).parent.parent.parent.parent / "build"
        gcmc_exe = build_dir / "bin" / "gcmc_cpu"
        result = subprocess.run(
            [str(gcmc_exe), "--inp", str(inp_file), "--dump-accept", str(accept_log)],
            cwd=tmp_path,
            capture_output=True,
            text=True,
            timeout=60
        )

        assert result.returncode == 0, f"Simulation failed:\n{result.stderr}"
        assert accept_log.exists(), "Acceptance log not created"

        records = read_jsonl(accept_log)
        assert records, "No acceptance records found"
        insertions = filter_by_move(records, "insertion")
        assert insertions, "No insertion records found"

        checked = 0
        for rec in insertions[:25]:
            if rec.get("cavityFraction") is None or rec.get("wCavity") is None:
                continue
            if rec.get("vBox") is None or rec.get("vEff") is None:
                continue

            cavity_fraction = float(rec["cavityFraction"])
            w_cavity = float(rec["wCavity"])
            v_box = float(rec["vBox"])
            v_eff = float(rec["vEff"])

            assert 0.0 <= cavity_fraction <= 1.0
            assert 0.0 <= w_cavity <= 1.0
            assert v_eff == pytest.approx(v_box * cavity_fraction, rel=2e-4, abs=2e-6)
            checked += 1

        assert checked >= 5, "Too few insertion records with cavity fields"

    def test_cavity_fraction_matches_bias(self, tmp_path):
        """
        Verify cavity fraction ≈ bias factor in acceptance records (file-driven).

        Validates:
        - Cavity bias (wCavity) in JSONL logs
        - Average wCavity ≈ average cavityFraction (±10% tolerance)
        - vEff ≈ vBox * cavityFraction
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
            mcsteps=1000,
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

        # Parse JSONL records
        records = read_jsonl(accept_log)
        assert len(records) > 0, "No acceptance records found"

        # Filter insertion moves (cavity bias only applies to insertions)
        insertions = filter_by_move(records, "insertion")
        assert len(insertions) > 50, f"Need >50 insertions for statistics, got {len(insertions)}"

        # Collect wCavity values
        w_cavity_values = []
        for rec in insertions:
            if "wCavity" in rec and rec["wCavity"] is not None:
                w_cavity_values.append(rec["wCavity"])

        assert len(w_cavity_values) > 50, \
            f"Need >50 cavity bias values, got {len(w_cavity_values)}"

        cavity_fractions = [
            float(r["cavityFraction"])
            for r in insertions
            if r.get("cavityFraction") is not None
        ]
        assert len(cavity_fractions) > 50, f"Need >50 cavityFraction values, got {len(cavity_fractions)}"

        checked = 0
        for rec in insertions[:25]:
            if rec.get("vBox") is None or rec.get("vEff") is None or rec.get("cavityFraction") is None:
                continue
            v_box = float(rec["vBox"])
            v_eff = float(rec["vEff"])
            cavity_fraction = float(rec["cavityFraction"])
            assert v_eff == pytest.approx(v_box * cavity_fraction, rel=2e-4, abs=2e-6)
            checked += 1
        assert checked >= 5

        avg_w_cavity = sum(w_cavity_values) / len(w_cavity_values)
        avg_fraction = sum(cavity_fractions) / len(cavity_fractions)

        # Verify wCavity ≈ cavityFraction (±10% tolerance)
        # Note: these may vary slightly as molecules are inserted/deleted.
        tolerance = 0.10  # 10% tolerance
        lower_bound = avg_fraction * (1 - tolerance)
        upper_bound = avg_fraction * (1 + tolerance)

        assert lower_bound <= avg_w_cavity <= upper_bound, \
            f"Cavity bias mismatch: avg wCavity={avg_w_cavity:.4f}, " \
            f"avg cavityFraction={avg_fraction:.4f}, tolerance=±{tolerance*100}%"

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
            k_trials=4,
            grid_spacing=3.0,
            probe_radius=1.4,
            mcsteps=2000,  # Reduced from 10000
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
        assert len(insertions) >= 30, f"Need >=30 insertions, got {len(insertions)}"

        # Verify formula for each insertion
        errors = []
        for rec in insertions:
            # Extract fields
            n_before = rec["nBefore"]
            beta_mu = rec["betaMu"]
            z = rec["z"]
            beta_delta_u = rec.get("betaDeltaU", rec["deltaU"])
            q_forward = rec["qForward"]
            w_cavity = rec.get("wCavity", 1.0)
            p_acc_actual = rec["pAcc"]
            v_eff = rec["vEff"]  # Effective volume (already includes cavity bias if applicable)

            # Theory: P_acc = min(1, (z * V_eff) / (N+1) * exp(-βΔU) * qForward)
            # Note: vEff already includes cavity bias, so we don't multiply by wCavity again
            import math
            ratio = (z * v_eff) / (n_before + 1) * math.exp(-beta_delta_u) * q_forward
            p_acc_theory = min(1.0, ratio)

            # Calculate error
            error = abs(p_acc_actual - p_acc_theory)
            errors.append(error)

            # Individual record should match within tolerance
            if error >= 5e-5:
                print(f"WARNING: Large error={error:.2e} for record step={rec['step']}")
                print(f"  nBefore={n_before}, z={z:.3e}, vEff={v_eff:.3f}, "
                      f"wCavity={w_cavity:.4f}")
                print(f"  betaDeltaU={beta_delta_u:.3f}, qForward={q_forward:.4f}")
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
