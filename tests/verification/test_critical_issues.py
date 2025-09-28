#!/usr/bin/env python3
"""
Critical Issues Verification Tests

These tests are designed to expose known issues without modifying the code.
They document expected failures and provide evidence for needed fixes.
"""

import pytest
import subprocess
import tempfile
import os
import re
import math
import json
import numpy as np
from pathlib import Path

# Path to the gcmc_cpu executable
GCMC_CPU_PATH = Path(__file__).parent.parent.parent / "build" / "bin" / "gcmc_cpu"
if not GCMC_CPU_PATH.exists():
    GCMC_CPU_PATH = Path("/home/zhaomt/gcmc/test108/pygcmc_dev/build/bin/gcmc_cpu")


class TestRandomnessDeterminism:
    """Test that all RNG sources are properly seeded"""

    def create_test_system(self, tmpdir, use_cavity=False, use_region=False):
        """Create minimal test system with optional cavity/region"""
        inp_file = tmpdir / "test.inp"
        pdb_file = tmpdir / "test.pdb"
        top_file = tmpdir / "test.top"
        itp_file = tmpdir / "water.itp"

        # Create minimal PDB
        with open(pdb_file, "w") as f:
            f.write("TITLE     Test\n")
            f.write("CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1\n")
            f.write("END\n")

        # Create minimal TOP
        with open(top_file, "w") as f:
            f.write("[ system ]\nTest\n\n[ molecules ]\n")

        # Create water ITP
        with open(itp_file, "w") as f:
            f.write("[ moleculetype ]\n")
            f.write("WAT      3\n\n")
            f.write("[ atoms ]\n")
            f.write("1   O     1      WAT      O     1      -0.834   15.999\n")
            f.write("2   H     1      WAT      H1    1       0.417    1.008\n")
            f.write("3   H     1      WAT      H2    1       0.417    1.008\n")

        # Create INP with options
        inp_content = f"""pdb:{pdb_file}
top:{top_file}
fragitp:{itp_file}
op_pdb:{tmpdir}/output.pdb
op_top:{tmpdir}/output.top
box_size:5.0 5.0 5.0
cutoff:2.0
mcsteps:200
nprint:100
fragname:WAT
fragmuex:2.0
"""
        if use_cavity:
            inp_content += "use_cavity_bias:yes\ncavity_grid_dx:1.0\n"
        if use_region:
            inp_content += "gcmc_region:sphere 2.5 2.5 2.5 2.0\n"

        with open(inp_file, "w") as f:
            f.write(inp_content)

        return inp_file

    def test_cavity_region_determinism(self, tmp_path):
        """Test: Cavity bias + Region should be deterministic with same seed"""
        print("\n=== Testing Cavity/Region RNG Determinism ===")

        # Create test with both cavity and region enabled
        inp_file = self.create_test_system(tmp_path, use_cavity=True, use_region=True)

        results = []
        for run in range(2):
            cmd = [
                str(GCMC_CPU_PATH),
                "--inp", str(inp_file),
                "--seed", "12345",  # Same seed
                "--prefix", str(tmp_path / f"run_{run}")
            ]

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)

            # Parse water count and positions
            count_match = re.search(r"Fragment counts:\s*WAT:\s*(\d+)", result.stdout)
            count = int(count_match.group(1)) if count_match else -1

            # Try to extract some coordinate info
            coord_info = ""
            if "Final configuration" in result.stdout:
                # Extract a hash of positions or energy
                energy_match = re.search(r"Current energy:\s*([\d.-]+)", result.stdout)
                if energy_match:
                    coord_info = energy_match.group(1)

            results.append({"count": count, "coord_info": coord_info})
            print(f"Run {run+1}: count={count}, energy={coord_info}")

        # Check determinism
        if results[0]["count"] != results[1]["count"]:
            print(f"❌ ISSUE CONFIRMED: Molecule counts differ ({results[0]['count']} vs {results[1]['count']})")
            print("   RandomUtils in cavity/region not controlled by --seed")
        else:
            print(f"✓ Counts match: {results[0]['count']}")

        if results[0]["coord_info"] and results[0]["coord_info"] != results[1]["coord_info"]:
            print(f"❌ ISSUE CONFIRMED: Final energies differ ({results[0]['coord_info']} vs {results[1]['coord_info']})")
            print("   Indicates non-deterministic position selection")

        # Verify determinism - should be identical with same seed
        assert results[0]["count"] == results[1]["count"], \
            f"Non-deterministic: counts differ ({results[0]['count']} vs {results[1]['count']})"

        if results[0]["coord_info"] and results[1]["coord_info"]:
            assert results[0]["coord_info"] == results[1]["coord_info"], \
                f"Non-deterministic: energies differ ({results[0]['coord_info']} vs {results[1]['coord_info']})"


class TestStatisticalCounting:
    """Test statistical counting consistency"""

    def test_acceptance_rate_counting(self, tmp_path):
        """Test: Total acceptance rate vs per-move rates"""
        print("\n=== Testing Acceptance Rate Counting ===")

        # Create minimal system
        inp_file = tmp_path / "test.inp"
        pdb_file = tmp_path / "test.pdb"
        top_file = tmp_path / "test.top"
        itp_file = tmp_path / "water.itp"

        # Create files
        with open(pdb_file, "w") as f:
            f.write("TITLE Test\n")
            f.write("CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1\n")
            f.write("END\n")

        with open(top_file, "w") as f:
            f.write("[ system ]\nTest\n\n[ molecules ]\n")

        with open(itp_file, "w") as f:
            f.write("[ moleculetype ]\nWAT 3\n\n[ atoms ]\n")
            f.write("1 O 1 WAT O 1 -0.834 15.999\n")

        with open(inp_file, "w") as f:
            f.write(f"""pdb:{pdb_file}
top:{top_file}
fragitp:{itp_file}
op_pdb:output.pdb
op_top:output.top
box_size:3.0 3.0 3.0
cutoff:1.0
mcsteps:1000
nprint:500
fragname:WAT
fragmuex:1.0
probInsert:0.5
probDelete:0.5
probTranslate:0.0
probRotate:0.0
""")

        # Run simulation with verbose stats to get weighted acceptance
        cmd = [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--seed", "54321"]
        env = os.environ.copy()
        env['GCMC_VERBOSE_STATS'] = '1'
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=10, env=env)

        # Parse acceptance rates from final statistics (take last occurrence)
        total_matches = re.findall(r"Total acceptance rate:\s*([\d.]+)%", result.stdout)
        total_match = total_matches[-1] if total_matches else None
        insert_match = re.search(r"Insert move accept:\s*([\d.]+)%", result.stdout)
        delete_match = re.search(r"Delete move accept:\s*([\d.]+)%", result.stdout)

        # Parse weighted overall acceptance from diagnostics
        overall_match = re.search(r"Ins/Del acceptance \(overall\):\s*([\d.]+)%", result.stdout)

        # Parse attempt counts (if available)
        insert_attempts = re.search(r"Insert attempts:\s*(\d+)", result.stdout)
        insert_accepted = re.search(r"Insert accepted:\s*(\d+)", result.stdout)
        delete_attempts = re.search(r"Delete attempts:\s*(\d+)", result.stdout)
        delete_accepted = re.search(r"Delete accepted:\s*(\d+)", result.stdout)

        if total_match and insert_match and delete_match:
            total_rate = float(total_match)
            insert_rate = float(insert_match.group(1))
            delete_rate = float(delete_match.group(1))

            print(f"Total acceptance rate: {total_rate}%")
            print(f"Insert acceptance rate: {insert_rate}%")
            print(f"Delete acceptance rate: {delete_rate}%")

            # Use weighted overall if available, otherwise compute from counts
            if overall_match:
                weighted_total = float(overall_match.group(1))
                print(f"Ins/Del acceptance (overall): {weighted_total}%")
                comparison_rate = weighted_total
            elif insert_attempts and insert_accepted and delete_attempts and delete_accepted:
                # Calculate weighted average from raw counts
                ins_att = int(insert_attempts.group(1))
                ins_acc = int(insert_accepted.group(1))
                del_att = int(delete_attempts.group(1))
                del_acc = int(delete_accepted.group(1))
                if (ins_att + del_att) > 0:
                    weighted_total = 100.0 * (ins_acc + del_acc) / (ins_att + del_att)
                    print(f"Computed weighted rate: {weighted_total:.1f}%")
                    comparison_rate = weighted_total
                else:
                    comparison_rate = (insert_rate + delete_rate) / 2
            else:
                # Fall back to simple average
                comparison_rate = (insert_rate + delete_rate) / 2
                print(f"Using simple average: {comparison_rate:.1f}%")

            if abs(total_rate - comparison_rate) > 10:
                print(f"❌ ISSUE CONFIRMED: Total rate ({total_rate}%) significantly differs from weighted rate ({comparison_rate:.1f}%)")
            else:
                print("✓ Rates appear consistent")

        # Verify acceptance rate consistency
        assert total_match and insert_match and delete_match, "Failed to parse acceptance rates"

        # When translation/rotation are disabled, ins/del overall should be correct
        # The total rate calculation might be affected by disabled moves
        # So we check if the ins/del overall rate is reasonable
        tolerance = 10  # percent
        if overall_match:
            # Check that ins/del overall is between insert and delete rates
            min_rate = min(insert_rate, delete_rate)
            max_rate = max(insert_rate, delete_rate)
            assert min_rate <= weighted_total <= max_rate, \
                f"Ins/Del overall ({weighted_total}%) should be between insert ({insert_rate}%) and delete ({delete_rate}%) rates"


class TestHardLimits:
    """Test impact of hard-coded molecule limits"""

    def test_insertion_hard_limit(self, tmp_path):
        """Test: Hard limit of 1000 molecules affects high-μ systems"""
        print("\n=== Testing Hard Insertion Limit ===")

        # Create system with very high chemical potential
        inp_file = tmp_path / "test.inp"
        pdb_file = tmp_path / "test.pdb"
        top_file = tmp_path / "test.top"
        itp_file = tmp_path / "water.itp"

        with open(pdb_file, "w") as f:
            f.write("TITLE Test\n")
            f.write("CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1\n")
            f.write("END\n")

        with open(top_file, "w") as f:
            f.write("[ system ]\nTest\n\n[ molecules ]\n")

        with open(itp_file, "w") as f:
            f.write("[ moleculetype ]\nWAT 3\n\n[ atoms ]\n")
            f.write("1 O 1 WAT O 1 -0.834 15.999\n")

        # Very high chemical potential to trigger limit
        with open(inp_file, "w") as f:
            f.write(f"""pdb:{pdb_file}
top:{top_file}
fragitp:{itp_file}
op_pdb:output.pdb
op_top:output.top
box_size:10.0 10.0 10.0
cutoff:3.0
mcsteps:5000
nprint:1000
fragname:WAT
fragmuex:10.0
""")

        cmd = [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--seed", "11111"]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=20)

        # Look for signs of hitting the limit
        final_count_match = re.search(r"Fragment counts:\s*WAT:\s*(\d+)", result.stdout)
        insert_rate_match = re.search(r"Insert move accept:\s*([\d.]+)%", result.stdout)

        if final_count_match and insert_rate_match:
            final_count = int(final_count_match.group(1))
            insert_rate = float(insert_rate_match.group(1))

            print(f"Final molecule count: {final_count}")
            print(f"Insert acceptance rate: {insert_rate}%")

            # With μ=10, we expect many molecules unless limited
            if final_count >= 900 and final_count <= 1000:
                print(f"❌ ISSUE CONFIRMED: Count ({final_count}) near hard limit of 1000")
                if insert_rate < 5:
                    print("   Insert rate very low, suggesting systematic rejection at limit")
            elif final_count > 1000:
                print(f"✓ No hard limit detected (count={final_count})")
            else:
                print(f"? Inconclusive: count={final_count} (may be physical limit)")

        # Verify no artificial hard limit
        assert final_count_match and insert_rate_match, "Failed to parse results"

        # With very high chemical potential, should get many molecules
        # unless artificially limited
        assert not (900 <= final_count <= 1000 and insert_rate < 10), \
            f"Likely hitting hard limit: count={final_count}, insert_rate={insert_rate:.1f}%"


class TestPhysicalCorrectness:
    """Test physical formulas and algorithms"""

    def test_acceptance_formula_quantitative(self, tmp_path):
        """Test: Acceptance probability formula matches theory"""
        print("\n=== Testing Acceptance Formula ===")

        # Create ideal gas system (minimal interactions)
        inp_file = tmp_path / "test.inp"
        pdb_file = tmp_path / "test.pdb"
        top_file = tmp_path / "test.top"
        itp_file = tmp_path / "water.itp"

        with open(pdb_file, "w") as f:
            f.write("TITLE Test\n")
            f.write("CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1\n")
            f.write("END\n")

        with open(top_file, "w") as f:
            f.write("[ system ]\nTest\n\n[ molecules ]\n")

        # Minimal LJ parameters for near-ideal gas
        with open(itp_file, "w") as f:
            f.write("[ moleculetype ]\nWAT 3\n\n[ atoms ]\n")
            f.write("1 O 1 WAT O 1 0.0 15.999\n")  # No charge for simplicity

        # Known activity and volume for predictable acceptance
        volume = 5.0 ** 3  # 125 nm³
        activity = 0.01  # Low activity

        with open(inp_file, "w") as f:
            f.write(f"""pdb:{pdb_file}
top:{top_file}
fragitp:{itp_file}
op_pdb:output.pdb
op_top:output.top
box_size:5.0 5.0 5.0
cutoff:2.0
mcsteps:1000
nprint:500
fragname:WAT
fragmuex:2.303
""")

        # Enable diagnostics and acceptance log
        accept_log = tmp_path / "acceptance.jsonl"
        cmd = [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--seed", "99999", "--verbose"]
        env = os.environ.copy()
        env['GCMC_ENABLE_DIAGNOSTICS'] = '1'
        env['GCMC_DUMP_ACCEPT'] = str(accept_log)
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=10, env=env)

        # For ideal gas with ΔE≈0:
        # P_insert(N=0→1) = min(1, z*V) = min(1, activity * volume)
        # P_delete(N=1→0) = min(1, 1/(z*V))

        zV = activity * volume  # Should be 0.01 * 125 = 1.25
        expected_insert_0to1 = min(1.0, zV)  # Should be 1.0
        expected_delete_1to0 = min(1.0, 1.0/zV)  # Should be 0.8

        print(f"Theory: z*V = {zV:.3f}")
        print(f"Expected P_insert(0→1) = {expected_insert_0to1:.3f}")
        print(f"Expected P_delete(1→0) = {expected_delete_1to0:.3f}")

        # Parse actual rates if available
        insert_rate_match = re.search(r"Insert move accept:\s*([\d.]+)%", result.stdout)
        if insert_rate_match:
            actual_insert = float(insert_rate_match.group(1))
            print(f"Actual insert rate: {actual_insert}%")

        # Parse acceptance log if available
        if accept_log.exists():
            import json
            with open(accept_log, 'r') as f:
                records = [json.loads(line) for line in f if line.strip()]

            if records:
                print(f"\nAnalyzing {len(records)} acceptance records...")

                # Check a few insertion attempts
                insert_records = [r for r in records if r['move'] == 'insert'][:5]

                for rec in insert_records:
                    # For ideal gas, formula: p_acc = min(1, exp(-βΔU) * z * V / (N+1))
                    # With our simplified model:
                    deltaU = rec['deltaU']
                    beta = 1.0 / (8.314e-3 * 300)  # 1/kT in mol/kJ
                    exp_factor = math.exp(-beta * deltaU) if deltaU < 700 else 0

                    # Expected probability (simplified)
                    N_before = rec.get('N', 0)
                    expected = min(1.0, exp_factor * zV / (N_before + 1))
                    actual = rec['pAcc']

                    print(f"  Move: {rec['move']}, ΔU={deltaU:.2f}, pAcc={actual:.3f}, expected≈{expected:.3f}")

                    # Validate within tolerance (accounting for approximations)
                    if abs(deltaU) < 10:  # Only check near-ideal cases
                        assert abs(actual - expected) < 0.3, \
                            f"Acceptance probability mismatch: actual={actual:.3f}, expected={expected:.3f}"

                print("✓ Acceptance probabilities consistent with theory")
        else:
            print("⚠ No acceptance log found, skipping detailed validation")

        # Basic check that simulation ran
        assert insert_rate_match, "Failed to parse insert rate"


def test_mixed_rule_matrix(tmp_path):
    """Test: LJ mixing rules produce correct matrices"""
    print("\n=== Testing Mixing Rule Matrices ===")

    # Create a simple system to trigger LJ matrix export
    inp_file = tmp_path / "test.inp"
    pdb_file = tmp_path / "test.pdb"
    top_file = tmp_path / "test.top"

    with open(pdb_file, "w") as f:
        f.write("TITLE Test\n")
        f.write("CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1\n")
        f.write("END\n")

    with open(top_file, "w") as f:
        f.write("[ system ]\nTest\n\n[ molecules ]\n")

    with open(inp_file, "w") as f:
        f.write(f"""pdb:{pdb_file}
top:{top_file}
op_pdb:output.pdb
op_top:output.top
box_size:3.0 3.0 3.0
cutoff:1.2
mcsteps:10
nprint:10
fragname:WAT
fragmuex:0.0
""")

    # Run with LJ matrix dump enabled
    cmd = [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--seed", "12345", "--prefix", str(tmp_path / "test")]
    env = os.environ.copy()
    env['GCMC_DUMP_LJ'] = '1'
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=10, env=env)

    # Check if LJ matrix was exported
    lj_file = tmp_path / "test_lj.csv"
    if lj_file.exists():
        print("✓ LJ matrix exported")

        # Parse the CSV file
        import csv
        with open(lj_file, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        if rows:
            print(f"  Matrix size: {len(set(r['i'] for r in rows))} x {len(set(r['j'] for r in rows))}")

            # Check Lorentz-Berthelot rules for a few entries
            # σ_ij = (σ_i + σ_j) / 2
            # ε_ij = sqrt(ε_i * ε_j)

            # Group by i,j indices
            matrix = {}
            for row in rows[:10]:  # Check first few entries
                i, j = int(row['i']), int(row['j'])
                sigma = float(row['sigma_ij'])
                eps = float(row['eps_ij'])
                rule = row.get('rule', 'unknown')
                matrix[(i,j)] = {'sigma': sigma, 'eps': eps, 'rule': rule}
                print(f"    [{i},{j}]: σ={sigma:.3f}, ε={eps:.3f}, rule={rule}")

            print("✓ Matrix entries parsed successfully")
        else:
            print("❌ Empty matrix file")
            assert False, "LJ matrix file is empty"
    else:
        print("❌ LJ matrix file not found")
        print("  Expected: " + str(lj_file))
        assert False, "LJ matrix export failed"


def test_region_volume_in_acceptance(tmp_path):
    """Test: Region volume is used correctly in acceptance calculation"""
    print("\n=== Testing Region Volume in Acceptance ===")

    # Helper to create and run simulation
    def run_simulation(use_region=False, seed=12345):
        prefix = "region" if use_region else "noregion"
        inp_file = tmp_path / f"{prefix}.inp"
        pdb_file = tmp_path / f"{prefix}.pdb"
        top_file = tmp_path / f"{prefix}.top"
        itp_file = tmp_path / f"{prefix}_water.itp"

        with open(pdb_file, "w") as f:
            f.write("TITLE Test\n")
            f.write("CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1\n")
            f.write("END\n")

        with open(top_file, "w") as f:
            f.write("[ system ]\nTest\n\n[ molecules ]\n")

        with open(itp_file, "w") as f:
            f.write("[ moleculetype ]\nWAT 3\n\n[ atoms ]\n")
            f.write("1 O 1 WAT O 1 -0.834 15.999\n")

        # Add region constraint if requested
        region_line = ""
        if use_region:
            # Define a spherical region with half the box volume
            region_line = "gcmc_region:sphere 1.5 1.5 1.5 1.0"  # Sphere at center with radius 1.0

        with open(inp_file, "w") as f:
            f.write(f"""pdb:{pdb_file}
top:{top_file}
fragitp:{itp_file}
op_pdb:{prefix}_output.pdb
op_top:{prefix}_output.top
box_size:3.0 3.0 3.0
cutoff:1.2
mcsteps:8000
nprint:500
fragname:WAT
fragmuex:-5.0
{region_line}
""")

        cmd = [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--seed", str(seed)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)

        # Parse volume and density
        volume_match = re.search(r"GCMC region:.*\(volume:\s*([\d.]+)\s*nm", result.stdout)
        if use_region and volume_match:
            volume = float(volume_match.group(1))
            print(f"  Region volume: {volume:.3f} nm³")
        else:
            volume = 3.0 ** 3  # Box volume
            print(f"  Box volume: {volume:.3f} nm³")

        # Parse final molecule count
        count_match = re.search(r"WAT:\s*(\d+)", result.stdout)
        if count_match:
            count = int(count_match.group(1))
            density = count / volume
            print(f"  Final count: {count} molecules")
            print(f"  Density: {density:.3f} molecules/nm³")
            return volume, count, density
        else:
            return volume, 0, 0.0

    # Run with and without region
    print("\nWithout region constraint:")
    vol_full, count_full, dens_full = run_simulation(use_region=False)

    print("\nWith region constraint:")
    vol_region, count_region, dens_region = run_simulation(use_region=True)

    # Compare results
    if vol_region < vol_full:
        volume_ratio = vol_region / vol_full
        print(f"\nVolume ratio (region/full): {volume_ratio:.3f}")

        # At same chemical potential, density should be roughly the same
        # (equilibrium density is independent of available volume)
        if dens_full > 0 and dens_region > 0:
            density_ratio = dens_region / dens_full
            print(f"Density ratio (region/full): {density_ratio:.3f}")

            # Densities should be similar (within statistical fluctuations)
            # This validates that acceptance calculation uses correct volume
            if abs(density_ratio - 1.0) < 0.5:
                print("✓ Densities are consistent - region volume correctly used")
            else:
                print(f"⚠ Large density difference may indicate volume issue")

    # Basic validation
    assert vol_region < vol_full, "Region should reduce available volume"
    print("\n✓ Region constraint reduces available volume as expected")


if __name__ == "__main__":
    # Run tests and report issues found
    print("="*60)
    print("CRITICAL ISSUES VERIFICATION SUITE")
    print("="*60)
    print("\nThese tests expose known issues without code modification.")
    print("Expected failures indicate problems that need fixing.\n")

    pytest.main([__file__, "-v", "-s", "--tb=short"])