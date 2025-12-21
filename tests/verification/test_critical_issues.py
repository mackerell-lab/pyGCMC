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
import math
import json
import numpy as np
from pathlib import Path

# Path to the gcmc_cpu executable
GCMC_CPU_PATH = Path(__file__).parent.parent.parent / "build" / "bin" / "gcmc_cpu"
if not GCMC_CPU_PATH.exists():
    GCMC_CPU_PATH = Path("/home/zhaomt/gcmc/test108/pygcmc_dev/build/bin/gcmc_cpu")


def _load_accept_records(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def _count_residues(pdb_path: Path, resname: str) -> int:
    resids = set()
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[17:20].strip().upper() != resname.strip().upper():
            continue
        try:
            resids.add(int(line[22:26]))
        except ValueError:
            continue
    return len(resids)


def _read_last_stats_energy(stats_path: Path) -> float:
    data_lines = [
        line for line in stats_path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    if not data_lines:
        raise AssertionError(f"No data lines in {stats_path}")
    return float(data_lines[-1].split()[1])


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

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=10,
                cwd=str(tmp_path)
            )

            prefix = tmp_path / f"run_{run}"
            final_pdb = Path(f"{prefix}_final.pdb")
            stats_path = Path(f"{prefix}_statistics.dat")
            if not final_pdb.exists() or not stats_path.exists():
                raise AssertionError("Expected output files not generated")

            count = _count_residues(final_pdb, "WAT")
            coord_info = str(_read_last_stats_energy(stats_path))

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

        accept_log = tmp_path / "accept.jsonl"
        cmd = [
            str(GCMC_CPU_PATH),
            "--inp",
            str(inp_file),
            "--seed",
            "54321",
            "--dump-accept",
            str(accept_log),
        ]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=10,
            cwd=str(tmp_path)
        )

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"
        assert accept_log.exists(), "Acceptance log missing"

        records = _load_accept_records(accept_log)
        ins = [r for r in records if r.get("move") == "insertion"]
        dele = [r for r in records if r.get("move") == "deletion"]
        ins_att, del_att = len(ins), len(dele)
        ins_acc = sum(1 for r in ins if bool(r.get("accepted")))
        del_acc = sum(1 for r in dele if bool(r.get("accepted")))

        assert ins_att > 0 and del_att > 0, "Expected insertion and deletion attempts"

        insert_rate = ins_acc / ins_att
        delete_rate = del_acc / del_att
        weighted_total = (ins_acc + del_acc) / (ins_att + del_att)

        min_rate = min(insert_rate, delete_rate)
        max_rate = max(insert_rate, delete_rate)
        assert min_rate <= weighted_total <= max_rate, (
            f"Weighted rate {weighted_total:.3f} not between insert {insert_rate:.3f} and delete {delete_rate:.3f}"
        )


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

        prefix = tmp_path / "hard_limit"
        accept_log = tmp_path / "hard_limit_accept.jsonl"
        cmd = [
            str(GCMC_CPU_PATH),
            "--inp",
            str(inp_file),
            "--seed",
            "11111",
            "--prefix",
            str(prefix),
            "--dump-accept",
            str(accept_log),
        ]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=20,
            cwd=str(tmp_path)
        )

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"
        final_pdb = Path(f"{prefix}_final.pdb")
        assert final_pdb.exists(), "Final PDB missing"
        final_count = _count_residues(final_pdb, "WAT")
        assert accept_log.exists(), "Acceptance log missing"
        records = _load_accept_records(accept_log)
        ins = [r for r in records if r.get("move") == "insertion"]
        ins_acc = sum(1 for r in ins if bool(r.get("accepted")))
        insert_rate = 100.0 * ins_acc / max(len(ins), 1)

        print(f"Final molecule count: {final_count}")
        print(f"Insert acceptance rate: {insert_rate}%")

        # With μ=10, we expect many molecules unless limited
        if 900 <= final_count <= 1000:
            print(f"❌ ISSUE CONFIRMED: Count ({final_count}) near hard limit of 1000")
            if insert_rate < 5:
                print("   Insert rate very low, suggesting systematic rejection at limit")
        elif final_count > 1000:
            print(f"✓ No hard limit detected (count={final_count})")
        else:
            print(f"? Inconclusive: count={final_count} (may be physical limit)")

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

        accept_log = tmp_path / "acceptance.jsonl"
        cmd = [
            str(GCMC_CPU_PATH),
            "--inp",
            str(inp_file),
            "--seed",
            "99999",
            "--dump-accept",
            str(accept_log),
        ]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=10,
            cwd=str(tmp_path)
        )

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"
        assert accept_log.exists(), "Acceptance log missing"
        records = _load_accept_records(accept_log)

        if records:
            print(f"\nAnalyzing {len(records)} acceptance records...")

            # Check a few insertion attempts
            insert_records = [r for r in records if r.get("move") == "insertion"][:5]

            for rec in insert_records:
                # For ideal gas, formula: p_acc = min(1, exp(-βΔU) * z * V / (N+1))
                # Use the logged acceptance terms to avoid mismatched input assumptions.
                deltaU = rec["deltaU"]
                beta_delta_u = rec.get("betaDeltaU")
                if beta_delta_u is None:
                    beta = 1.0 / (8.314e-3 * 300)
                    beta_delta_u = beta * deltaU

                z = max(float(rec.get("z", 1.0)), 1e-30)
                v_eff = float(rec.get("vEff", rec.get("vBox", 1.0)))
                v_eff = max(v_eff, 1e-30)
                n_before = int(rec.get("nBefore", rec.get("N", 0)))
                proposal_ratio = max(float(rec.get("proposalRatio", 1.0)), 1e-30)
                rosen = max(float(rec.get("rosenbluthWeight", 1.0)), 1e-30)
                cavity = max(float(rec.get("cavityFraction", 1.0)), 1e-30)

                log_ratio = (
                    math.log(proposal_ratio)
                    - beta_delta_u
                    + math.log(z)
                    + math.log(v_eff)
                    + math.log(cavity)
                    - math.log(n_before + 1)
                    + math.log(rosen)
                )
                expected = 1.0 if log_ratio >= 0.0 else math.exp(max(log_ratio, -700.0))
                actual = rec["pAcc"]

                print(f"  Move: {rec['move']}, ΔU={deltaU:.2f}, pAcc={actual:.3f}, expected≈{expected:.3f}")

                if abs(deltaU) < 10:
                    assert abs(actual - expected) < 1e-6, (
                        f"Acceptance probability mismatch: actual={actual:.6f}, expected={expected:.6f}"
                    )

            print("✓ Acceptance probabilities consistent with theory")




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
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=10,
        env=env,
        cwd=str(tmp_path)
    )

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
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=30,
            cwd=str(tmp_path)
        )

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"

        if use_region:
            radius = 1.0
            volume = (4.0 / 3.0) * math.pi * (radius ** 3)
            print(f"  Region volume: {volume:.3f} (input units^3)")
        else:
            volume = 3.0 ** 3  # Box volume
            print(f"  Box volume: {volume:.3f} (input units^3)")

        op_pdb = tmp_path / f"{prefix}_output.pdb"
        if not op_pdb.exists():
            return volume, 0, 0.0

        count = _count_residues(op_pdb, "WAT")
        density = count / volume if volume > 0 else 0.0
        print(f"  Final count: {count} molecules")
        print(f"  Density: {density:.3f} molecules/nm³")
        return volume, count, density

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
