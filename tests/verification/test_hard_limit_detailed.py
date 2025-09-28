#!/usr/bin/env python3
"""
Detailed test for hard molecule limit issue
"""

import subprocess
import tempfile
import re
from pathlib import Path

GCMC_CPU_PATH = Path("/home/zhaomt/gcmc/test108/pygcmc_dev/build/bin/gcmc_cpu")


def test_hard_limit_with_extreme_mu():
    """Test with extremely high chemical potential to hit the 1000 molecule limit"""
    print("\n=== DETAILED HARD LIMIT TEST ===")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)

        # Create minimal system
        pdb_file = tmppath / "test.pdb"
        top_file = tmppath / "test.top"
        itp_file = tmppath / "water.itp"
        inp_file = tmppath / "test.inp"

        with open(pdb_file, "w") as f:
            f.write("TITLE Test\n")
            f.write("CRYST1  200.000  200.000  200.000  90.00  90.00  90.00 P 1           1\n")
            f.write("END\n")

        with open(top_file, "w") as f:
            f.write("[ system ]\nTest\n\n[ molecules ]\n")

        with open(itp_file, "w") as f:
            f.write("[ moleculetype ]\nWAT 3\n\n[ atoms ]\n")
            f.write("1 O 1 WAT O 1 -0.01 15.999\n")  # Very weak interaction

        # Test with different chemical potentials
        mu_values = [5.0, 10.0, 15.0, 20.0]
        results = []

        for mu in mu_values:
            print(f"\nTesting with μ = {mu} kJ/mol")

            with open(inp_file, "w") as f:
                f.write(f"""pdb:{pdb_file}
top:{top_file}
fragitp:{itp_file}
op_pdb:output.pdb
op_top:output.top
box_size:20.0 20.0 20.0
cutoff:5.0
mcsteps:10000
nprint:2000
fragname:WAT
fragmuex:{mu}
""")

            cmd = [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--seed", "42"]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)

            # Parse results
            count_match = re.search(r"Fragment counts:\s*WAT:\s*(\d+)", result.stdout)
            insert_rate_match = re.search(r"Insert move accept:\s*([\d.]+)%", result.stdout)
            insert_attempts_match = re.search(r"Insert attempts:\s*(\d+)", result.stdout)
            insert_accepted_match = re.search(r"Insert accepted:\s*(\d+)", result.stdout)

            if count_match:
                count = int(count_match.group(1))
                insert_rate = float(insert_rate_match.group(1)) if insert_rate_match else 0
                attempts = int(insert_attempts_match.group(1)) if insert_attempts_match else 0
                accepted = int(insert_accepted_match.group(1)) if insert_accepted_match else 0

                results.append({
                    "mu": mu,
                    "count": count,
                    "insert_rate": insert_rate,
                    "attempts": attempts,
                    "accepted": accepted
                })

                print(f"  Count: {count}")
                print(f"  Insert rate: {insert_rate:.1f}%")
                print(f"  Insert attempts: {attempts}, accepted: {accepted}")

                # Check if we're hitting the limit
                if 990 <= count <= 1010:
                    print(f"  ⚠️  Count near 1000 limit!")
                    if insert_rate < 10:
                        print(f"  ❌ HARD LIMIT CONFIRMED: Low acceptance despite high μ")

        # Analyze trend
        print("\n=== ANALYSIS ===")
        print("μ (kJ/mol) | Count | Insert Rate")
        print("-" * 35)
        for r in results:
            marker = "⚠️" if 990 <= r["count"] <= 1010 else ""
            print(f"{r['mu']:5.1f}      | {r['count']:5d} | {r['insert_rate']:6.1f}%  {marker}")

        # Check for saturation at 1000
        high_mu_counts = [r["count"] for r in results if r["mu"] >= 10]
        if high_mu_counts and all(990 <= c <= 1010 for c in high_mu_counts):
            print("\n❌ HARD LIMIT CONFIRMED: All high-μ runs saturate near 1000")
            print("   This is non-physical and due to hard-coded limit in GCMCEngine.cpp:106")
        elif any(c > 1010 for c in high_mu_counts):
            print("\n✓ No hard limit detected - molecules can exceed 1000")
        else:
            print("\n? Results inconclusive")


if __name__ == "__main__":
    test_hard_limit_with_extreme_mu()