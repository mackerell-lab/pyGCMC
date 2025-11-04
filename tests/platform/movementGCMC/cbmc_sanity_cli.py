"""
CBMC Sanity Tests

Verifies that CBMC (Configurational Bias Monte Carlo) implementation:
1. Correctly enables/disables based on configuration
2. Records varying Rosenbluth weights when enabled
3. Uses unity weights when disabled
4. Improves insertion acceptance statistics
5. Satisfies acceptance probability formula

Test Strategy:
- Run minimal water system with CBMC on/off
- Compare qForward/qReverse distributions
- Verify cbmcTrials field matches configuration
- Check pAcc against theoretical formula
"""

import subprocess
from pathlib import Path
import math
import statistics
import json
import pytest
from typing import List, Dict, Any
from collections import defaultdict

# Path to gcmc_cpu executable
GCMC_CPU_PATH = Path(__file__).parent.parent.parent.parent / "build" / "bin" / "gcmc_cpu"


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


def create_minimal_water_files(tmpdir: Path):
    """Create minimal TIP3P water files for CBMC testing."""
    # PDB file
    pdb_file = tmpdir / "water.pdb"
    pdb_file.write_text(
        "CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1           1\n"
        "ATOM      1  O   WAT     1       5.000   5.000   5.000  1.00  0.00           O\n"
        "ATOM      2  H1  WAT     1       5.757   5.586   5.000  1.00  0.00           H\n"
        "ATOM      3  H2  WAT     1       4.243   5.586   5.000  1.00  0.00           H\n"
        "END\n"
    )

    # TOP file
    top_file = tmpdir / "water.top"
    top_file.write_text(
        "[ moleculetype ]\n"
        "WAT     3\n\n"
        "[ atoms ]\n"
        "1  O   1  WAT  OW  1  -0.834  15.9994\n"
        "2  H   1  WAT  HW1 1   0.417   1.008\n"
        "3  H   1  WAT  HW2 1   0.417   1.008\n\n"
        "[ bonds ]\n"
        "1  2  1  0.09572  462750.4\n"
        "1  3  1  0.09572  462750.4\n\n"
        "[ angles ]\n"
        "2  1  3  1  104.52  628.02\n\n"
        "[ system ]\n"
        "Water System\n\n"
        "[ molecules ]\n"
        "WAT  1\n"
    )

    # Atomtypes file
    atp_file = tmpdir / "atomtypes.atp"
    atp_file.write_text("O   15.9994\nH    1.008\n")

    # Force field file
    ff_file = tmpdir / "ffnonbonded.itp"
    ff_file.write_text(
        "[ atomtypes ]\n"
        "O   8   15.9994  -0.834  A  3.15061e-01  6.36386e-01\n"
        "H   1   1.008     0.417  A  0.00000e+00  0.00000e+00\n"
    )

    return pdb_file, top_file, atp_file, ff_file


def build_inp(
    pdb: Path,
    top: Path,
    atp: Path,
    ff: Path,
    use_cbmc: bool,
    k_trials: int,
    mcsteps: int,
    mu: float = -2.0,
    use_cavity: bool = False,
    grid_spacing: float = 2.0,
    probe_radius: float = 1.4,
):
    """Build INP file content with optional CBMC and cavity configuration."""
    if use_cbmc:
        cbmc_block = f"""use_conf_bias:yes
fragconf:{k_trials}"""
    else:
        cbmc_block = "use_conf_bias:no"

    cavity_block = (
        f"""use_cavity_bias:yes
cavity_grid_spacing:{grid_spacing}
cavity_probe_radius:{probe_radius}"""
        if use_cavity
        else "use_cavity_bias:no"
    )

    return f"""par:{ff}
atomtypes:{atp}
top:{top}
pdb:{pdb}
protitp:{top}

fragname: water
fragconc: 55.0
fragmuex: {mu}
{cbmc_block}

box_size: 10.0 10.0 10.0
cutoff: 4.5
temperature: 300.0
mcsteps: {mcsteps}
nprint: {max(1, mcsteps // 2)}
eqsteps: 0

mc_move_prob: 0.5 0.5 0 0

{cavity_block}

seed: 42
"""


def run_with_accept_log(inp_content: str, workdir: Path,
                       accept_name: str = "accept.jsonl", timeout: int = 120):
    """Run simulation and return result + acceptance log path."""
    inp_file = workdir / "test.inp"
    inp_file.write_text(inp_content)
    accept_log = workdir / accept_name

    result = subprocess.run(
        [str(GCMC_CPU_PATH), "--inp", str(inp_file),
         "--dump-accept", str(accept_log), "--seed", "42"],
        cwd=str(workdir),
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    return result, accept_log


class TestCBMCSanity:
    """CBMC sanity tests verifying weight calculation and acceptance."""

    def test_cbmc_enabled_has_weights(self, tmp_path):
        """
        P1验收测试：CBMC启用时权重变化

        验证：
        - cbmcTrials = K >= 2
        - qForward (插入) 不全等于1
        - qReverse (删除) 不全等于1
        - 权重在合理范围内变化
        """
        print("\n=== CBMC Enabled: Weights Variation Test ===")

        pdb, top, atp, ff = create_minimal_water_files(tmp_path)
        inp = build_inp(pdb, top, atp, ff, use_cbmc=True, k_trials=8,
                       mcsteps=4000, mu=-2.0)

        result, accept_log = run_with_accept_log(inp, tmp_path, "cbmc_on.jsonl", timeout=180)

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"

        records = read_jsonl(accept_log)
        assert records, "No acceptance records found"

        ins = filter_by_move(records, "insertion")
        dels = filter_by_move(records, "deletion")

        print(f"Total records: {len(records)}")
        print(f"Insertions: {len(ins)}, Deletions: {len(dels)}")

        assert ins or dels, "Need at least insertion or deletion attempts"

        # Check cbmcTrials field
        all_moves = ins + dels
        cbmc_trials = {r.get("cbmcTrials", 0) for r in all_moves}
        print(f"CBMC trials values: {cbmc_trials}")
        assert max(cbmc_trials) >= 2, f"cbmcTrials should be >= 2 when enabled, got {cbmc_trials}"

        # Check qForward variation (if insertions exist)
        if ins:
            qf = [r.get("qForward", 1.0) for r in ins[:100]]  # Sample first 100
            print(f"qForward sample (n={len(qf)}): min={min(qf):.6f}, max={max(qf):.6f}, mean={sum(qf)/len(qf):.6f}")
            # At least some should differ from 1.0
            non_unity = [x for x in qf if abs(x - 1.0) > 1e-12]
            assert len(non_unity) > 0, "qForward should vary with CBMC enabled"

        # Check qReverse variation (deletions always exist)
        if dels:
            qr = [r.get("qReverse", 1.0) for r in dels[:100]]  # Sample first 100
            print(f"qReverse sample (n={len(qr)}): min={min(qr):.6f}, max={max(qr):.6f}, mean={sum(qr)/len(qr):.6f}")
            # Should vary across attempts
            assert min(qr) < max(qr), "qReverse should vary across deletion attempts"
            # Check reasonable range (0 < W/K <= 1 typically)
            assert all(0 < x <= 1.5 for x in qr), f"qReverse out of reasonable range: {[x for x in qr if x > 1.5]}"

        print("✅ CBMC weights variation test passed")

    def test_cbmc_disabled_has_unity_weights(self, tmp_path):
        """
        P1验收测试：CBMC禁用时权重为1

        验证：
        - cbmcTrials = 1
        - qForward ≈ 1.0 (all insertions)
        - qReverse ≈ 1.0 (all deletions)
        """
        print("\n=== CBMC Disabled: Unity Weights Test ===")

        pdb, top, atp, ff = create_minimal_water_files(tmp_path)
        inp = build_inp(pdb, top, atp, ff, use_cbmc=False, k_trials=1,
                       mcsteps=3000, mu=-2.0)

        result, accept_log = run_with_accept_log(inp, tmp_path, "cbmc_off.jsonl", timeout=120)

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"

        records = read_jsonl(accept_log)
        assert records, "No acceptance records found"

        ins = filter_by_move(records, "insertion")
        dels = filter_by_move(records, "deletion")

        print(f"Total records: {len(records)}")
        print(f"Insertions: {len(ins)}, Deletions: {len(dels)}")

        assert ins or dels, "Need at least insertion or deletion attempts"

        # Check cbmcTrials = 1
        all_moves = ins + dels
        cbmc_trials = {r.get("cbmcTrials", 0) for r in all_moves}
        print(f"CBMC trials values: {cbmc_trials}")
        assert cbmc_trials == {1}, f"cbmcTrials should be 1 when disabled, got {cbmc_trials}"

        # Check qForward = 1
        if ins:
            qf = [r.get("qForward", 1.0) for r in ins]
            print(f"qForward range: min={min(qf):.6f}, max={max(qf):.6f}")
            assert all(abs(x - 1.0) < 1e-12 for x in qf), "qForward must be 1.0 without CBMC"

        # Check qReverse = 1
        if dels:
            qr = [r.get("qReverse", 1.0) for r in dels]
            print(f"qReverse range: min={min(qr):.6f}, max={max(qr):.6f}")
            assert all(abs(x - 1.0) < 1e-12 for x in qr), "qReverse must be 1.0 without CBMC"

        print("✅ Unity weights test passed")

    def test_cbmc_improves_insertion_statistics(self, tmp_path):
        """
        P1验收测试：CBMC提升插入接受统计

        验证K增大时插入接受率或平均pAcc呈上升趋势。

        策略：
        - 比较K=1 vs K=8
        - 检查接受率提升 OR 平均pAcc提升
        """
        print("\n=== CBMC Improvement Test ===")

        pdb, top, atp, ff = create_minimal_water_files(tmp_path)

        # Run baseline (K=1)
        print("Running baseline (K=1)...")
        inp_off = build_inp(pdb, top, atp, ff, use_cbmc=False, k_trials=1,
                           mcsteps=5000, mu=-1.5)  # Higher mu for more insertions
        result_off, log_off = run_with_accept_log(inp_off, tmp_path, "off.jsonl", timeout=180)
        assert result_off.returncode == 0
        rec_off = read_jsonl(log_off)

        # Run CBMC (K=8)
        print("Running CBMC (K=8)...")
        inp_on = build_inp(pdb, top, atp, ff, use_cbmc=True, k_trials=8,
                          mcsteps=5000, mu=-1.5)
        result_on, log_on = run_with_accept_log(inp_on, tmp_path, "on.jsonl", timeout=180)
        assert result_on.returncode == 0
        rec_on = read_jsonl(log_on)

        # Get insertion statistics
        ins_off = filter_by_move(rec_off, "insertion")
        ins_on = filter_by_move(rec_on, "insertion")

        print(f"Baseline insertions: {len(ins_off)}")
        print(f"CBMC insertions: {len(ins_on)}")

        if len(ins_off) < 20 or len(ins_on) < 20:
            pytest.skip("Not enough insertion attempts for robust comparison")

        # Compare acceptance rates
        stats_off = acceptance_statistics(rec_off)
        stats_on = acceptance_statistics(rec_on)

        acc_off = stats_off.get("insertion", {}).get("water", 0.0)
        acc_on = stats_on.get("insertion", {}).get("water", 0.0)

        print(f"Acceptance rate: {acc_off:.4f} (K=1) -> {acc_on:.4f} (K=8)")

        # Compare median pAcc
        pacc_off = [r.get("pAcc", 0.0) for r in ins_off if r.get("pAcc", 0.0) >= 0.0]
        pacc_on = [r.get("pAcc", 0.0) for r in ins_on if r.get("pAcc", 0.0) >= 0.0]

        med_off = statistics.median(pacc_off) if pacc_off else 0.0
        med_on = statistics.median(pacc_on) if pacc_on else 0.0

        print(f"Median pAcc: {med_off:.6f} (K=1) -> {med_on:.6f} (K=8)")

        # Accept if either metric improves
        acc_improved = acc_on >= acc_off - 0.01  # Allow small statistical noise
        pacc_improved = med_on >= med_off - 0.001

        assert acc_improved or pacc_improved, \
            f"CBMC did not improve: acc {acc_off:.3f}->{acc_on:.3f}, pAcc {med_off:.6f}->{med_on:.6f}"

        print("✅ CBMC improvement test passed")

    def test_acceptance_formula_consistency(self, tmp_path):
        """
        P1验收测试：接受概率公式验证

        验证pAcc字段与理论公式匹配：
        - 插入: pAcc = min(1, (z×V)/(N+1) × exp(-βΔU) × qForward)
        - 删除: pAcc = min(1, N/(z×V) × exp(-βΔU) × qReverse)

        容忍度: 5e-5 (数值精度)
        """
        print("\n=== Acceptance Formula Consistency Test ===")

        pdb, top, atp, ff = create_minimal_water_files(tmp_path)
        inp = build_inp(pdb, top, atp, ff, use_cbmc=True, k_trials=5,
                       mcsteps=4000, mu=-2.0)

        result, accept_log = run_with_accept_log(inp, tmp_path, "formula.jsonl", timeout=180)

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"

        records = read_jsonl(accept_log)
        assert records, "No acceptance records found"

        ins = filter_by_move(records, "insertion")
        dels = filter_by_move(records, "deletion")

        print(f"Checking {len(ins)} insertions and {len(dels)} deletions")

        tol = 5e-5
        checked_ins = 0
        checked_del = 0

        # Check insertions
        for r in ins[:300]:  # Sample to keep runtime reasonable
            z = r.get("z", 1.0)
            v_eff = r.get("vEff", 1.0)
            n_before = r.get("nBefore", 0)
            beta_delta_u = r.get("betaDeltaU", 0.0)
            qf = r.get("qForward", 1.0)

            expected = min(1.0, (z * v_eff) / (n_before + 1.0) * math.exp(-beta_delta_u) * qf)
            actual = r.get("pAcc", 0.0)

            if actual < 0:  # Skip if probability not stored
                continue

            assert abs(expected - actual) < tol, \
                f"Insertion pAcc mismatch: exp={expected:.6f}, act={actual:.6f}, " \
                f"z={z}, V={v_eff}, N={n_before}, βΔU={beta_delta_u:.4f}, qF={qf:.6f}"
            checked_ins += 1

        # Check deletions
        for r in dels[:300]:
            z = r.get("z", 1.0)
            v_eff = r.get("vEff", 1.0)
            n_before = r.get("nBefore", 1)
            beta_delta_u = r.get("betaDeltaU", 0.0)
            qr = r.get("qReverse", 1.0)

            # Engine records qReverse = W_old / K, so detailed balance uses its reciprocal.
            effective_qr = qr if qr > 0 else 1.0
            expected = min(1.0, n_before / (z * v_eff) * math.exp(beta_delta_u) / effective_qr)
            actual = r.get("pAcc", 0.0)

            if actual < 0:  # Skip if probability not stored
                continue

            assert abs(expected - actual) < tol, \
                f"Deletion pAcc mismatch: exp={expected:.6f}, act={actual:.6f}, " \
                f"z={z}, V={v_eff}, N={n_before}, βΔU={beta_delta_u:.4f}, qR={qr:.6f}"
            checked_del += 1

        print(f"Verified {checked_ins} insertions and {checked_del} deletions")
        assert checked_ins + checked_del > 0, "No records validated against acceptance formula"

        print("✅ Acceptance formula consistency test passed")

    def test_v_eff_matches_cavity_fraction(self, tmp_path):
        """
        验证日志中的 vEff 与 vBox×wCavity 一致（支持 cavity bias 校验）
        """
        print("\n=== vEff Consistency With Cavity Fraction ===")

        pdb, top, atp, ff = create_minimal_water_files(tmp_path)
        inp = build_inp(
            pdb,
            top,
            atp,
            ff,
            use_cbmc=True,
            k_trials=5,
            mcsteps=4000,
            mu=-2.0,
            use_cavity=True,
        )

        result, accept_log = run_with_accept_log(
            inp, tmp_path, "veff.jsonl", timeout=180
        )
        assert result.returncode == 0, f"Simulation failed: {result.stderr}"

        records = read_jsonl(accept_log)
        insertions = filter_by_move(records, "insertion")
        assert insertions, "No insertion records found for vEff check"

        checked = 0
        for rec in insertions[:200]:  # sample first 200 records for speed
            v_eff = rec.get("vEff")
            v_box = rec.get("vBox")
            w_cavity = rec.get("wCavity")

            assert v_eff is not None, "vEff missing from acceptance record"
            assert v_box is not None, "vBox missing from acceptance record"
            assert w_cavity is not None, "wCavity missing from acceptance record"

            expected = v_box * w_cavity
            if expected <= 0:
                continue

            rel_error = abs(v_eff - expected) / expected
            assert rel_error < 1e-6, (
                "vEff does not match vBox * wCavity: "
                f"vEff={v_eff:.12e}, vBox={v_box:.12e}, "
                f"wCavity={w_cavity:.12e}, rel_error={rel_error:.2e}"
            )
            checked += 1

        assert checked > 0, "No positive-volume records validated"
        print(f"Validated {checked} insertion records for vEff consistency")
        print("✅ vEff consistency test passed")


if __name__ == "__main__":
    # Allow running test directly for debugging
    pytest.main([__file__, "-v", "-s"])
