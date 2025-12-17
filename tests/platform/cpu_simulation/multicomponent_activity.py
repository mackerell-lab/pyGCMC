"""
Test multi-component chemical potential consistency.

Verifies P1: 多组分化学势一致性
- In ideal gas limit, N_i/N_j ≈ exp[β(μ_i−μ_j)]
- Activity ratio matches theoretical prediction
"""

import pytest
import subprocess
from pathlib import Path
import json
import math


# Path to gcmc_cpu executable
GCMC_CPU_PATH = Path(__file__).parent.parent.parent.parent / "build" / "bin" / "gcmc_cpu"


class TestMulticomponentActivity:
    """Test activity ratios in multi-component systems"""

    @staticmethod
    def _first_accept_record(path: Path, *, move: str, species: str) -> dict:
        records = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
        want_move = move.strip().lower()
        want_species = species.strip().upper()
        for rec in records:
            if str(rec.get("move", "")).strip().lower() != want_move:
                continue
            if str(rec.get("species", "")).strip().upper() != want_species:
                continue
            return rec
        raise AssertionError(f"No {move}/{species} record found in {path}; first records: {records[:3]}")

    def test_activity_ratio_ideal(self, tmp_path):
        """
        P1验收测试：多组分化学势一致性（理想气体近似）

        验证：
        - 使用 --dump-accept 的结构化字段（不依赖 stdout 日志文本）
        - 两组分在 μVT（理想气体近似）下：z_i / z_j = exp[β(μ_i - μ_j)]
          其中 z 来自 acceptance JSONL 的字段（内部单位：kJ/mol 与 nm）
        """
        work = tmp_path / "multi_component_activity"
        work.mkdir(parents=True, exist_ok=True)

        test_data_dir = Path(__file__).parent.parent.parent / "data"
        na_itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
        cl_itp = test_data_dir / "charmm36.ff" / "mol" / "cl.itp"
        assert na_itp.exists()
        assert cl_itp.exists()

        # μ in kJ/mol (internal); INP uses kcal/mol by default (gcmc_gpu/opencl convention).
        mu_na_kj = -5.00
        mu_cl_kj = -7.00
        temperature_k = 300.0
        beta = 1.0 / (8.314e-3 * temperature_k)  # mol/kJ
        expected_ratio = math.exp(beta * (mu_na_kj - mu_cl_kj))

        kj_to_kcal = 1.0 / 4.184
        mu_na_kcal = mu_na_kj * kj_to_kcal
        mu_cl_kcal = mu_cl_kj * kj_to_kcal

        out_prefix = work / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        accept_log = work / "out" / "acceptance.jsonl"

        # Default INP semantics: Å for lengths and kcal/mol for fragmuex (no need to declare inp_units:gcmc_gpu).
        inp_file = work / "test.inp"
        inp_file.write_text(
            f"""
fragitp:{na_itp}
fragitp:{cl_itp}
fragname:NA CL
fragmuex:{mu_na_kcal:.8f} {mu_cl_kcal:.8f}
mctime:1 1

box_size:100.0 100.0 100.0
cutoff:12.0
temperature:{temperature_k}
moves_per_step:1
mcsteps:200
nprint:1000
mc_move_prob:1 0 0 0
""".strip()
            + "\n"
        )

        result = subprocess.run(
            [
                str(GCMC_CPU_PATH),
                "--inp",
                str(inp_file),
                "--prefix",
                str(out_prefix),
                "--dump-accept",
                str(accept_log),
                "--seed",
                "123",
            ],
            cwd=str(work),
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        assert accept_log.exists()

        rec_na = self._first_accept_record(accept_log, move="insertion", species="NA")
        rec_cl = self._first_accept_record(accept_log, move="insertion", species="CL")

        # Unit conversion check (kcal -> kJ in acceptance log).
        assert float(rec_na["mu"]) == pytest.approx(mu_na_kj, rel=1e-6, abs=1e-6)
        assert float(rec_cl["mu"]) == pytest.approx(mu_cl_kj, rel=1e-6, abs=1e-6)

        ratio = float(rec_na["z"]) / float(rec_cl["z"])
        assert ratio == pytest.approx(expected_ratio, rel=1e-6, abs=1e-12)

    def test_activity_ratio_conc_and_muex(self, tmp_path):
        """
        When both concentration and excess chemical potential are supplied, acceptance uses:
            z = conc(M) * NA_CONV * exp(beta * mu_ex)
        so z_i/z_j = (conc_i/conc_j) * exp(beta*(mu_i-mu_j)).
        """
        work = tmp_path / "multi_component_activity_conc_mu"
        work.mkdir(parents=True, exist_ok=True)

        test_data_dir = Path(__file__).parent.parent.parent / "data"
        na_itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
        cl_itp = test_data_dir / "charmm36.ff" / "mol" / "cl.itp"
        assert na_itp.exists()
        assert cl_itp.exists()

        mu_na_kj = -3.0
        mu_cl_kj = -1.0
        conc_na = 0.20
        conc_cl = 0.10
        temperature_k = 300.0
        beta = 1.0 / (8.314e-3 * temperature_k)  # mol/kJ

        expected_ratio = (conc_na / conc_cl) * math.exp(beta * (mu_na_kj - mu_cl_kj))

        kj_to_kcal = 1.0 / 4.184
        mu_na_kcal = mu_na_kj * kj_to_kcal
        mu_cl_kcal = mu_cl_kj * kj_to_kcal

        out_prefix = work / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        accept_log = work / "out" / "acceptance.jsonl"

        inp_file = work / "test.inp"
        inp_file.write_text(
            f"""
fragitp:{na_itp}
fragitp:{cl_itp}
fragname:NA CL
fragconc:{conc_na:.6f} {conc_cl:.6f}
fragmuex:{mu_na_kcal:.8f} {mu_cl_kcal:.8f}
mctime:1 1

box_size:100.0 100.0 100.0
cutoff:12.0
temperature:{temperature_k}
moves_per_step:1
mcsteps:200
nprint:1000
mc_move_prob:1 0 0 0
""".strip()
            + "\n"
        )

        result = subprocess.run(
            [
                str(GCMC_CPU_PATH),
                "--inp",
                str(inp_file),
                "--prefix",
                str(out_prefix),
                "--dump-accept",
                str(accept_log),
                "--seed",
                "321",
            ],
            cwd=str(work),
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        assert accept_log.exists()

        rec_na = self._first_accept_record(accept_log, move="insertion", species="NA")
        rec_cl = self._first_accept_record(accept_log, move="insertion", species="CL")

        assert float(rec_na["mu"]) == pytest.approx(mu_na_kj, rel=1e-6, abs=1e-6)
        assert float(rec_cl["mu"]) == pytest.approx(mu_cl_kj, rel=1e-6, abs=1e-6)

        ratio = float(rec_na["z"]) / float(rec_cl["z"])
        assert ratio == pytest.approx(expected_ratio, rel=1e-6, abs=1e-12)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
