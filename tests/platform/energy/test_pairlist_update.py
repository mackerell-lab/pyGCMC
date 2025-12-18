"""
Test neighbor list (pairlist) functionality.

gcmc_cpu parses these keys for compatibility, but does not implement
pairlist scheduling yet. The keys must be surfaced as ignored_inp_keys
to avoid silent behavior drift.
"""

import json
import subprocess
from pathlib import Path

import pytest

GCMC_CPU_PATH = Path(__file__).parent.parent.parent.parent / "build" / "bin" / "gcmc_cpu"


class TestPairlistUpdate:
    """Test neighbor list INP compatibility (ignored_inp_keys)."""

    @staticmethod
    def create_test_system(
        tmpdir,
        box_size=30.0,
        n_molecules=10,
        mcsteps=100,
        seed=42,
        pairlist_cutoff=None,
        pairlist_cutoff_protein=None,
        pairlist_freq=None,
        use_group_cutoff=None,
    ):
        """Create test system for pairlist compatibility checks."""

        pdb_file = tmpdir / "test.pdb"
        pdb_content = (
            f"CRYST1   {box_size:.3f}   {box_size:.3f}   {box_size:.3f}  90.00  90.00  90.00 P 1           1\n"
        )

        atom_id = 1
        for i in range(n_molecules):
            x = (i % 3) * 10.0 + 5.0
            y = ((i // 3) % 3) * 10.0 + 5.0
            z = (i // 9) * 10.0 + 5.0

            pdb_content += f"ATOM  {atom_id:5d}  O   WAT  {i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00\n"
            atom_id += 1
            pdb_content += f"ATOM  {atom_id:5d}  H1  WAT  {i+1:4d}    {x+0.757:8.3f}{y+0.586:8.3f}{z:8.3f}  1.00  0.00\n"
            atom_id += 1
            pdb_content += f"ATOM  {atom_id:5d}  H2  WAT  {i+1:4d}    {x-0.757:8.3f}{y+0.586:8.3f}{z:8.3f}  1.00  0.00\n"
            atom_id += 1

        pdb_content += "END\n"
        pdb_file.write_text(pdb_content)

        top_file = tmpdir / "test.top"
        top_content = """[ defaults ]
1 2 yes 0.5 0.8333

[ atomtypes ]
O   8   15.9994  -0.834  A  3.15061e-01  6.36386e-01
H   1   1.008     0.417  A  0.00000e+00  0.00000e+00

[ moleculetype ]
WAT   3

[ atoms ]
1  O   1  WAT  O   1  -0.834  15.9994
2  H   1  WAT  H1  2   0.417   1.008
3  H   1  WAT  H2  3   0.417   1.008

[ bonds ]
1  2  1  0.09572  502416.0
1  3  1  0.09572  502416.0

[ angles ]
2  1  3  1  104.52  628.02

[ system ]
Test Water

[ molecules ]
"""
        for _ in range(n_molecules):
            top_content += "WAT  1\n"

        top_file.write_text(top_content)

        atp_file = tmpdir / "atomtypes.atp"
        atp_file.write_text("O   15.9994\nH    1.008\n")

        ff_file = tmpdir / "ffnonbonded.itp"
        ff_file.write_text("""[ atomtypes ]
O   8   15.9994  -0.834  A  3.15061e-01  6.36386e-01
H   1   1.008     0.417  A  0.00000e+00  0.00000e+00
""")

        inp_file = tmpdir / "test.inp"
        inp_content = f"""par:{ff_file}
atomtypes:{atp_file}
top:{top_file}
pdb:{pdb_file}
protitp:{top_file}

fragname: water
fragconc: 55.0
fragmuex: -5.60

box_size: {box_size} {box_size} {box_size}
cutoff: {box_size/2 - 1.0}
temperature: 300
mcsteps: {mcsteps}
nprint: 100
eqsteps: 0

mc_move_prob: 0 0 0.5 0.5

seed: {seed}
"""

        if pairlist_cutoff is not None:
            inp_content += f"pairlist_cutoff: {pairlist_cutoff}\n"
        if pairlist_cutoff_protein is not None:
            inp_content += f"pairlist_cutoff_protein: {pairlist_cutoff_protein}\n"
        if pairlist_freq is not None:
            inp_content += f"pairlist_freq: {pairlist_freq}\n"
        if use_group_cutoff is not None:
            if isinstance(use_group_cutoff, str):
                value = use_group_cutoff
            else:
                value = "yes" if use_group_cutoff else "no"
            inp_content += f"use_group_cutoff: {value}\n"

        inp_file.write_text(inp_content)

        return inp_file

    @staticmethod
    def run_dump_params(inp_file, workdir):
        params_json = workdir / "params.json"
        out_prefix = workdir / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)

        result = subprocess.run(
            [
                str(GCMC_CPU_PATH),
                "--inp",
                str(inp_file),
                "--prefix",
                str(out_prefix),
                "--dump-params",
                str(params_json),
            ],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(workdir),
        )
        assert result.returncode == 0, result.stdout + result.stderr
        return json.loads(params_json.read_text())

    @staticmethod
    def assert_pairlist_ignored(params, keys):
        ignored = set(params["basic"]["ignored_inp_keys"])
        unknown = set(params["basic"]["unknown_inp_keys"])
        for key in keys:
            assert key in ignored
            assert key not in unknown

    def test_pairlist_cutoff_effect(self, tmp_path):
        """
        pairlist_cutoff 作为兼容输入存在，但当前应被标记为 ignored。
        """
        if not GCMC_CPU_PATH.exists():
            pytest.skip(f"gcmc_cpu not built: {GCMC_CPU_PATH}")

        run1_dir = tmp_path / "run_large_cutoff"
        run1_dir.mkdir()
        inp_file1 = self.create_test_system(
            run1_dir,
            box_size=30.0,
            n_molecules=5,
            mcsteps=0,
            pairlist_cutoff=2.0,
            seed=42,
        )

        params1 = self.run_dump_params(inp_file1, run1_dir)
        self.assert_pairlist_ignored(params1, ["pairlist_cutoff"])

        run2_dir = tmp_path / "run_small_cutoff"
        run2_dir.mkdir()
        inp_file2 = self.create_test_system(
            run2_dir,
            box_size=30.0,
            n_molecules=5,
            mcsteps=0,
            pairlist_cutoff=0.8,
            seed=42,
        )

        params2 = self.run_dump_params(inp_file2, run2_dir)
        self.assert_pairlist_ignored(params2, ["pairlist_cutoff"])

    def test_pairlist_update_frequency(self, tmp_path):
        """
        pairlist_freq 作为兼容输入存在，但当前应被标记为 ignored。
        """
        if not GCMC_CPU_PATH.exists():
            pytest.skip(f"gcmc_cpu not built: {GCMC_CPU_PATH}")

        run1_dir = tmp_path / "run_freq_high"
        run1_dir.mkdir()
        inp_file1 = self.create_test_system(
            run1_dir,
            box_size=30.0,
            n_molecules=5,
            mcsteps=0,
            pairlist_freq=10,
            seed=42,
        )

        params1 = self.run_dump_params(inp_file1, run1_dir)
        self.assert_pairlist_ignored(params1, ["pairlist_freq"])

        run2_dir = tmp_path / "run_freq_low"
        run2_dir.mkdir()
        inp_file2 = self.create_test_system(
            run2_dir,
            box_size=30.0,
            n_molecules=5,
            mcsteps=0,
            pairlist_freq=1000,
            seed=42,
        )

        params2 = self.run_dump_params(inp_file2, run2_dir)
        self.assert_pairlist_ignored(params2, ["pairlist_freq"])

    def test_basic_simulation_with_pairlist(self, tmp_path):
        """
        兼容性 smoke：带有 pairlist 相关参数的模拟可正常运行，且这些参数被标记 ignored。
        """
        if not GCMC_CPU_PATH.exists():
            pytest.skip(f"gcmc_cpu not built: {GCMC_CPU_PATH}")

        inp_file = self.create_test_system(
            tmp_path,
            box_size=30.0,
            n_molecules=5,
            mcsteps=50,
            pairlist_cutoff=1.5,
            pairlist_cutoff_protein=1.8,
            pairlist_freq=50,
            use_group_cutoff=False,
            seed=42,
        )

        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_file), "--seed", "42"],
            capture_output=True,
            text=True,
            timeout=60,
            cwd=str(tmp_path),
        )

        assert result.returncode == 0, f"Simulation failed: {result.stderr}"

        final_pdb = tmp_path / "gcmc_final.pdb"
        assert final_pdb.exists(), "Output PDB not found"

        stats_files = list(tmp_path.glob("*statistics.dat"))
        assert len(stats_files) > 0, "No statistics.dat file"

        params = self.run_dump_params(inp_file, tmp_path)
        self.assert_pairlist_ignored(
            params,
            [
                "pairlist_cutoff",
                "pairlist_cutoff_protein",
                "pairlist_freq",
                "use_group_cutoff",
            ],
        )
