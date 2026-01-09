"""
Multi-fragment physical contract checks (CLI-driven).
"""

from __future__ import annotations

import math
import statistics
from collections import Counter
from pathlib import Path

import pytest

from .physical_contracts import (
    _build_state_from_pdb_top,
    _coulomb_energy_kj_mol,
    _first_accept_record,
    _first_atom_xyz_angstrom,
    _lj_energy_kj_mol,
    _load_accept_records,
    _min_image_delta_nm,
    _read_cryst1_box_angstrom,
    _run_gcmc_cpu,
    _write_inp,
    _write_text,
)


def test_multifragment_cross_species_deltaU_matches_energy_components(gcmc_cpu, temp_dir):
    """
    Cross-fragment insertion: deltaU should match analytic LJ+Coulomb and backend energy diff.
    """
    work = Path(temp_dir) / "physical_contracts" / "multi_fragment_deltaU"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
ATOM      1  A   FIX A   1      15.000  15.000  15.000  1.00  0.00           C
END
""",
    )

    top_ins = work / "sys_ins.top"
    _write_text(
        top_ins,
        """
[ defaults ]
1 2 yes 0.5 0.8333

[ moleculetype ]
FIX  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   A     1      FIX       A     1     0.500   12.011

[ system ]
MultiFrag

[ molecules ]
FIX  1
""",
    )

    top_del = work / "sys_del.top"
    _write_text(
        top_del,
        """
[ defaults ]
1 2 yes 0.5 0.8333

[ moleculetype ]
FIX  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   A     1      FIX       A     1     0.500   12.011

[ moleculetype ]
FRB  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   B     1      FRB       B     1    -0.250   14.000

[ system ]
MultiFrag

[ molecules ]
FIX  1
FRB  1
""",
    )

    par = work / "par.itp"
    _write_text(
        par,
        """
[ defaults ]
1 2 yes 0.5 0.8333

[ atomtypes ]
; name  at.num  mass   charge  ptype  sigma   epsilon
A       0       12.011 0.000   A      0.300   0.200
B       0       14.000 0.000   A      0.280   0.100
""",
    )

    frag_a = work / "frag_a.itp"
    _write_text(
        frag_a,
        """
[ moleculetype ]
FRA  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   A     1      FRA       A     1     0.500   12.011
""",
    )

    frag_b = work / "frag_b.itp"
    _write_text(
        frag_b,
        """
[ moleculetype ]
FRB  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   B     1      FRB       B     1    -0.250   14.000
""",
    )

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "run.inp"
    _write_inp(
        inp,
        f"""
random_seed:13579
par:{par}
fragitp:{frag_a}
fragitp:{frag_b}
fragname:FRA FRB
fragconc:1.0 1.0
fragmuex:0.0 30.0
mctime:0 1

pdb:{pdb}
top:{top_ins}
box_size:30.0 30.0 30.0
gcmc_region:box 18.0 15.0 15.0 19.0 16.0 16.0
cutoff:6.0
temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:1 0 0 0
""",
    )

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work,
        inp=inp,
        out_prefix=out_prefix,
        extra_args=["--dump-accept", str(accept_log)],
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="insertion", species="FRB")
    assert bool(rec.get("accepted")) is True

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()

    box_ang = _read_cryst1_box_angstrom(final_pdb)
    a_pos = _first_atom_xyz_angstrom(final_pdb, resname="FIX")
    b_pos = _first_atom_xyz_angstrom(final_pdb, resname="FRB")
    dx, dy, dz = _min_image_delta_nm(a_pos, b_pos, box_ang=box_ang)
    r_nm = math.sqrt(dx * dx + dy * dy + dz * dz)

    cutoff_nm = 6.0 / 10.0
    assert r_nm < cutoff_nm

    sigma_mix = 0.5 * (0.300 + 0.280)
    eps_mix = math.sqrt(0.200 * 0.100)
    expected = _lj_energy_kj_mol(
        r_nm=r_nm, sigma_nm=sigma_mix, eps_kj_mol=eps_mix
    ) + _coulomb_energy_kj_mol(r_nm=r_nm, q1=0.500, q2=-0.250)
    assert float(rec["deltaU"]) == pytest.approx(expected, rel=1e-3, abs=1e-2)

    box_nm = (3.0, 3.0, 3.0)
    mc_before, state_before = _build_state_from_pdb_top(
        par,
        pdb,
        top_ins,
        cutoff_nm=cutoff_nm,
        box_nm=box_nm,
    )
    import pygcmc

    pygcmc.computeSystemEnergyPBCCutoff(state_before)
    elec_before, vdw_before = pygcmc.getTotalEnergyComponents(state_before)

    mc_after, state_after = _build_state_from_pdb_top(
        par,
        final_pdb,
        top_del,
        cutoff_nm=cutoff_nm,
        box_nm=box_nm,
    )
    pygcmc.computeSystemEnergyPBCCutoff(state_after)
    elec_after, vdw_after = pygcmc.getTotalEnergyComponents(state_after)

    delta_components = (elec_after - elec_before) + (vdw_after - vdw_before)
    assert float(rec["deltaU"]) == pytest.approx(delta_components, rel=5e-3, abs=0.15)


def _n_before_samples(records: list[dict], species: str) -> list[int]:
    samples: list[int] = []
    want = species.strip().upper()
    for rec in records:
        if str(rec.get("species", "")).strip().upper() != want:
            continue
        move = str(rec.get("move", "")).strip().lower()
        if move not in ("insertion", "deletion"):
            continue
        samples.append(int(rec.get("nBefore", 0)))
    return samples


def test_multifragment_poisson_distribution_ideal_gas(gcmc_cpu, temp_dir):
    """
    Two-fragment ideal gas: each species should follow Poisson statistics independently.
    """
    work = Path(temp_dir) / "physical_contracts" / "poisson_cli_multi"
    work.mkdir(parents=True, exist_ok=True)

    par = work / "par.itp"
    _write_text(
        par,
        """
[ defaults ]
1 2 yes 0.5 0.8333

[ atomtypes ]
; name  at.num  mass   charge  ptype  sigma   epsilon
X       0       1.000  0.000   A      0.300   0.000
Y       0       1.000  0.000   A      0.300   0.000
""",
    )

    frag_x = work / "frag_x.itp"
    _write_text(
        frag_x,
        """
[ moleculetype ]
X  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   X     1      X        X     1     0.000   1.000
""",
    )

    frag_y = work / "frag_y.itp"
    _write_text(
        frag_y,
        """
[ moleculetype ]
Y  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   Y     1      Y        Y     1     0.000   1.000
""",
    )

    box_ang = 10.0
    v_box_expected = (box_ang / 10.0) ** 3
    beta = 1.0 / (0.008314462618 * 300.0)
    target_mean_x = 4.0
    target_mean_y = 2.0
    target_z_x = target_mean_x / v_box_expected
    target_z_y = target_mean_y / v_box_expected
    mu_x_kj = math.log(target_z_x) / beta
    mu_y_kj = math.log(target_z_y) / beta
    mu_x_kcal = mu_x_kj / 4.184
    mu_y_kcal = mu_y_kj / 4.184

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "run.inp"
    mcsteps = 12000
    _write_inp(
        inp,
        f"""
random_seed:24680
par:{par}
fragitp:{frag_x}
fragitp:{frag_y}
fragname:X Y
fragconc:0.0 0.0
fragmuex:{mu_x_kcal} {mu_y_kcal}
# Non-uniform fragment selection weights must not change μVT equilibrium (proposal cancels).
mctime:2 1

box_size:10.0 10.0 10.0
gcmc_region:box 0 0 0 10 10 10
cutoff:4.0
temperature:300.0
moves_per_step:1
mcsteps:{mcsteps}
nprint:1
mc_move_prob:1 1 0 0
""",
    )

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work,
        inp=inp,
        out_prefix=out_prefix,
        extra_args=["--dump-accept", str(accept_log)],
        timeout=90,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    records = _load_accept_records(accept_log)
    assert records, "acceptance log unexpectedly empty"

    expected_by_species: dict[str, float] = {}
    for species, target_mean, target_z in (
        ("X", target_mean_x, target_z_x),
        ("Y", target_mean_y, target_z_y),
    ):
        first = next(
            r for r in records
            if str(r.get("species", "")).strip().upper() == species
        )
        assert float(first["vBox"]) == pytest.approx(v_box_expected, rel=1e-12, abs=1e-12)
        assert float(first["z"]) == pytest.approx(target_z, rel=2e-3, abs=1e-3)
        expected_by_species[species] = target_z * v_box_expected

    # Joint distribution should factorize for independent ideal-gas species (covariance ~ 0).
    counts = {"X": 0, "Y": 0}
    series_x: list[int] = []
    series_y: list[int] = []
    for rec in records:
        move = str(rec.get("move", "")).strip().lower()
        species = str(rec.get("species", "")).strip().upper()
        if move in ("insertion", "deletion") and species in counts:
            assert int(rec.get("nBefore", -1)) == counts[species]
        if move == "insertion" and species in counts and bool(rec.get("accepted")):
            counts[species] += 1
        elif move == "deletion" and species in counts and bool(rec.get("accepted")):
            counts[species] -= 1
        series_x.append(counts["X"])
        series_y.append(counts["Y"])

    burnin = max(300, mcsteps // 10)
    series_x = series_x[burnin:] if len(series_x) > burnin else series_x
    series_y = series_y[burnin:] if len(series_y) > burnin else series_y
    assert len(series_x) == len(series_y)
    assert len(series_x) >= 500

    mean_x = statistics.mean(series_x)
    mean_y = statistics.mean(series_y)
    var_x = statistics.pvariance(series_x)
    var_y = statistics.pvariance(series_y)

    for species, series, sample_mean, sample_var in (
        ("X", series_x, mean_x, var_x),
        ("Y", series_y, mean_y, var_y),
    ):
        expected_mean = float(expected_by_species[species])
        assert sample_mean == pytest.approx(expected_mean, rel=0.12, abs=0.25)
        assert abs(sample_var - expected_mean) / max(expected_mean, 1e-6) < 0.18

        hist = Counter(series)
        center = int(round(expected_mean))
        n_min = max(0, center - 2)
        n_max = center + 2
        min_bin_count = 30

        checked = 0
        for n in range(n_min, n_max + 1):
            c0 = hist.get(n, 0)
            c1 = hist.get(n + 1, 0)
            if c0 < min_bin_count or c1 < min_bin_count:
                continue
            empirical = c1 / c0
            expected = expected_mean / float(n + 1)
            assert empirical == pytest.approx(expected, rel=0.20, abs=0.10)
            checked += 1

        assert checked >= 3, f"Insufficient populated bins for {species}: checked={checked}, hist={hist}"

    cov_xy = statistics.mean(
        (x - mean_x) * (y - mean_y) for x, y in zip(series_x, series_y)
    )
    corr = cov_xy / math.sqrt(max(var_x, 1e-12) * max(var_y, 1e-12))
    assert abs(corr) < 0.15, f"Unexpected X/Y correlation in ideal gas: corr={corr:.3f}"
