"""
End-to-end tests for GROMACS ITP nonbonded parsing -> MCState LJ matrix wiring.

These tests provide an independent ("truth") check by reconstructing the expected
Lennard-Jones energy from:
  - the written final PDB coordinates (Å), and
  - known sigma/epsilon values from a minimal .itp PAR file (nm / kJ/mol),
then comparing against the structured acceptance log deltaU (kJ/mol).

This catches regressions where the simulation silently falls back to a placeholder/zero LJ matrix.
"""

from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import pytest


def _write_text(path: Path, content: str) -> None:
    path.write_text(content.strip() + "\n")


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


def _first_atom_xyz_angstrom(pdb_path: Path, *, resname: str) -> tuple[float, float, float]:
    want = resname.strip().upper()
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[17:20].strip().upper() != want:
            continue
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        return x, y, z
    raise AssertionError(f"Residue {resname} not found in {pdb_path}")


def _lj_energy_kj_mol(*, r_nm: float, sigma_nm: float, eps_kj_mol: float) -> float:
    r2 = r_nm * r_nm
    sigma_r2 = (sigma_nm * sigma_nm) / r2
    sigma_r6 = sigma_r2 * sigma_r2 * sigma_r2
    sigma_r12 = sigma_r6 * sigma_r6
    return 4.0 * eps_kj_mol * (sigma_r12 - sigma_r6)


def _c6_c12_from_sigma_eps(*, sigma_nm: float, eps_kj_mol: float) -> tuple[float, float]:
    c6 = 4.0 * eps_kj_mol * (sigma_nm**6)
    c12 = 4.0 * eps_kj_mol * (sigma_nm**12)
    return c6, c12


def _sigma_eps_from_c6_c12(*, c6: float, c12: float) -> tuple[float, float]:
    if c6 == 0.0 or c12 == 0.0:
        return 0.0, 0.0
    sigma = (c12 / c6) ** (1.0 / 6.0)
    eps = (c6 * c6) / (4.0 * c12)
    return sigma, eps


def test_itp_atomtypes_builds_lj_matrix_and_deltaU_matches_analytic_lj(
    gcmc_cpu, temp_dir
):
    """
    Independent energy check:
    - Initial system: one neutral atom type C at box center
    - Insert: one neutral atom type NA within a region near the initial atom
    - PAR: defines sigma/epsilon for both types (nm/kJ)
    Expect: acceptance log deltaU equals analytic LJ(C-NA) computed from final PDB coords.
    """
    work = Path(temp_dir) / "itp_nonbonded_energy" / "atomtypes_only"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
ATOM      1  C   MOL A   1      15.000  15.000  15.000  1.00  0.00           C
END
""",
    )

    top = work / "sys.top"
    _write_text(
        top,
        """
[ defaults ]
1 2 yes 0.5 0.8333

[ moleculetype ]
MOL  2

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge    mass
1   C     1      MOL      C     1     0.000   12.011

[ system ]
Minimal

[ molecules ]
MOL  1
""",
    )

    frag_itp = work / "na.itp"
    _write_text(
        frag_itp,
        """
[ moleculetype ]
NA  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   NA    1      NA       NA    1     0.000   22.9898
""",
    )

    # Sigma/epsilon are nm / kJ/mol (GROMACS convention)
    sigma_c = 0.300
    eps_c = 1.000
    sigma_na = 0.400
    eps_na = 2.000

    par = work / "ffnonbonded.itp"
    _write_text(
        par,
        f"""
[ atomtypes ]
; name  at.num  mass     charge   ptype    sigma      epsilon
C       6       12.011   0.000    A        {sigma_c:.6f}   {eps_c:.6f}
NA      11      22.990   0.000    A        {sigma_na:.6f}  {eps_na:.6f}
""",
    )

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "run.inp"
    _write_text(
        inp,
        f"""
random_seed:123
par:{par}
top:{top}
pdb:{pdb}
fragitp:{frag_itp}

box_size:30.0 30.0 30.0
cutoff:12.0
gcmc_region:box 20.0 12.5 12.5 25.0 17.5 17.5

temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1

fragname:NA
fragconc:0.01
fragmuex:10.0
mc_move_prob:1 0 0 0
""",
    )

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--dump-accept", str(accept_log)],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="insertion", species="NA")
    assert bool(rec.get("accepted")) is True

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()

    x0, y0, z0 = _first_atom_xyz_angstrom(final_pdb, resname="MOL")
    x1, y1, z1 = _first_atom_xyz_angstrom(final_pdb, resname="NA")

    dx = (x1 - x0) / 10.0
    dy = (y1 - y0) / 10.0
    dz = (z1 - z0) / 10.0
    r_nm = math.sqrt(dx * dx + dy * dy + dz * dz)

    sigma_mixed = 0.5 * (sigma_c + sigma_na)
    eps_mixed = math.sqrt(eps_c * eps_na)
    expected = _lj_energy_kj_mol(r_nm=r_nm, sigma_nm=sigma_mixed, eps_kj_mol=eps_mixed)

    assert float(rec["deltaU"]) == pytest.approx(expected, rel=5e-4, abs=5e-4)


def test_itp_nonbond_params_override_applies_nbfix_in_lj_matrix(
    gcmc_cpu, temp_dir
):
    """
    End-to-end check that [ nonbond_params ] entries are applied as NBFIX overrides.

    We set an explicit C-NA override (sigma/epsilon) and verify deltaU matches the
    override parameters (not the Lorentz-Berthelot mixed parameters).
    """
    work = Path(temp_dir) / "itp_nonbonded_energy" / "nonbond_params_override"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
ATOM      1  C   MOL A   1      15.000  15.000  15.000  1.00  0.00           C
END
""",
    )

    top = work / "sys.top"
    _write_text(
        top,
        """
[ defaults ]
1 2 yes 0.5 0.8333

[ moleculetype ]
MOL  2

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge    mass
1   C     1      MOL      C     1     0.000   12.011

[ system ]
Minimal

[ molecules ]
MOL  1
""",
    )

    frag_itp = work / "na.itp"
    _write_text(
        frag_itp,
        """
[ moleculetype ]
NA  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   NA    1      NA       NA    1     0.000   22.9898
""",
    )

    sigma_c = 0.300
    eps_c = 1.000
    sigma_na = 0.400
    eps_na = 2.000

    sigma_override = 0.400
    eps_override = 20.000

    par = work / "ffnonbonded.itp"
    _write_text(
        par,
        f"""
[ atomtypes ]
; name  at.num  mass     charge   ptype    sigma      epsilon
C       6       12.011   0.000    A        {sigma_c:.6f}   {eps_c:.6f}
NA      11      22.990   0.000    A        {sigma_na:.6f}  {eps_na:.6f}

[ nonbond_params ]
; type1  type2  func  sigma  epsilon
C   NA   1   {sigma_override:.6f}  {eps_override:.6f}
""",
    )

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "run.inp"
    _write_text(
        inp,
        f"""
random_seed:123
par:{par}
top:{top}
pdb:{pdb}
fragitp:{frag_itp}

box_size:30.0 30.0 30.0
cutoff:12.0
gcmc_region:box 20.0 12.5 12.5 25.0 17.5 17.5

temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1

fragname:NA
fragconc:0.01
fragmuex:10.0
mc_move_prob:1 0 0 0
""",
    )

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--dump-accept", str(accept_log)],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="insertion", species="NA")
    assert bool(rec.get("accepted")) is True

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()

    x0, y0, z0 = _first_atom_xyz_angstrom(final_pdb, resname="MOL")
    x1, y1, z1 = _first_atom_xyz_angstrom(final_pdb, resname="NA")

    dx = (x1 - x0) / 10.0
    dy = (y1 - y0) / 10.0
    dz = (z1 - z0) / 10.0
    r_nm = math.sqrt(dx * dx + dy * dy + dz * dz)

    expected = _lj_energy_kj_mol(r_nm=r_nm, sigma_nm=sigma_override, eps_kj_mol=eps_override)
    assert float(rec["deltaU"]) == pytest.approx(expected, rel=5e-4, abs=5e-4)


def test_itp_pairtypes_override_applies_nbfix_in_lj_matrix(gcmc_cpu, temp_dir):
    """
    End-to-end check that [ pairtypes ] entries are applied as LJ overrides.

    This mirrors the gcmc_gpu parser behavior where both [pairtypes] and [nonbond_params]
    populate the override map used for LJ.
    """
    work = Path(temp_dir) / "itp_nonbonded_energy" / "pairtypes_override"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
ATOM      1  C   MOL A   1      15.000  15.000  15.000  1.00  0.00           C
END
""",
    )

    top = work / "sys.top"
    _write_text(
        top,
        """
[ defaults ]
1 2 yes 0.5 0.8333

[ moleculetype ]
MOL  2

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge    mass
1   C     1      MOL      C     1     0.000   12.011

[ system ]
Minimal

[ molecules ]
MOL  1
""",
    )

    frag_itp = work / "na.itp"
    _write_text(
        frag_itp,
        """
[ moleculetype ]
NA  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   NA    1      NA       NA    1     0.000   22.9898
""",
    )

    sigma_c = 0.300
    eps_c = 1.000
    sigma_na = 0.400
    eps_na = 2.000

    sigma_override = 0.410
    eps_override = 20.000

    par = work / "ffnonbonded.itp"
    _write_text(
        par,
        f"""
[ atomtypes ]
; name  at.num  mass     charge   ptype    sigma      epsilon
C       6       12.011   0.000    A        {sigma_c:.6f}   {eps_c:.6f}
NA      11      22.990   0.000    A        {sigma_na:.6f}  {eps_na:.6f}

[ pairtypes ]
; type1  type2  func  sigma  epsilon
C   NA   1   {sigma_override:.6f}  {eps_override:.6f}
""",
    )

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "run.inp"
    _write_text(
        inp,
        f"""
random_seed:123
par:{par}
top:{top}
pdb:{pdb}
fragitp:{frag_itp}

box_size:30.0 30.0 30.0
cutoff:12.0
gcmc_region:box 20.0 12.5 12.5 25.0 17.5 17.5

temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1

fragname:NA
fragconc:0.01
fragmuex:10.0
mc_move_prob:1 0 0 0
""",
    )

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--dump-accept", str(accept_log)],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="insertion", species="NA")
    assert bool(rec.get("accepted")) is True

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()

    x0, y0, z0 = _first_atom_xyz_angstrom(final_pdb, resname="MOL")
    x1, y1, z1 = _first_atom_xyz_angstrom(final_pdb, resname="NA")

    dx = (x1 - x0) / 10.0
    dy = (y1 - y0) / 10.0
    dz = (z1 - z0) / 10.0
    r_nm = math.sqrt(dx * dx + dy * dy + dz * dz)

    expected = _lj_energy_kj_mol(r_nm=r_nm, sigma_nm=sigma_override, eps_kj_mol=eps_override)
    assert float(rec["deltaU"]) == pytest.approx(expected, rel=5e-4, abs=5e-4)


def test_itp_pairtypes_strict_mode_does_not_override_inter_residue(gcmc_cpu, temp_dir):
    """
    Strict GROMACS mode: [ pairtypes ] must not override inter-residue LJ.

    This regression ensures pairtypes are reserved for 1-4 interactions only.
    """
    work = Path(temp_dir) / "itp_nonbonded_energy" / "pairtypes_strict_mode"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
ATOM      1  C   MOL A   1      15.000  15.000  15.000  1.00  0.00           C
END
""",
    )

    top = work / "sys.top"
    _write_text(
        top,
        """
[ defaults ]
1 2 yes 0.5 0.8333

[ moleculetype ]
MOL  2

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge    mass
1   C     1      MOL      C     1     0.000   12.011

[ system ]
Minimal

[ molecules ]
MOL  1
""",
    )

    frag_itp = work / "na.itp"
    _write_text(
        frag_itp,
        """
[ moleculetype ]
NA  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   NA    1      NA       NA    1     0.000   22.9898
""",
    )

    sigma_c = 0.300
    eps_c = 1.000
    sigma_na = 0.400
    eps_na = 2.000

    sigma_override = 0.410
    eps_override = 50.000

    par = work / "ffnonbonded.itp"
    _write_text(
        par,
        f"""
[ atomtypes ]
; name  at.num  mass     charge   ptype    sigma      epsilon
C       6       12.011   0.000    A        {sigma_c:.6f}   {eps_c:.6f}
NA      11      22.990   0.000    A        {sigma_na:.6f}  {eps_na:.6f}

[ pairtypes ]
; type1  type2  func  sigma  epsilon
C   NA   1   {sigma_override:.6f}  {eps_override:.6f}
""",
    )

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "run.inp"
    _write_text(
        inp,
        f"""
random_seed:123
par:{par}
top:{top}
pdb:{pdb}
fragitp:{frag_itp}
itp_pairtypes_mode:strict

box_size:30.0 30.0 30.0
cutoff:12.0
gcmc_region:box 20.0 12.5 12.5 25.0 17.5 17.5

temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1

fragname:NA
fragconc:0.01
fragmuex:10.0
mc_move_prob:1 0 0 0
""",
    )

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--dump-accept", str(accept_log)],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="insertion", species="NA")
    assert bool(rec.get("accepted")) is True

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()

    x0, y0, z0 = _first_atom_xyz_angstrom(final_pdb, resname="MOL")
    x1, y1, z1 = _first_atom_xyz_angstrom(final_pdb, resname="NA")

    dx = (x1 - x0) / 10.0
    dy = (y1 - y0) / 10.0
    dz = (z1 - z0) / 10.0
    r_nm = math.sqrt(dx * dx + dy * dy + dz * dz)

    sigma_mixed = 0.5 * (sigma_c + sigma_na)
    eps_mixed = math.sqrt(eps_c * eps_na)
    expected = _lj_energy_kj_mol(r_nm=r_nm, sigma_nm=sigma_mixed, eps_kj_mol=eps_mixed)

    assert float(rec["deltaU"]) == pytest.approx(expected, rel=5e-4, abs=5e-4)


def test_dump_params_reports_itp_defaults_comb_rule(gcmc_cpu, temp_dir):
    """
    Ensure [defaults] from a separate ITP file is visible in --dump-params.
    """
    work = Path(temp_dir) / "itp_nonbonded_energy" / "defaults_dump_params"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
ATOM      1  C   MOL A   1      15.000  15.000  15.000  1.00  0.00           C
END
""",
    )

    top = work / "sys.top"
    _write_text(
        top,
        """
[ defaults ]
1 2 yes 0.5 0.8333

[ moleculetype ]
MOL  2

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge    mass
1   C     1      MOL      C     1     0.000   12.011

[ system ]
Minimal

[ molecules ]
MOL  1
""",
    )

    frag_itp = work / "na.itp"
    _write_text(
        frag_itp,
        """
[ moleculetype ]
NA  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   NA    1      NA       NA    1     0.000   22.9898
""",
    )

    defaults_itp = work / "forcefield.itp"
    _write_text(
        defaults_itp,
        """
[ defaults ]
1 3 yes 1.0 1.0
""",
    )

    par = work / "ffnonbonded.itp"
    _write_text(
        par,
        """
[ atomtypes ]
; name  at.num  mass     charge   ptype    sigma      epsilon
C       6       12.011   0.000    A        0.300000   1.000000
NA      11      22.990   0.000    A        0.400000   2.000000
""",
    )

    dump_params = work / "params.json"
    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    inp = work / "run.inp"
    _write_text(
        inp,
        f"""
random_seed:123
par:{par}
itp_defaults:{defaults_itp}
top:{top}
pdb:{pdb}
fragitp:{frag_itp}

box_size:30.0 30.0 30.0
cutoff:12.0
gcmc_region:box 20.0 12.5 12.5 25.0 17.5 17.5

temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1

fragname:NA
fragconc:0.01
fragmuex:10.0
mc_move_prob:1 0 0 0
""",
    )

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--dump-params", str(dump_params)],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    params = json.loads(dump_params.read_text())
    assert params["basic"]["gromacs_defaults_present"] is True
    assert params["basic"]["gromacs_nbfunc"] == 1
    assert params["basic"]["gromacs_comb_rule"] == 3


def test_itp_defaults_gen_pairs_and_fudge_are_visible_and_strict_fails(gcmc_cpu, temp_dir):
    """
    [defaults] gen-pairs/fudge values must be visible and trigger strict failure when unsupported.
    """
    work = Path(temp_dir) / "itp_nonbonded_energy" / "defaults_gen_pairs_fudge"
    work.mkdir(parents=True, exist_ok=True)

    frag_itp = work / "frag.itp"
    _write_text(
        frag_itp,
        """
[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   X     1      FRG      X     1     0.000   1.000
""",
    )

    defaults_itp = work / "forcefield.itp"
    _write_text(
        defaults_itp,
        """
[ defaults ]
1 2 no 0.5 0.5
""",
    )

    par = work / "ffnonbonded.itp"
    _write_text(
        par,
        """
[ atomtypes ]
; name  at.num  mass     charge   ptype    sigma      epsilon
X       0       1.000    0.000    A        0.300000   0.000000
""",
    )

    inp = work / "run.inp"
    _write_text(
        inp,
        f"""
random_seed:123
par:{par}
itp_defaults:{defaults_itp}
fragitp:{frag_itp}

fragname:FRG
fragconc:1.0
fragmuex:0.0

box_size:10.0 10.0 10.0
cutoff:8.0
temperature:300.0
mcsteps:0
nprint:1
mc_move_prob:1 0 0 0
""",
    )

    dump_params = work / "params.json"
    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--dump-params", str(dump_params)],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    params = json.loads(dump_params.read_text())
    assert params["basic"]["gromacs_defaults_present"] is True
    assert params["basic"]["gromacs_gen_pairs_present"] is True
    assert params["basic"]["gromacs_gen_pairs"].lower() == "no"
    assert params["basic"]["gromacs_fudge_present"] is True
    assert float(params["basic"]["gromacs_fudge_lj"]) == pytest.approx(0.5, abs=1e-6)
    assert float(params["basic"]["gromacs_fudge_qq"]) == pytest.approx(0.5, abs=1e-6)

    ignored = set(params["basic"]["ignored_inp_keys"])
    assert "gromacs_gen_pairs" in ignored
    assert "gromacs_fudge_lj" in ignored
    assert "gromacs_fudge_qq" in ignored

    dump_params_strict = work / "params_strict.json"
    result_strict = subprocess.run(
        [
            gcmc_cpu,
            "--inp",
            str(inp),
            "--prefix",
            str(out_prefix),
            "--dump-params",
            str(dump_params_strict),
            "--strict-inp-keys",
        ],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result_strict.returncode != 0


def test_itp_nonbond_params_take_precedence_over_pairtypes_for_same_pair(gcmc_cpu, temp_dir):
    """
    Regression for override precedence.

    When BOTH [ nonbond_params ] and [ pairtypes ] specify the same type pair, we treat
    nonbond_params as an NBFIX-like override and it must take precedence.

    The fixture places [pairtypes] *after* [nonbond_params] so a naïve "last wins"
    parser would be wrong.
    """
    work = Path(temp_dir) / "itp_nonbonded_energy" / "override_precedence"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
ATOM      1  C   MOL A   1      15.000  15.000  15.000  1.00  0.00           C
END
""",
    )

    top = work / "sys.top"
    _write_text(
        top,
        """
[ defaults ]
1 2 yes 0.5 0.8333

[ moleculetype ]
MOL  2

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge    mass
1   C     1      MOL      C     1     0.000   12.011

[ system ]
Minimal

[ molecules ]
MOL  1
""",
    )

    frag_itp = work / "na.itp"
    _write_text(
        frag_itp,
        """
[ moleculetype ]
NA  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   NA    1      NA       NA    1     0.000   22.9898
""",
    )

    sigma_c = 0.300
    eps_c = 1.000
    sigma_na = 0.400
    eps_na = 2.000

    sigma_nbfix = 0.390
    eps_nbfix = 30.000
    sigma_pair = 0.410
    eps_pair = 20.000

    par = work / "ffnonbonded.itp"
    _write_text(
        par,
        f"""
[ atomtypes ]
; name  at.num  mass     charge   ptype    sigma      epsilon
C       6       12.011   0.000    A        {sigma_c:.6f}   {eps_c:.6f}
NA      11      22.990   0.000    A        {sigma_na:.6f}  {eps_na:.6f}

[ nonbond_params ]
; type1  type2  func  sigma  epsilon
C   NA   1   {sigma_nbfix:.6f}  {eps_nbfix:.6f}

[ pairtypes ]
; type1  type2  func  sigma  epsilon
C   NA   1   {sigma_pair:.6f}  {eps_pair:.6f}
""",
    )

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "run.inp"
    _write_text(
        inp,
        f"""
random_seed:123
par:{par}
top:{top}
pdb:{pdb}
fragitp:{frag_itp}

box_size:30.0 30.0 30.0
cutoff:12.0
gcmc_region:box 20.0 12.5 12.5 25.0 17.5 17.5

temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1

fragname:NA
fragconc:0.01
fragmuex:10.0
mc_move_prob:1 0 0 0
""",
    )

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--dump-accept", str(accept_log)],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="insertion", species="NA")
    assert bool(rec.get("accepted")) is True

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()

    x0, y0, z0 = _first_atom_xyz_angstrom(final_pdb, resname="MOL")
    x1, y1, z1 = _first_atom_xyz_angstrom(final_pdb, resname="NA")

    dx = (x1 - x0) / 10.0
    dy = (y1 - y0) / 10.0
    dz = (z1 - z0) / 10.0
    r_nm = math.sqrt(dx * dx + dy * dy + dz * dz)

    expected = _lj_energy_kj_mol(r_nm=r_nm, sigma_nm=sigma_nbfix, eps_kj_mol=eps_nbfix)
    assert float(rec["deltaU"]) == pytest.approx(expected, rel=5e-4, abs=5e-4)


def test_itp_defaults_comb_rule_geometric_sigma_mixing(gcmc_cpu, temp_dir):
    """
    Regression for [defaults] comb-rule=3 (geometric sigma).

    With comb-rule 3, sigma and epsilon should both be mixed geometrically.
    """
    work = Path(temp_dir) / "itp_nonbonded_energy" / "comb_rule_geom_sigma"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
ATOM      1  C   MOL A   1      15.000  15.000  15.000  1.00  0.00           C
END
""",
    )

    top = work / "sys.top"
    _write_text(
        top,
        """
[ defaults ]
1 2 yes 0.5 0.8333

[ moleculetype ]
MOL  2

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge    mass
1   C     1      MOL      C     1     0.000   12.011

[ system ]
Minimal

[ molecules ]
MOL  1
""",
    )

    frag_itp = work / "na.itp"
    _write_text(
        frag_itp,
        """
[ moleculetype ]
NA  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   NA    1      NA       NA    1     0.000   22.9898
""",
    )

    sigma_c = 0.200
    eps_c = 1.200
    sigma_na = 0.600
    eps_na = 2.500

    par = work / "ffnonbonded.itp"
    _write_text(
        par,
        f"""
[ defaults ]
1 3 yes 1.0 1.0

[ atomtypes ]
; name  at.num  mass     charge   ptype    sigma      epsilon
C       6       12.011   0.000    A        {sigma_c:.6f}   {eps_c:.6f}
NA      11      22.990   0.000    A        {sigma_na:.6f}  {eps_na:.6f}
""",
    )

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "run.inp"
    _write_text(
        inp,
        f"""
random_seed:123
par:{par}
top:{top}
pdb:{pdb}
fragitp:{frag_itp}

box_size:30.0 30.0 30.0
cutoff:12.0
gcmc_region:box 20.0 12.5 12.5 25.0 17.5 17.5

temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1

fragname:NA
fragconc:0.01
fragmuex:30.0
mc_move_prob:1 0 0 0
""",
    )

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--dump-accept", str(accept_log)],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="insertion", species="NA")
    assert bool(rec.get("accepted")) is True

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()

    x0, y0, z0 = _first_atom_xyz_angstrom(final_pdb, resname="MOL")
    x1, y1, z1 = _first_atom_xyz_angstrom(final_pdb, resname="NA")

    dx = (x1 - x0) / 10.0
    dy = (y1 - y0) / 10.0
    dz = (z1 - z0) / 10.0
    r_nm = math.sqrt(dx * dx + dy * dy + dz * dz)

    sigma_mixed = math.sqrt(sigma_c * sigma_na)
    eps_mixed = math.sqrt(eps_c * eps_na)
    expected = _lj_energy_kj_mol(r_nm=r_nm, sigma_nm=sigma_mixed, eps_kj_mol=eps_mixed)

    assert float(rec["deltaU"]) == pytest.approx(expected, rel=5e-4, abs=5e-4)


def test_itp_defaults_comb_rule_one_converts_c6_c12_and_overrides(gcmc_cpu, temp_dir):
    """
    Regression for [defaults] comb-rule=1 (C6/C12 -> sigma/epsilon).

    Ensures both atomtypes and nonbond_params are converted from C6/C12
    into sigma/epsilon before building the LJ matrix.
    """
    work = Path(temp_dir) / "itp_nonbonded_energy" / "comb_rule_c6_c12"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
ATOM      1  C   MOL A   1      15.000  15.000  15.000  1.00  0.00           C
END
""",
    )

    top = work / "sys.top"
    _write_text(
        top,
        """
[ defaults ]
1 2 yes 0.5 0.8333

[ moleculetype ]
MOL  2

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge    mass
1   C     1      MOL      C     1     0.000   12.011

[ system ]
Minimal

[ molecules ]
MOL  1
""",
    )

    frag_itp = work / "na.itp"
    _write_text(
        frag_itp,
        """
[ moleculetype ]
NA  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   NA    1      NA       NA    1     0.000   22.9898
""",
    )

    sigma_c = 0.280
    eps_c = 1.100
    sigma_na = 0.450
    eps_na = 2.200
    sigma_override = 0.520
    eps_override = 3.700

    c6_c, c12_c = _c6_c12_from_sigma_eps(sigma_nm=sigma_c, eps_kj_mol=eps_c)
    c6_na, c12_na = _c6_c12_from_sigma_eps(sigma_nm=sigma_na, eps_kj_mol=eps_na)
    c6_override, c12_override = _c6_c12_from_sigma_eps(
        sigma_nm=sigma_override,
        eps_kj_mol=eps_override,
    )

    par = work / "ffnonbonded.itp"
    _write_text(
        par,
        f"""
[ defaults ]
1 1 yes 1.0 1.0

[ atomtypes ]
; name  at.num  mass     charge   ptype    C6        C12
C       6       12.011   0.000    A        {c6_c:.8e}   {c12_c:.8e}
NA      11      22.990   0.000    A        {c6_na:.8e}  {c12_na:.8e}

[ nonbond_params ]
; type1  type2  func  C6  C12
C   NA   1   {c6_override:.8e}  {c12_override:.8e}
""",
    )

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "run.inp"
    _write_text(
        inp,
        f"""
random_seed:123
par:{par}
top:{top}
pdb:{pdb}
fragitp:{frag_itp}

box_size:30.0 30.0 30.0
cutoff:12.0
gcmc_region:box 20.0 12.5 12.5 25.0 17.5 17.5

temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1

fragname:NA
fragconc:0.01
fragmuex:30.0
mc_move_prob:1 0 0 0
""",
    )

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--dump-accept", str(accept_log)],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="insertion", species="NA")
    assert bool(rec.get("accepted")) is True

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()

    x0, y0, z0 = _first_atom_xyz_angstrom(final_pdb, resname="MOL")
    x1, y1, z1 = _first_atom_xyz_angstrom(final_pdb, resname="NA")

    dx = (x1 - x0) / 10.0
    dy = (y1 - y0) / 10.0
    dz = (z1 - z0) / 10.0
    r_nm = math.sqrt(dx * dx + dy * dy + dz * dz)

    expected = _lj_energy_kj_mol(
        r_nm=r_nm,
        sigma_nm=sigma_override,
        eps_kj_mol=eps_override,
    )
    assert float(rec["deltaU"]) == pytest.approx(expected, rel=5e-4, abs=5e-4)
