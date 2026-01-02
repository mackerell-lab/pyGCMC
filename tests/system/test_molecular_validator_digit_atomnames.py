from __future__ import annotations

from pathlib import Path


def test_molecular_combine_accepts_pdb_atom_names_starting_with_digit(tmp_path: Path):
    """
    Regression: PDB atom names like '1HD2' must still be treated as element 'H' when
    validating structure vs topology, otherwise protein systems silently fall back to
    INP-only initialization.
    """
    import pygcmc

    pdb = tmp_path / "sys.pdb"
    pdb.write_text(
        "\n".join(
            [
                "CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1",
                # Atom name starts with digit; element field intentionally left empty.
                "ATOM      1 1HD2 MOL A   1       0.000   0.000   0.000  1.00  0.00",
                "END",
                "",
            ]
        )
    )

    top = tmp_path / "sys.top"
    top.write_text(
        "\n".join(
            [
                "[ defaults ]",
                "1 2 yes 1.0 1.0",
                "",
                "[ moleculetype ]",
                "MOL  1",
                "",
                "[ atoms ]",
                "; nr  type  resnr  residue  atom  cgnr  charge  mass",
                "1   H     1      MOL      1HD2   1     0.000   1.008",
                "",
                "[ system ]",
                "Minimal",
                "",
                "[ molecules ]",
                "MOL  1",
                "",
            ]
        )
    )

    structure = pygcmc.PDBParser.parse_file(str(pdb))
    topology = pygcmc.TOPParser.parse_file(str(top))
    molecular = pygcmc.MolecularSystem().combine(structure, topology)
    assert molecular is not None

