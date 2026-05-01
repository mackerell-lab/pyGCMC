# CDK2 PGP/PME Example

This example runs the CDK2 manuscript benchmark system (56,470 atoms; 14,193 waters; 40.269 × 61.062 × 61.817 Å box) through the paper energy modes C/D/E.

## Command

```bash
./examples/build_gcmc_cpu.sh
./examples/cdk2_pgp_pme/run.sh
```

Outputs are written to `examples/cdk2_pgp_pme/runs/`.

This uses the CMake-built C++ executable directly. It does not require `setup.py` or pybind imports.
Input coordinates, topology, and force-field files are read from `data/paper_pgp_pme/`.

## Expected Result

`paper_summary.log` should show:

- `C-D` and `E-D` paired trial-energy differences against the slow raw PME reference path.
- `E-C` as the mesh-self-restoration term inside the PGP implementation for the same trial.
- PGP modes much faster than standard PME. Mode D is intentionally slow for this full CDK2 system; use `EXAMPLE_MODES=pgp_full,pgp_full_pme` for a fast PGP-only smoke run.
