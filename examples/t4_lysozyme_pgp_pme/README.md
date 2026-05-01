# T4 Lysozyme PGP/PME Example

This example runs the T4 lysozyme manuscript benchmark system (37,871 atoms; 10,436 waters; 36.736 × 40.850 × 49.379 Å box) through the paper energy modes C/D/E.

## Command

```bash
./examples/build_gcmc_cpu.sh
./examples/t4_lysozyme_pgp_pme/run.sh
```

Outputs are written to `examples/t4_lysozyme_pgp_pme/runs/`.

This uses the CMake-built C++ executable directly. It does not require `setup.py` or pybind imports.
Input coordinates, topology, and force-field files are read from `data/paper_pgp_pme/`.

## Expected Result

`paper_summary.log` should show:

- `C-D` and `E-D` paired trial-energy differences against the slow raw PME reference path.
- `E-C` as the mesh-self-restoration term inside the PGP implementation for the same trial.
- PGP modes much faster than standard PME, consistent with the manuscript speed comparison for the T4 system.
