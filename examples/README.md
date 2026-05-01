# Paper PGP/PME Examples

These examples reproduce the two manuscript-scale PGP/PME comparison cases using the current `gcmc_cpu` executable directly.

The simulation path is CMake-only: it does not require `setup.py`, `PYTHONPATH`, or importing the pybind11 `pygcmc` module. Python is used only after the C++ runs to summarize JSONL logs.
The examples use standalone input assets under `data/paper_pgp_pme/`.

## Cases

- `t4_lysozyme_pgp_pme/`: T4 lysozyme benchmark system, 37,871 atoms.
- `cdk2_pgp_pme/`: CDK2 benchmark system, 56,470 atoms.

Both cases run the manuscript energy modes:

- `pgp_full` = Mode C, cross-only PGP with the PME mesh-self artifact excluded.
- `pme` = Mode D, standard PME reference.
- `pgp_full_pme` = Mode E, PGP with mesh-self restored to reproduce PME.

## Run

Build `gcmc_cpu` first from the repository root:

```bash
./examples/build_gcmc_cpu.sh
```

Equivalent explicit CMake commands:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target gcmc_cpu -j4
```

Run one case from the repository root:

```bash
./examples/t4_lysozyme_pgp_pme/run.sh
./examples/cdk2_pgp_pme/run.sh
```

The scripts write per-mode inputs, raw stdout/stderr, acceptance JSONL files, timing data, `summary.json`, and `paper_summary.log` under each case's `runs/` directory.

The scripts look for `gcmc_cpu` in this order: `GCMC_CPU`, `PATH`, then `build/bin/gcmc_cpu`.

## Notes

The default is `EXAMPLE_STEPS=1` because Mode D (`pme`) is intentionally slow on full manuscript systems. Increase `EXAMPLE_STEPS` for timing averages. To run only the two fast PGP modes:

```bash
EXAMPLE_MODES=pgp_full,pgp_full_pme ./examples/cdk2_pgp_pme/run.sh
```
