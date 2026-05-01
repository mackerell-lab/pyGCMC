# pygcmc

## Minimal Build

From the repository root, build the standalone C++ runner:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target gcmc_cpu -j4
```

## Run One Example

Run one paper example from the repository `examples/` directory:

```bash
./examples/t4_lysozyme_pgp_pme/run.sh
cat examples/t4_lysozyme_pgp_pme/runs/paper_summary.log
```

The script uses `gcmc_cpu` and writes inputs, logs, timing data, `summary.json`, and `paper_summary.log` under the example's `runs/` directory.
