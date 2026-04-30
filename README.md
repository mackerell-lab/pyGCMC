# pygcmc

## Minimal Build

From the workspace root, build the standalone C++ runner:

```bash
cmake -S pygcmc_dev -B pygcmc_dev/build -DCMAKE_BUILD_TYPE=Release
cmake --build pygcmc_dev/build --target gcmc_cpu -j4
```

## Run One Example

Run one paper example from the workspace-level `example/` directory:

```bash
./example/t4_lysozyme_pgp_pme/run.sh
cat example/t4_lysozyme_pgp_pme/runs/paper_summary.log
```

The script uses `gcmc_cpu` and writes inputs, logs, timing data, `summary.json`, and `paper_summary.log` under the example's `runs/` directory.
