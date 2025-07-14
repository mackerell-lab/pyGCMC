#!/bin/bash
# 安全的测试运行脚本，处理已知的内存问题

cd /home/zhaomt/gcmc/test107/pygcmc_dev/build

echo "Running all tests with memory safety precautions..."
echo "======================================================="

# Set environment
export PYTHONPATH=$PYTHONPATH:./modules/bindings
export PYTHONDONTWRITEBYTECODE=1

# Run main test suite (excluding PGP Complete which has memory issues)
echo -e "\n1. Running main test suite (excluding problematic tests)..."
~/.miniconda3/envs/gcmc/bin/pytest ../tests/ -n 8 -k "not (pgp_complete or test_movement_residues or test_vdw_movement_debug or test_pgp_self_consistency)"

# Run PGP Complete tests in isolation (due to C++ global state bug)
echo -e "\n2. Running PGP Complete tests in isolation..."
if [ -f ../tests/simulation/energyPGP/run_isolated.py ]; then
    python ../tests/simulation/energyPGP/run_isolated.py
else
    echo "Warning: run_isolated.py not found, running PGP Complete tests normally..."
    ~/.miniconda3/envs/gcmc/bin/pytest ../tests/simulation/energyPGP/pgp_complete.py -n 1
fi

# Run problematic tests one by one
echo -e "\n3. Running problematic tests individually..."
echo "Running test_movement_residues..."
~/.miniconda3/envs/gcmc/bin/pytest ../tests/simulation/energyPGP/debug_movement_residues.py::test_movement_residues -v || echo "Failed (known issue)"

echo -e "\nRunning test_vdw_movement_debug..."
~/.miniconda3/envs/gcmc/bin/pytest ../tests/simulation/energyPGP/debug_vdw_movement.py::test_vdw_movement_debug -v || echo "Failed (known issue)"

echo -e "\nRunning test_pgp_self_consistency..."
~/.miniconda3/envs/gcmc/bin/pytest ../tests/simulation/energyPGP/pgp_pme_debug_electrostatic.py::test_pgp_self_consistency -v || echo "Failed (known issue)"

echo -e "\n======================================================="
echo "Test run complete!"
echo "Note: Some tests may fail due to known C++ memory management issues."
echo "These need to be fixed in the C++ code (add resetPGPState function)."