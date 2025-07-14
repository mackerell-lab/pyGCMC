#!/bin/bash
# 运行不稳定的测试脚本

cd /home/zhaomt/gcmc/test107/pygcmc_dev/build

echo "运行不稳定的测试（串行模式）..."
echo "这些测试在并行执行时会导致段错误"
echo "================================================"

export PYTHONPATH=$PYTHONPATH:./modules/bindings
export PYTHONDONTWRITEBYTECODE=1

# 不稳定的测试列表
tests=(
    "../tests/simulation/energyPGP/debug_movement_residues.py::test_movement_residues"
    "../tests/simulation/energyPGP/debug_vdw_movement.py::test_vdw_movement_debug"
    "../tests/simulation/energyPGP/pgp_pme_debug_electrostatic.py::test_pgp_self_consistency"
)

passed=0
failed=0

for test in "${tests[@]}"; do
    echo -e "\n运行: $test"
    if PYTEST_RUN_UNSTABLE_TESTS=true ~/.miniconda3/envs/gcmc/bin/pytest "$test" -v --tb=short; then
        ((passed++))
        echo "✓ 通过"
    else
        ((failed++))
        echo "✗ 失败"
    fi
done

echo -e "\n================================================"
echo "结果: $passed 通过, $failed 失败"

if [ $failed -gt 0 ]; then
    echo -e "\n注意: 这些测试失败是由于C++代码中的内存管理问题"
    echo "需要在C++层面修复全局状态管理"
    exit 1
fi