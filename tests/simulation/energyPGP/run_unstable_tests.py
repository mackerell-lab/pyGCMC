#!/usr/bin/env python
"""
运行不稳定的测试，每个测试在独立的进程中
"""

import subprocess
import sys
import os

# 设置环境变量
env = os.environ.copy()
env['PYTHONPATH'] = env.get('PYTHONPATH', '') + ':./modules/bindings'
env['PYTHONDONTWRITEBYTECODE'] = '1'

# 不稳定的测试列表
unstable_tests = [
    'debug_movement_residues.py::test_movement_residues',
    'debug_vdw_movement.py::test_vdw_movement_debug',
    'pgp_pme_debug_electrostatic.py::test_pgp_self_consistency'
]

print("运行不稳定的测试（每个在独立进程中）...")
print("=" * 60)

passed = 0
failed = 0
failed_tests = []

for test in unstable_tests:
    print(f"\n运行: {test}")
    cmd = [
        sys.executable, '-m', 'pytest',
        f'../tests/simulation/energyPGP/{test}',
        '-v', '--tb=short', '--no-header'
    ]
    
    try:
        result = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=30)
        if result.returncode == 0:
            print(f"✓ 通过")
            passed += 1
        else:
            print(f"✗ 失败")
            failed += 1
            failed_tests.append(test)
            if result.stderr:
                print(f"错误: {result.stderr}")
    except subprocess.TimeoutExpired:
        print(f"✗ 超时")
        failed += 1
        failed_tests.append(test)

print("\n" + "=" * 60)
print(f"结果: {passed} 通过, {failed} 失败")
if failed_tests:
    print("\n失败的测试:")
    for test in failed_tests:
        print(f"  - {test}")
    sys.exit(1)
else:
    print("\n所有不稳定的测试都通过了！")
    sys.exit(0)