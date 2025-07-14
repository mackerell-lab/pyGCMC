#!/usr/bin/env python3
"""
运行不稳定测试的隔离脚本
每个测试在完全独立的进程中运行，避免内存污染
"""

import subprocess
import sys
import os
import time
import json
from pathlib import Path

class IsolatedTestRunner:
    def __init__(self):
        self.build_dir = '/home/zhaomt/gcmc/test107/pygcmc_dev/build'
        self.test_dir = '/home/zhaomt/gcmc/test107/pygcmc_dev/tests'
        
        # 不稳定的测试列表
        self.unstable_tests = [
            ('simulation/energyPGP/pgp_pme_debug_electrostatic.py', 'test_pgp_self_consistency'),
            ('simulation/energyPGP/debug_movement_residues.py', 'test_movement_residues'),
            ('simulation/energyPGP/debug_vdw_movement.py', 'test_vdw_movement_debug'),
            ('simulation/energyOpenmm/naive_energy_components.py', 'test_compare_openmm_naive_nonbonded'),
        ]
        
    def run_single_test(self, test_file, test_function):
        """在独立进程中运行单个测试"""
        test_path = os.path.join(self.test_dir, test_file)
        
        # 创建运行脚本
        runner_script = f"""
import sys
import os
import gc
import subprocess

# 设置环境
os.chdir('{self.build_dir}')
sys.path.insert(0, './modules/bindings')
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
os.environ['MALLOC_CHECK_'] = '0'
os.environ['MALLOC_PERTURB_'] = '0'

# 使用pytest运行以确保正确的导入行为
test_path = '{test_path}::{test_function}'

# 构建pytest命令
cmd = [
    sys.executable, '-m', 'pytest',
    test_path,
    '-v', '-s', '--tb=short',
    '--no-header'
]

# 运行测试
result = subprocess.run(cmd, capture_output=True, text=True)

print(result.stdout)
if result.stderr:
    print("STDERR:", result.stderr)

sys.exit(result.returncode)
"""
        
        # 使用subprocess运行，完全隔离
        env = os.environ.copy()
        env['PYTHONPATH'] = ''  # 清空PYTHONPATH，避免干扰
        
        try:
            result = subprocess.run(
                [sys.executable, '-c', runner_script],
                capture_output=True,
                text=True,
                timeout=60,
                cwd=self.build_dir,
                env=env
            )
            
            return {
                'success': result.returncode == 0,
                'stdout': result.stdout,
                'stderr': result.stderr,
                'returncode': result.returncode
            }
            
        except subprocess.TimeoutExpired:
            return {
                'success': False,
                'stdout': '',
                'stderr': f'测试超时 (60秒)',
                'returncode': -1
            }
        except Exception as e:
            return {
                'success': False,
                'stdout': '',
                'stderr': f'运行错误: {e}',
                'returncode': -2
            }
    
    def run_all_tests(self, num_iterations=1):
        """运行所有不稳定的测试"""
        print(f"隔离运行不稳定的测试 (迭代次数: {num_iterations})")
        print("="*70)
        
        all_results = []
        
        for iteration in range(num_iterations):
            if num_iterations > 1:
                print(f"\n\n迭代 {iteration + 1}/{num_iterations}")
                print("-"*70)
            
            iteration_results = []
            
            for test_file, test_function in self.unstable_tests:
                test_name = f"{test_file}::{test_function}"
                print(f"\n测试: {test_name}")
                
                # 在每个测试之间等待一小段时间，让系统清理资源
                time.sleep(0.5)
                
                result = self.run_single_test(test_file, test_function)
                
                iteration_results.append({
                    'test': test_name,
                    'success': result['success'],
                    'iteration': iteration + 1
                })
                
                if result['success']:
                    print("✅ 通过")
                else:
                    print("❌ 失败")
                    if result['stderr']:
                        print(f"错误输出:\n{result['stderr'][:500]}")
                
                # 如果需要，打印详细输出
                if not result['success'] and result['stdout']:
                    print(f"标准输出:\n{result['stdout'][:500]}")
            
            all_results.extend(iteration_results)
        
        # 统计结果
        self.print_summary(all_results)
        
        return all_results
    
    def print_summary(self, results):
        """打印测试结果摘要"""
        print("\n\n" + "="*70)
        print("测试结果摘要")
        print("="*70)
        
        # 按测试分组统计
        test_stats = {}
        for result in results:
            test_name = result['test']
            if test_name not in test_stats:
                test_stats[test_name] = {'passed': 0, 'failed': 0, 'total': 0}
            
            test_stats[test_name]['total'] += 1
            if result['success']:
                test_stats[test_name]['passed'] += 1
            else:
                test_stats[test_name]['failed'] += 1
        
        # 打印每个测试的统计
        for test_name, stats in test_stats.items():
            success_rate = (stats['passed'] / stats['total']) * 100
            print(f"\n{test_name}:")
            print(f"  通过: {stats['passed']}/{stats['total']} ({success_rate:.1f}%)")
            if stats['failed'] > 0:
                print(f"  失败: {stats['failed']}")
        
        # 总体统计
        total_tests = len(results)
        total_passed = sum(1 for r in results if r['success'])
        total_failed = total_tests - total_passed
        overall_success_rate = (total_passed / total_tests) * 100 if total_tests > 0 else 0
        
        print(f"\n总体结果:")
        print(f"  总测试数: {total_tests}")
        print(f"  通过: {total_passed} ({overall_success_rate:.1f}%)")
        print(f"  失败: {total_failed}")
        
        if total_failed > 0:
            print("\n⚠️  注意: 这些测试失败是由于C++全局状态管理问题导致的内存错误")
            print("需要在C++层面修复才能完全解决")

def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(description='隔离运行不稳定的测试')
    parser.add_argument('-n', '--iterations', type=int, default=1,
                        help='运行迭代次数 (默认: 1)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='显示详细输出')
    
    args = parser.parse_args()
    
    runner = IsolatedTestRunner()
    results = runner.run_all_tests(num_iterations=args.iterations)
    
    # 如果有失败，返回非零退出码
    failed_count = sum(1 for r in results if not r['success'])
    sys.exit(1 if failed_count > 0 else 0)

if __name__ == "__main__":
    main()