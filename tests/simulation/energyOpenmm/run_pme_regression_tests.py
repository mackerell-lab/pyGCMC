#!/usr/bin/env python
"""
独立运行 PME 回归测试脚本

由于 PME 使用全局状态，这些测试不能在并行测试环境中运行。
使用此脚本单独验证 PME Total 能量计算的正确性。
"""

import sys
import os

# 添加路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# 确保 pygcmc 可以被导入
if 'PYTHONPATH' not in os.environ:
    os.environ['PYTHONPATH'] = ''
if 'modules/bindings' not in os.environ['PYTHONPATH']:
    os.environ['PYTHONPATH'] = os.path.abspath('modules/bindings') + ':' + os.environ['PYTHONPATH']
    
# 重新设置 sys.path
bindings_path = os.path.abspath("modules/bindings")
if bindings_path not in sys.path:
    sys.path.insert(0, bindings_path)

# 导入测试函数（移除 @pytest.mark.skip 装饰器的影响）
from pme_residue_offset_pattern import test_residue_offset_pattern
from analyze_pme_total_bug import test_pme_total_configurations

def run_tests():
    """运行 PME 回归测试"""
    print("=" * 80)
    print("运行 PME Total 能量计算回归测试")
    print("=" * 80)
    
    # 测试 1：残基偏移模式测试
    print("\n测试 1: test_residue_offset_pattern")
    print("-" * 40)
    try:
        # 获取原始函数（跳过装饰器）
        if hasattr(test_residue_offset_pattern, '__wrapped__'):
            test_func = test_residue_offset_pattern.__wrapped__
        else:
            test_func = test_residue_offset_pattern
        
        test_func()
        print("✅ test_residue_offset_pattern 通过")
    except Exception as e:
        print(f"❌ test_residue_offset_pattern 失败: {e}")
        return False
    
    # 测试 2：PME Total 配置测试
    print("\n测试 2: test_pme_total_configurations")
    print("-" * 40)
    try:
        # 获取原始函数（跳过装饰器）
        if hasattr(test_pme_total_configurations, '__wrapped__'):
            test_func = test_pme_total_configurations.__wrapped__
        else:
            test_func = test_pme_total_configurations
            
        test_func()
        print("✅ test_pme_total_configurations 通过")
    except Exception as e:
        print(f"❌ test_pme_total_configurations 失败: {e}")
        return False
    
    print("\n" + "=" * 80)
    print("✅ 所有 PME 回归测试通过！")
    print("PME Total 能量计算公式已正确修复：total = real_space + reciprocal + self")
    print("=" * 80)
    return True

if __name__ == "__main__":
    # 运行测试
    success = run_tests()
    sys.exit(0 if success else 1)