#!/usr/bin/env python3
"""修复所有 PGP 测试文件的导入问题"""

import os
import re
from pathlib import Path

def fix_imports_in_file(file_path):
    """修复单个文件中的导入"""
    with open(file_path, 'r') as f:
        content = f.read()
    
    original_content = content
    
    # 需要从 pgp_wrapper 导入的函数列表
    pgp_functions = [
        'setPGPParameters', 'initializePMEParameters', 'setPMEParameters',
        'precomputeGridPotential', 'computeSystemEnergyPGP', 'computeMovementEnergyPGP',
        'calculateMoleculeEnergy', 'computeSystemEnergyPME', 'computeMovementEnergyPME',
        'computeSystemEnergyPMEComplete', 'computeMovementEnergyPMEFixed',
        'computeSystemEnergyPMEFixed', 'computeSystemEnergyEwald',
        'computeSystemVdwEnergyCutoff', 'interpolateMoleculeEnergy',
        'resetPGPState', 'computeSystemEnergyPGPComplete', 'computeMovementEnergyPGPComplete'
    ]
    
    # 查找所有从 pygcmc 导入的行
    import_lines = re.findall(r'^from pygcmc import (.+)$', content, re.MULTILINE)
    
    # 收集需要从 pygcmc 和 pgp_wrapper 导入的内容
    pygcmc_imports = []
    wrapper_imports = []
    
    for line in import_lines:
        # 分割导入的项目
        items = [item.strip() for item in line.split(',')]
        for item in items:
            if item in pgp_functions:
                wrapper_imports.append(item)
            else:
                pygcmc_imports.append(item)
    
    # 删除所有现有的 from pygcmc import 行
    content = re.sub(r'^from pygcmc import .+$\n?', '', content, flags=re.MULTILINE)
    
    # 删除现有的 pgp_wrapper 导入（避免重复）
    content = re.sub(r'^from \.pgp_wrapper import .+$\n?', '', content, flags=re.MULTILINE)
    content = re.sub(r'^from \.?pgp_wrapper import .+$\n?', '', content, flags=re.MULTILINE)
    
    # 构建新的导入部分
    new_imports = []
    
    # 保留 import pygcmc
    if 'import pygcmc' in content:
        pass  # 已经有了
    elif pygcmc_imports or wrapper_imports:
        new_imports.append('import pygcmc')
    
    # 添加 pgp_wrapper 导入
    if wrapper_imports or 'from . import pgp_wrapper' in original_content:
        new_imports.append('from . import pgp_wrapper')
    
    # 添加 pygcmc 导入
    if pygcmc_imports:
        # 将导入分组，每行最多3个
        for i in range(0, len(pygcmc_imports), 3):
            group = pygcmc_imports[i:i+3]
            new_imports.append(f"from pygcmc import {', '.join(group)}")
    
    # 添加 pgp_wrapper 导入
    if wrapper_imports:
        # 将导入分组，每行最多3个
        for i in range(0, len(wrapper_imports), 3):
            group = wrapper_imports[i:i+3]
            new_imports.append(f"from .pgp_wrapper import {', '.join(group)}")
    
    # 找到插入点（在其他导入之后）
    # 查找最后一个 import 语句的位置
    import_pattern = re.compile(r'^(import |from .+ import )', re.MULTILINE)
    matches = list(import_pattern.finditer(content))
    
    if matches:
        # 找到最后一个导入语句的结束位置
        last_import_end = max(m.end() for m in matches)
        # 找到该行的结束
        next_newline = content.find('\n', last_import_end)
        if next_newline == -1:
            next_newline = last_import_end
        insert_pos = next_newline + 1
    else:
        # 如果没有导入，在文件开头插入
        insert_pos = 0
    
    # 插入新的导入
    if new_imports:
        import_block = '\n'.join(new_imports) + '\n'
        content = content[:insert_pos] + import_block + content[insert_pos:]
    
    # 清理多余的空行
    content = re.sub(r'\n{3,}', '\n\n', content)
    
    if content != original_content:
        with open(file_path, 'w') as f:
            f.write(content)
        return True
    return False

def main():
    test_dir = Path("/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation/energyPGP")
    
    # 获取所有 Python 文件
    py_files = list(test_dir.glob("*.py"))
    
    # 排除特定文件
    exclude_files = ['__init__.py', 'pgp_wrapper.py', 'fix_all_pgp_imports.py', 
                     'fix_pgp_to_context.py', 'fix_test_energy_pgp_imports.py',
                     'clean_all_imports.py', 'fix_all_imports.py', 'fix_remaining_imports.py']
    
    fixed_count = 0
    for py_file in py_files:
        if py_file.name in exclude_files:
            continue
            
        print(f"检查: {py_file.name}")
        if fix_imports_in_file(py_file):
            print(f"  ✓ 已修复")
            fixed_count += 1
        else:
            print(f"  - 无需修改")
    
    print(f"\n总共修复了 {fixed_count} 个文件")

if __name__ == "__main__":
    main()
