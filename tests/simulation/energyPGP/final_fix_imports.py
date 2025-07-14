#!/usr/bin/env python3
"""最终修复所有 PGP 测试文件的导入问题"""

import os
import re
from pathlib import Path

def fix_imports_in_file(file_path):
    """修复单个文件中的导入"""
    with open(file_path, 'r') as f:
        content = f.read()
    
    original_content = content
    
    # 1. 删除重复的 "from . import pgp_wrapper"
    lines = content.split('\n')
    cleaned_lines = []
    seen_pgp_wrapper_import = False
    
    for line in lines:
        if line.strip() == 'from . import pgp_wrapper':
            if not seen_pgp_wrapper_import:
                cleaned_lines.append(line)
                seen_pgp_wrapper_import = True
        else:
            cleaned_lines.append(line)
    
    content = '\n'.join(cleaned_lines)
    
    # 2. 删除重复的导入
    # 找到所有导入行并去重
    import_pattern = re.compile(r'^from \.pgp_wrapper import (.+)$', re.MULTILINE)
    imports = []
    for match in import_pattern.finditer(content):
        items = [item.strip() for item in match.group(1).split(',')]
        imports.extend(items)
    
    # 去重
    unique_imports = list(dict.fromkeys(imports))
    
    # 删除所有 pgp_wrapper 导入
    content = re.sub(r'^from \.pgp_wrapper import .+$\n?', '', content, flags=re.MULTILINE)
    
    # 在 "from . import pgp_wrapper" 后面添加去重后的导入
    if unique_imports and 'from . import pgp_wrapper' in content:
        # 将导入分组，每行最多4个
        import_lines = []
        for i in range(0, len(unique_imports), 4):
            group = unique_imports[i:i+4]
            import_lines.append(f"from .pgp_wrapper import {', '.join(group)}")
        
        import_block = '\n'.join(import_lines)
        content = content.replace(
            'from . import pgp_wrapper',
            f'from . import pgp_wrapper\n{import_block}'
        )
    
    # 3. 修复遗留的 "No newline at end of file"
    if content and not content.endswith('\n'):
        content += '\n'
    
    # 4. 删除文件末尾的 "No newline at end of file" 文本
    content = re.sub(r'\s*\s*$', '', content)
    
    # 5. 清理多余的空行
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
    exclude_files = ['__init__.py', 'pgp_wrapper.py', 'final_fix_imports.py']
    
    fixed_count = 0
    for py_file in py_files:
        if py_file.name in exclude_files:
            continue
            
        print(f"修复: {py_file.name}")
        if fix_imports_in_file(py_file):
            print(f"  ✓ 已修复")
            fixed_count += 1
        else:
            print(f"  - 无需修改")
    
    print(f"\n总共修复了 {fixed_count} 个文件")

if __name__ == "__main__":
    main()