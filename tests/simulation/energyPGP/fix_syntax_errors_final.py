#!/usr/bin/env python3
"""Fix syntax errors caused by corrupted lines"""

import os
import re
from pathlib import Path

def fix_syntax_errors(file_path):
    """Fix syntax errors in a file"""
    with open(file_path, 'r') as f:
        content = f.read()
    
    original_content = content
    
    # Pattern to match corrupted lines like:
    precomputeGridPotential(state)
    pattern = r'precomputeGridPotential\(([^)]+)\):\.(\d+)f\} kJ/mol[^"]*"\)'
    
    def replace_corrupted(match):
        arg = match.group(1)
        return f'precomputeGridPotential({arg})'
    
    content = re.sub(pattern, replace_corrupted, content)
    
    # Also fix lines that have print statements merged with function calls
    pattern2 = r'(precomputeGridPotential\([^)]+\))([^\n]+print\()'
    content = re.sub(pattern2, r'\1\n    \2', content)
    
    if content != original_content:
        with open(file_path, 'w') as f:
            f.write(content)
        return True
    return False

def main():
    test_dir = Path("/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation/energyPGP")
    
    # Get all Python files
    py_files = list(test_dir.glob("*.py"))
    
    fixed_count = 0
    for py_file in py_files:
        if py_file.name in ['fix_syntax_errors_final.py']:
            continue
            
        print(f"Checking: {py_file.name}")
        if fix_syntax_errors(py_file):
            print(f"  ✓ Fixed syntax errors")
            fixed_count += 1
        else:
            print(f"  - No syntax errors found")
    
    print(f"\nFixed {fixed_count} files")

if __name__ == "__main__":
    main()