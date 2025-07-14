#!/usr/bin/env python3
"""Fix setPGPParameters calls with keyword arguments"""

import os
import re
from pathlib import Path

def fix_keyword_args(file_path):
    """Fix setPGPParameters calls with keyword arguments"""
    with open(file_path, 'r') as f:
        content = f.read()
    
    original_content = content
    
    # Pattern to match setPGPParameters calls with alpha= keyword
    pattern = r'setPGPParameters\(alpha=([^,]+),'
    
    def replace_keyword(match):
        alpha_value = match.group(1).strip()
        return f'setPGPParameters({alpha_value},'
    
    content = re.sub(pattern, replace_keyword, content)
    
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
        if py_file.name in ['fix_keyword_args.py']:
            continue
            
        print(f"Checking: {py_file.name}")
        if fix_keyword_args(py_file):
            print(f"  ✓ Fixed setPGPParameters keyword arguments")
            fixed_count += 1
        else:
            print(f"  - No changes needed")
    
    print(f"\nFixed {fixed_count} files")

if __name__ == "__main__":
    main()