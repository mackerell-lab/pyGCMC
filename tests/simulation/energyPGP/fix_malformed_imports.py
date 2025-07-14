#!/usr/bin/env python3
"""Fix malformed import statements in PGP test files"""

import os
import re
from pathlib import Path

def fix_malformed_imports(file_path):
    """Fix malformed import statements in a file"""
    with open(file_path, 'r') as f:
        content = f.read()
    
    original_content = content
    
    # Pattern to find malformed imports where 'from pygcmc' appears inside another import
    # This matches patterns like:
    # from .module import (
    # from pygcmc import Something
    #     func1,
    #     func2
    # )
    pattern = re.compile(
        r'(from \.[^\s]+ import \(\s*\n)(from pygcmc import [^\n]+\n)((?:\s+[^\n]+(?:,\s*\n)?)+\))',
        re.MULTILINE
    )
    
    def fix_import(match):
        helper_import = match.group(1)
        pygcmc_import = match.group(2)
        helper_functions = match.group(3)
        # Return the pygcmc import first, then the helper import
        return pygcmc_import + helper_import + helper_functions
    
    content = pattern.sub(fix_import, content)
    
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
        if py_file.name in ['fix_malformed_imports.py']:
            continue
            
        print(f"Checking: {py_file.name}")
        if fix_malformed_imports(py_file):
            print(f"  ✓ Fixed malformed imports")
            fixed_count += 1
        else:
            print(f"  - No malformed imports found")
    
    print(f"\nFixed {fixed_count} files")

if __name__ == "__main__":
    main()