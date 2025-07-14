#!/usr/bin/env python3
"""Fix all remaining import issues in PGP test files."""

import os
import re
from pathlib import Path

def fix_imports_in_file(file_path):
    """Fix imports in a single file."""
    with open(file_path, 'r') as f:
        content = f.read()
    
    original_content = content
    
    # Fix non-relative imports of pgp_wrapper
    patterns = [
        (r'^from pgp_wrapper import (.+)$', r'from .pgp_wrapper import \1')
        (r'^import pgp_wrapper$', r'from . import pgp_wrapper')
    ]
    
    for pattern, replacement in patterns:
        content = re.sub(pattern, replacement, content, flags=re.MULTILINE)
    
    if content != original_content:
        with open(file_path, 'w') as f:
            f.write(content)
        return True
    return False

def main():
    test_dir = Path("/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation/energyPGP")
    fixed_count = 0
    
    for py_file in test_dir.glob("*.py"):
        if py_file.name in ['__init__.py', 'pgp_wrapper.py', 'fix_all_imports.py']:
            continue
            
        if fix_imports_in_file(py_file):
            print(f"Fixed imports in: {py_file.name}")
            fixed_count += 1
    
    print(f"\nFixed {fixed_count} files")

if __name__ == "__main__":
    main()
