#!/usr/bin/env python3
"""Fix initializePMEParameters calls to use correct number of arguments"""

import os
import re
from pathlib import Path

def fix_pme_init_calls(file_path):
    """Fix initializePMEParameters calls in a file"""
    with open(file_path, 'r') as f:
        content = f.read()
    
    original_content = content
    
    # Pattern to match initializePMEParameters calls with 5 arguments
    # This matches calls like: initializePMEParameters(cutoff, box, alpha, mesh_size, spline_order)
    pattern = r'initializePMEParameters\s*\(\s*([^,]+),\s*([^,]+),\s*([^,]+),\s*[^,]+,\s*[^)]+\)'
    
    # Replace with 3-argument version
    def replace_call(match):
        cutoff = match.group(1).strip()
        box = match.group(2).strip()
        alpha = match.group(3).strip()
        return f'initializePMEParameters({cutoff}, {box}, {alpha})'
    
    content = re.sub(pattern, replace_call, content)
    
    # Also fix calls with 4 arguments
    pattern2 = r'initializePMEParameters\s*\(\s*([^,]+),\s*([^,]+),\s*([^,]+),\s*[^)]+\)'
    content = re.sub(pattern2, replace_call, content)
    
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
        if py_file.name in ['fix_pme_init_calls.py']:
            continue
            
        print(f"Checking: {py_file.name}")
        if fix_pme_init_calls(py_file):
            print(f"  ✓ Fixed initializePMEParameters calls")
            fixed_count += 1
        else:
            print(f"  - No changes needed")
    
    print(f"\nFixed {fixed_count} files")

if __name__ == "__main__":
    main()