#!/usr/bin/env python3
"""Fix setPGPParameters calls to use correct number of arguments"""

import os
import re
from pathlib import Path

def fix_pgp_params_calls(file_path):
    """Fix setPGPParameters calls in a file"""
    with open(file_path, 'r') as f:
        content = f.read()
    
    original_content = content
    
    # Pattern to match setPGPParameters calls with only 1 argument (alpha)
    # This matches calls like: setPGPParameters(alpha)
    pattern1 = r'setPGPParameters\s*\(\s*([^,\)]+)\s*\)'
    
    def replace_single_arg(match):
        alpha = match.group(1).strip()
        # Use default values for mesh_size and other parameters
        return f'setPGPParameters({alpha}, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)'
    
    content = re.sub(pattern1, replace_single_arg, content)
    
    # Pattern to match setPGPParameters calls with wrong number of arguments
    # Fix calls that look like: setPGPParameters(alpha, mesh_size, cutoff, mesh_size, spline_order, tolerance)
    # but are missing some arguments
    pattern2 = r'setPGPParameters\s*\(\s*([^,]+),\s*([^,]+),\s*([^,]+),\s*([^,]+),\s*([^,]+)\s*\)'
    
    def replace_five_args(match):
        args = [match.group(i).strip() for i in range(1, 6)]
        # Add missing tolerance parameter
        return f'setPGPParameters({args[0]}, {args[1]}, {args[2]}, {args[3]}, {args[4]}, 1e-6)'
    
    content = re.sub(pattern2, replace_five_args, content)
    
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
        if py_file.name in ['fix_pgp_params_calls.py', 'pgp_wrapper.py']:
            continue
            
        print(f"Checking: {py_file.name}")
        if fix_pgp_params_calls(py_file):
            print(f"  ✓ Fixed setPGPParameters calls")
            fixed_count += 1
        else:
            print(f"  - No changes needed")
    
    print(f"\nFixed {fixed_count} files")

if __name__ == "__main__":
    main()