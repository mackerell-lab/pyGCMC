#!/usr/bin/env python3
"""Fix all syntax errors in PGP test files"""

import os
import re
from pathlib import Path

def fix_all_syntax_errors(file_path):
    """Fix all syntax errors in a file"""
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    fixed = False
    new_lines = []
    i = 0
    
    while i < len(lines):
        line = lines[i]
        
        # Fix lines that have just "elec_initial)" or similar
        if re.match(r'^\s*\w+_(initial|final|pme|pgp|complete)\s*\)\s*$', line):
            # This is likely a truncated function call
            # Replace with a proper function call
            var_match = re.match(r'^\s*(\w+)_(initial|final|pme|pgp|complete)\s*\)', line)
            if var_match:
                indent = len(line) - len(line.lstrip())
                prefix = var_match.group(1)
                suffix = var_match.group(2)
                if 'elec' in prefix:
                    new_lines.append(' ' * indent + f'elec_{suffix}, vdw_{suffix}, total_{suffix} = computeSystemEnergyPMEComplete(state)\n')
                else:
                    new_lines.append(line)
                fixed = True
                i += 1
                continue
        
        # Fix corrupted precomputeGridPotential lines
        if 'precomputeGridPotential(state):' in line and '.6f}' in line:
            # Extract the argument
            indent = len(line) - len(line.lstrip())
            new_lines.append(' ' * indent + 'precomputeGridPotential(state)\n')
            fixed = True
            i += 1
            continue
        
        # Fix lines with "No newline at end of file" in the middle
        if 'No newline at end of file' in line and not line.strip().startswith('#'):
            # Remove this text
            line = line.replace('No newline at end of file', '')
            fixed = True
        
        new_lines.append(line)
        i += 1
    
    if fixed:
        with open(file_path, 'w') as f:
            f.writelines(new_lines)
    
    return fixed

def main():
    test_dir = Path("/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation/energyPGP")
    
    # Get all Python files
    py_files = list(test_dir.glob("*.py"))
    
    fixed_count = 0
    for py_file in py_files:
        if py_file.name in ['fix_all_syntax_errors.py']:
            continue
            
        print(f"Checking: {py_file.name}")
        if fix_all_syntax_errors(py_file):
            print(f"  ✓ Fixed syntax errors")
            fixed_count += 1
        else:
            print(f"  - No syntax errors found")
    
    print(f"\nFixed {fixed_count} files")

if __name__ == "__main__":
    main()