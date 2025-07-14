#!/usr/bin/env python3
"""Fix indentation issues in PGP test files"""

import os
import re
from pathlib import Path

def fix_indentation(file_path):
    """Fix indentation issues"""
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    fixed_lines = []
    for i, line in enumerate(lines):
        # Fix lines that start with setPGPParameters without proper indentation
        if line.strip().startswith('setPGPParameters(') and not line.startswith('    '):
            fixed_lines.append('        ' + line.strip() + '\n')
        else:
            fixed_lines.append(line)
    
    with open(file_path, 'w') as f:
        f.writelines(fixed_lines)

def main():
    test_dir = Path("/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation/energyPGP")
    
    files_to_fix = [
        "pgp_grid_convergence.py",
        "pgp_cutoff_continuity.py"
    ]
    
    for file_name in files_to_fix:
        file_path = test_dir / file_name
        if file_path.exists():
            print(f"Fixing indentation in {file_name}...")
            fix_indentation(file_path)
    
    print("Done!")

if __name__ == "__main__":
    main()