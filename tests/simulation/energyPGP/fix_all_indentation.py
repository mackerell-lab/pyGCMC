#!/usr/bin/env python3
"""Fix all indentation issues in PGP test files"""

import os
import re
from pathlib import Path

def fix_indentation_issues(file_path):
    """Fix common indentation patterns"""
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    fixed_lines = []
    i = 0
    while i < len(lines):
        line = lines[i]
        
        # Pattern: line that starts with "    setPGPParameters" (4 spaces) should be 8 spaces
        if line.startswith('    setPGPParameters(') and not line.startswith('        '):
            fixed_lines.append('        ' + line.strip() + '\n')
        # Pattern: line after setPGPParameters that has wrong indentation
        elif i > 0 and 'setPGPParameters(' in lines[i-1] and line.strip() and not line.startswith('#'):
            # Check if it's a continuation of function arguments
            if not line.strip().startswith(')'):
                # It's probably the next statement, ensure proper indentation
                stripped = line.strip()
                if stripped.startswith('precomputeGridPotential(') or stripped.startswith('computeSystemEnergyPGP('):
                    fixed_lines.append('        ' + stripped + '\n')
                else:
                    fixed_lines.append(line)
            else:
                fixed_lines.append(line)
        else:
            fixed_lines.append(line)
        i += 1
    
    with open(file_path, 'w') as f:
        f.writelines(fixed_lines)

def main():
    test_dir = Path("/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation/energyPGP")
    
    # Files that need fixing
    files_to_fix = [
        "pgp_grid_convergence.py",
        "pgp_cutoff_continuity.py"
    ]
    
    for file_name in files_to_fix:
        file_path = test_dir / file_name
        if file_path.exists():
            print(f"Fixing {file_name}...")
            fix_indentation_issues(file_path)
    
    print("Done!")

if __name__ == "__main__":
    main()