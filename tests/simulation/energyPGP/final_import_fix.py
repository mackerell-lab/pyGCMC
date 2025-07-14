#!/usr/bin/env python3
"""Final comprehensive fix for all malformed imports"""

import os
import re
from pathlib import Path

def fix_malformed_imports_comprehensive(file_path):
    """Fix all types of malformed import statements"""
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    fixed = False
    new_lines = []
    i = 0
    
    while i < len(lines):
        line = lines[i]
        
        # Check if this is the start of a malformed import block
        # Pattern: from .module import (
        if (line.strip().startswith('from .') and 'import (' in line and 
            i + 1 < len(lines) and lines[i+1].strip().startswith('from pygcmc import')):
            
            # Found malformed import
            fixed = True
            
            # Extract the helper import start
            helper_import_start = line
            
            # Extract the pygcmc import
            pygcmc_import = lines[i+1]
            new_lines.append(pygcmc_import)
            
            # Now add the helper import start
            new_lines.append(helper_import_start)
            
            # Skip past the pygcmc import line
            i += 2
            
            # Continue adding the rest of the helper imports until we find the closing )
            while i < len(lines):
                if ')' in lines[i] and not lines[i].strip().startswith('from'):
                    new_lines.append(lines[i])
                    i += 1
                    break
                new_lines.append(lines[i])
                i += 1
        else:
            new_lines.append(line)
            i += 1
    
    if fixed:
        # Write back the fixed content
        with open(file_path, 'w') as f:
            f.writelines(new_lines)
    
    return fixed

def main():
    test_dir = Path("/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation/energyPGP")
    
    # Get all Python files
    py_files = list(test_dir.glob("*.py"))
    
    fixed_count = 0
    for py_file in py_files:
        if py_file.name in ['final_import_fix.py']:
            continue
            
        print(f"Checking: {py_file.name}")
        if fix_malformed_imports_comprehensive(py_file):
            print(f"  ✓ Fixed malformed imports")
            fixed_count += 1
        else:
            print(f"  - No malformed imports found")
    
    print(f"\nFixed {fixed_count} files")

if __name__ == "__main__":
    main()