#!/usr/bin/env python3
"""Find and fix remaining syntax errors in PGP test files"""

import os
import re
from pathlib import Path

def find_syntax_errors(file_path):
    """Find common syntax error patterns"""
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    errors_found = []
    
    # Pattern 1: for loops with missing iterator
    for i, line in enumerate(lines):
        if re.match(r'^\s*for\s+\w+\):', line):
            errors_found.append((i, 'incomplete_for_loop', line))
    
    # Pattern 2: assert statements with extra parenthesis
    for i, line in enumerate(lines):
        if re.match(r'^\s*assert.*\)\s*\):', line):
            errors_found.append((i, 'extra_parenthesis', line))
    
    # Pattern 3: Lines ending with colon but no content
    for i, line in enumerate(lines):
        if re.match(r'^.*[^"\']\s*:\s*$', line) and i+1 < len(lines):
            next_line = lines[i+1]
            if not re.match(r'^\s+', next_line) and not next_line.strip().startswith('"""'):
                errors_found.append((i, 'missing_body', line))
    
    return errors_found

def main():
    test_dir = Path("/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation/energyPGP")
    
    # Get all Python files
    py_files = list(test_dir.glob("*.py"))
    
    print("Scanning for syntax errors...")
    
    files_with_errors = []
    for py_file in py_files:
        if py_file.name.startswith('fix_'):
            continue
            
        errors = find_syntax_errors(py_file)
        if errors:
            files_with_errors.append((py_file, errors))
            print(f"\n{py_file.name}:")
            for line_num, error_type, line in errors:
                print(f"  Line {line_num+1}: {error_type}")
                print(f"    {line.strip()}")
    
    if not files_with_errors:
        print("No obvious syntax errors found.")
    else:
        print(f"\nFound potential errors in {len(files_with_errors)} files")

if __name__ == "__main__":
    main()