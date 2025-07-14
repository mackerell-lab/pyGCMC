#!/usr/bin/env python3
"""Check syntax of all Python files."""

import ast
import sys
from pathlib import Path

def check_file(file_path):
    """Check syntax of a single file."""
    try:
        with open(file_path, 'r') as f:
            content = f.read()
        ast.parse(content)
        return True, None
    except SyntaxError as e:
        return False, str(e)

def main():
    test_dir = Path("/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation/energyPGP")
    
    errors = []
    for py_file in sorted(test_dir.glob("*.py")):
        if py_file.name.startswith(('fix_', 'clean_', 'check_')):
            continue
        ok, error = check_file(py_file)
        if not ok:
            errors.append((py_file.name, error))
            print(f"✗ {py_file.name}: {error}")
    
    if errors:
        print(f"\nFound {len(errors)} files with syntax errors")
        sys.exit(1)
    else:
        print("✓ All files have valid syntax")

if __name__ == "__main__":
    main()
