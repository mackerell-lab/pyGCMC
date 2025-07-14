#!/usr/bin/env python3
"""Fix common syntax errors in test files."""

import re
from pathlib import Path

def fix_unterminated_docstring(content):
    """Fix unterminated triple-quoted strings."""
    # Count triple quotes
    single_count = content.count("'''")
    double_count = content.count('"""')
    
    # If odd number, add one at the end
    if single_count % 2 == 1:
        content = content.rstrip() + "\n'''"
    if double_count % 2 == 1:
        content = content.rstrip() + '\n"""'
    
    return content

def fix_file(file_path):
    """Fix syntax errors in a file."""
    try:
        with open(file_path, 'r') as f:
            content = f.read()
    except:
        return False
    
    original = content
    
    # Fix unterminated docstrings
    content = fix_unterminated_docstring(content)
    
    # Ensure file ends with newline
    if not content.endswith('\n'):
        content += '\n'
    
    if content != original:
        with open(file_path, 'w') as f:
            f.write(content)
        return True
    return False

def main():
    test_dir = Path("/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation/energyPGP")
    
    fixed = 0
    for py_file in sorted(test_dir.glob("*.py")):
        if py_file.name.startswith(('fix_', 'clean_', 'check_')):
            continue
        if fix_file(py_file):
            print(f"Fixed: {py_file.name}")
            fixed += 1
    
    print(f"\nFixed {fixed} files")

if __name__ == "__main__":
    main()
