#!/usr/bin/env python3
"""Clean all duplicate and malformed imports in PGP test files."""

import os
import re
from pathlib import Path

def clean_imports_in_file(file_path):
    """Clean imports in a single file."""
    with open(file_path, 'r') as f:
        content = f.read()
    
    original_content = content
    
    # Fix duplicate imports on same line
    content = re.sub(r'from \.pgp_wrapper import ([^,\n]+), \1', r'from .pgp_wrapper import \1', content)
    
    # Fix lines with trailing commas after import statements
    content = re.sub(r'(from \.pgp_wrapper import [^,\n]+)\s*,\s*$', r'\1', content, flags=re.MULTILINE)
    
    # Fix duplicate import lines
    lines = content.split('\n')
    seen_imports = set()
    cleaned_lines = []
    
    for line in lines:
        if line.strip().startswith('from .pgp_wrapper import'):
            if line not in seen_imports:
                seen_imports.add(line)
                cleaned_lines.append(line)
        else:
            cleaned_lines.append(line)
    
    content = '\n'.join(cleaned_lines)
    
    # Remove empty import statements
    content = re.sub(r'^from \.pgp_wrapper import\s*$', '', content, flags=re.MULTILINE)
    
    # Fix multiple commas
    content = re.sub(r',\s*,', ',', content)
    
    # Remove trailing commas in imports
    content = re.sub(r'(from [^,\n]+import [^,\n]+),\n', r'\1\n', content)
    
    # Fix syntax errors with duplicate imports on consecutive lines
    content = re.sub(r'(from \.pgp_wrapper import [^\n]+)\n,\s*([a-zA-Z_]+)', r'\1, \2', content)
    
    # Remove standalone comma lines after imports
    content = re.sub(r'(from \.pgp_wrapper import [^\n]+)\n,\s*$', r'\1', content, flags=re.MULTILINE)
    
    if content != original_content:
        with open(file_path, 'w') as f:
            f.write(content)
        return True
    return False

def main():
    test_dir = Path("/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation/energyPGP")
    fixed_count = 0
    
    for py_file in test_dir.glob("*.py"):
        if py_file.name in ['__init__.py', 'pgp_wrapper.py', 'clean_all_imports.py']:
            continue
            
        if clean_imports_in_file(py_file):
            print(f"Cleaned imports in: {py_file.name}")
            fixed_count += 1
    
    print(f"\nCleaned {fixed_count} files")

if __name__ == "__main__":
    main()
