#!/usr/bin/env python3
"""
Final cleanup of broken imports in test files.
"""

import re
import os
from pathlib import Path
from . import pgp_wrapper

def fix_file(file_path):
    """Fix a single file."""
    with open(file_path, 'r') as f:
        content = f.read()
    
    original_content = content
    
    # Fix broken import patterns
    content = re.sub(r'from pgp_wrapper import\s*from pgp_wrapper import', 'from pgp_wrapper import', content)
    content = re.sub(r'from pgp_wrapper importfrom pgp_wrapper import', 'from pgp_wrapper import', content)
    
    # Remove duplicate imports
    lines = content.split('\n')
    seen_imports = set()
    cleaned_lines = []
    
    for line in lines:
        if line.startswith('from pgp_wrapper import') or line.startswith('from . import pgp_wrapper'):
            if line not in seen_imports:
                cleaned_lines.append(line)
                seen_imports.add(line)
        else:
            cleaned_lines.append(line)
    
    content = '\n'.join(cleaned_lines)
    
    # Write back if changed
    if content != original_content:
        with open(file_path, 'w') as f:
            f.write(content)
        return True
    return False

def main():
    """Fix all PGP test files."""
    test_dir = Path("/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation/energyPGP")
    
    fixed_count = 0
    for py_file in test_dir.glob("*.py"):
        if py_file.name not in ['__init__.py', 'fix_pgp_to_context.py', 'fix_remaining_imports.py', 'clean_imports.py', 'final_clean.py']:
            if fix_file(py_file):
                print(f"Fixed: {py_file.name}")
                fixed_count += 1
    
    print(f"\nFixed {fixed_count} files")

if __name__ == "__main__":
    main()
