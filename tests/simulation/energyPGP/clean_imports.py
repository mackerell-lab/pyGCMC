#!/usr/bin/env python3
"""
Clean up duplicate and malformed imports in PGP test files.
"""

import re
import os
from pathlib import Path
from . import pgp_wrapper

def clean_imports_in_file(file_path):
    """Clean imports in a single file."""
    with open(file_path, 'r') as f:
        content = f.read()
    
    original_content = content
    
    # Remove duplicate lines
    lines = content.split('\n')
    cleaned_lines = []
    seen_lines = set()
    
    for line in lines:
        # Skip malformed import lines
        if 'from pgp_wrapper importfrom pgp_wrapper import' in line:
            # Extract the actual imports
            match = re.search(r'from pgp_wrapper import\s*([^\n]+)$', line)
            if match:
                line = f'from pgp_wrapper import {match.group(1)}'
        
        # Skip duplicate lines
        if line.strip() and line not in seen_lines:
            cleaned_lines.append(line)
            seen_lines.add(line)
        elif not line.strip():
            # Keep empty lines for formatting
            cleaned_lines.append(line)
    
    content = '\n'.join(cleaned_lines)
    
    # Remove multiple consecutive empty lines
    content = re.sub(r'\n\n\n+', '\n\n', content)
    
    # Remove duplicate "from . import pgp_wrapper" lines
    content = re.sub(r'(from \. import pgp_wrapper\n)+', r'from . import pgp_wrapper\n', content)
    
    # Fix specific import patterns
    content = re.sub(r'from \. import pgp_wrapper\nfrom pgp_wrapper import initializePMEParameters\n.*?from pgp_wrapper import initializePMEParameters\n', 
                    r'from . import pgp_wrapper\nfrom pgp_wrapper import initializePMEParameters\n', content, flags=re.DOTALL)
    
    # Write back if changed
    if content != original_content:
        with open(file_path, 'w') as f:
            f.write(content)
        return True
    return False

def main():
    """Clean all PGP test files."""
    test_dir = Path("/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation/energyPGP")
    
    fixed_count = 0
    for py_file in test_dir.glob("*.py"):
        if py_file.name not in ['__init__.py', 'fix_pgp_to_context.py', 'fix_remaining_imports.py', 'clean_imports.py']:
            if clean_imports_in_file(py_file):
                print(f"Cleaned: {py_file.name}")
                fixed_count += 1
    
    print(f"\nCleaned {fixed_count} files")

if __name__ == "__main__":
    main()
