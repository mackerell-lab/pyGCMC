#!/usr/bin/env python3
"""Add missing imports to all PGP test files"""

import os
import re
from pathlib import Path

# Functions that should be imported from pgp_wrapper instead of being used directly
WRAPPER_FUNCTIONS = [
    'setPMEParameters',
    'initializePMEParameters', 
    'setPGPParameters',
    'precomputeGridPotential',
    'computeSystemEnergyPGP',
    'computeMovementEnergyPGP',
    'calculateMoleculeEnergy',
    'computeSystemEnergyPME',
    'computeMovementEnergyPME',
    'computeSystemEnergyPMEComplete',
    'computeMovementEnergyPMEFixed',
    'computeSystemEnergyPMEFixed',
    'computeSystemEnergyPGPComplete',
    'computeMovementEnergyPGPComplete',
    'resetPGPState',
    'computeSystemEnergyEwald',
    'computeSystemVdwEnergyCutoff',
    'setPGPParametersIndependent',
    'initializePGPParametersIndependent',
    'computeSystemEnergyPGPIndependent',
    'computeMovementEnergyPGPIndependent'
]

def find_undefined_functions(content):
    """Find functions that are used but not defined or imported"""
    # Find all function calls
    function_calls = re.findall(r'\b([a-zA-Z_]\w*)\s*\(', content)
    function_calls = set(function_calls)
    
    # Find all defined functions
    defined_funcs = re.findall(r'def\s+([a-zA-Z_]\w*)\s*\(', content)
    defined_funcs = set(defined_funcs)
    
    # Find all imported functions
    imported_funcs = set()
    # Direct imports: from module import func
    direct_imports = re.findall(r'from\s+[\w.]+\s+import\s+([^\n]+)', content)
    for imp in direct_imports:
        # Handle multiple imports
        funcs = [f.strip() for f in imp.split(',')]
        for f in funcs:
            # Remove parentheses if present
            f = f.strip('()')
            if ' as ' in f:
                # Handle aliased imports
                f = f.split(' as ')[1].strip()
            imported_funcs.add(f)
    
    # Check which wrapper functions are used but not imported
    missing = []
    for func in WRAPPER_FUNCTIONS:
        if func in function_calls and func not in defined_funcs and func not in imported_funcs:
            missing.append(func)
    
    return missing

def add_missing_imports(file_path):
    """Add missing imports to a file"""
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Check if file already imports from pgp_wrapper
    if 'from . import pgp_wrapper' not in content:
        return False
    
    # Find missing functions
    missing = find_undefined_functions(content)
    if not missing:
        return False
    
    # Find the line after "from . import pgp_wrapper"
    lines = content.split('\n')
    insert_idx = None
    for i, line in enumerate(lines):
        if line.strip() == 'from . import pgp_wrapper':
            insert_idx = i + 1
            break
    
    if insert_idx is None:
        return False
    
    # Check if there's already an import from pgp_wrapper
    existing_import_idx = None
    for i in range(insert_idx, len(lines)):
        if lines[i].strip().startswith('from .pgp_wrapper import'):
            existing_import_idx = i
            break
        # Stop if we hit a non-import line
        if lines[i].strip() and not lines[i].strip().startswith(('from', 'import')):
            break
    
    # Format the import
    import_lines = []
    for i in range(0, len(missing), 3):
        batch = missing[i:i+3]
        import_lines.append(f"from .pgp_wrapper import {', '.join(batch)}")
    
    if existing_import_idx is not None:
        # Merge with existing import
        existing_line = lines[existing_import_idx]
        existing_funcs = re.search(r'from .pgp_wrapper import (.+)', existing_line)
        if existing_funcs:
            current_funcs = [f.strip() for f in existing_funcs.group(1).split(',')]
            all_funcs = list(set(current_funcs + missing))
            # Format nicely
            new_import_lines = []
            for i in range(0, len(all_funcs), 3):
                batch = all_funcs[i:i+3]
                new_import_lines.append(f"from .pgp_wrapper import {', '.join(batch)}")
            
            # Replace the existing import
            lines[existing_import_idx] = new_import_lines[0]
            for j, new_line in enumerate(new_import_lines[1:], 1):
                lines.insert(existing_import_idx + j, new_line)
    else:
        # Add new import after "from . import pgp_wrapper"
        for i, import_line in enumerate(import_lines):
            lines.insert(insert_idx + i, import_line)
    
    # Write back
    with open(file_path, 'w') as f:
        f.write('\n'.join(lines))
    
    return True

def main():
    test_dir = Path("/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation/energyPGP")
    
    # Get all Python test files
    py_files = list(test_dir.glob("*.py"))
    
    # Exclude certain files
    exclude_files = ['__init__.py', 'pgp_wrapper.py', 'add_missing_imports.py', 
                     'conftest.py', 'helpers.py', 'run_isolated.py']
    
    fixed_count = 0
    for py_file in py_files:
        if py_file.name in exclude_files or py_file.name.endswith('_helpers.py'):
            continue
            
        print(f"Checking: {py_file.name}")
        if add_missing_imports(py_file):
            print(f"  ✓ Added missing imports")
            fixed_count += 1
        else:
            print(f"  - No missing imports")
    
    print(f"\nFixed {fixed_count} files")

if __name__ == "__main__":
    main()