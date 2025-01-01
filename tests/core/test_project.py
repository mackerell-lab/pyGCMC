#!/usr/bin/env python3

import os
import sys
from pygcmc import Project

def get_test_data_dir():
    """Get the path to the test data directory."""
    print("Getting test data directory...")
    current_dir = os.path.dirname(os.path.abspath(__file__))
    test_data_dir = os.path.join(current_dir, "..", "data")
    print(f"Test data directory: {test_data_dir}")
    return test_data_dir

def check_test_files(test_data_dir):
    """Check if all required test files exist."""
    print("\nChecking test files:")
    files = ["test.pdb", "test.top", "silcs.str"]
    for file in files:
        file_path = os.path.join(test_data_dir, file)
        print(f"Checking {file_path}...")
        if os.path.exists(file_path):
            print(f"Found: {file_path}")
        else:
            print(f"Error: {file_path} not found!")
            sys.exit(1)
    return {name: os.path.join(test_data_dir, name) for name in files}

def main():
    print("\n=== Starting main ===\n")
    test_project()

def test_project():
    print("=== Starting test_project ===")
    
    # Get test data directory
    print("Getting test data directory...")
    test_data_dir = get_test_data_dir()
    print(f"Test data directory: {test_data_dir}\n")

    # Define test files
    test_files = {
        'test.pdb': os.path.join(test_data_dir, 'test.pdb'),
        'test.top': os.path.join(test_data_dir, 'test.top'),
        'silcs.str': os.path.join(test_data_dir, 'silcs.str'),
        'par_all36m_prot.prm': os.path.join(test_data_dir, 'par_all36m_prot.prm'),
        'par_all36_cgenff.prm': os.path.join(test_data_dir, 'par_all36_cgenff.prm'),
        'toppar_water_ions.str': os.path.join(test_data_dir, 'toppar_water_ions.str')
    }

    # Check if test files exist
    print("Checking test files:")
    for name, path in test_files.items():
        print(f"Checking {path}...")
        if os.path.exists(path):
            print(f"Found: {path}")
        else:
            print(f"Error: {path} not found!")
            return
    print()

    # Create project
    print("=== Creating project ===")
    print("About to create Project object...")
    project = Project("test_project")
    print("Project object created successfully")
    print(f"Project name: {project.get_name()}\n")

    # Load structure
    print("=== Loading structure ===")
    print("About to load structure from:")
    print(f"PDB: {test_files['test.pdb']}")
    print(f"TOP: {test_files['test.top']}")
    print("Calling load_structure...")
    structure = project.load_structure(test_files['test.pdb'], test_files['test.top'])
    print("Structure loaded successfully")
    print(f"Structure contains {len(structure.atoms)} atoms\n")

    # Load force field
    print("=== Loading force field ===")
    param_files = [
        test_files['par_all36m_prot.prm'],
        test_files['par_all36_cgenff.prm'],
        test_files['toppar_water_ions.str'],
        test_files['silcs.str']
    ]
    print("About to load force field parameters:")
    for param_file in param_files:
        print(f"  {param_file}")
    print("Calling load_forcefield...")
    forcefield = project.load_forcefield(param_files)
    print("Force field loaded successfully")
    
    # Print force field parameters before applying
    print("\n=== Force field parameters before applying ===")
    forcefield.print_nonbonded_params()
    
    # Apply force field to structure
    print("\n=== Applying force field to structure ===")
    print("Calling apply_forcefield...")
    structure.apply_forcefield(forcefield)
    print("Force field applied successfully\n")

    # Print detailed information
    print("=== Printing detailed information ===")
    project.print_all_atoms()
    project.print_forcefield_info()

if __name__ == '__main__':
    main() 