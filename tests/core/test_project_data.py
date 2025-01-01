#!/usr/bin/env python3

import os
import sys
from pygcmc import Project

def get_test_data_dir():
    """Get the path to the test data directory."""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    print("Getting test data directory...")
    test_dir = os.path.join(current_dir, "..", "data")
    print(f"Test data directory: {test_dir}\n")
    return test_dir

def test_structure_data_access():
    """Test the new data access methods for structure information."""
    print("\n=== Testing Structure Data Access ===")
    
    # Setup
    test_data_dir = get_test_data_dir()
    print("Creating project...")
    project = Project("test_project")
    print("Project created successfully")
    
    # Load structure and force field
    pdb_file = os.path.join(test_data_dir, "test.pdb")
    top_file = os.path.join(test_data_dir, "test.top")
    print(f"\nLoading structure from:")
    print(f"PDB: {pdb_file}")
    print(f"TOP: {top_file}")
    structure = project.load_structure(pdb_file, top_file)
    print("Structure loaded successfully")
    
    # Load force field
    print("\nLoading force field parameters:")
    param_files = [
        os.path.join(test_data_dir, "par_all36m_prot.prm"),
        os.path.join(test_data_dir, "par_all36_cgenff.prm"),
        os.path.join(test_data_dir, "toppar_water_ions.str"),
        os.path.join(test_data_dir, "silcs.str")
    ]
    for param_file in param_files:
        print(f"  {param_file}")
    ff = project.load_forcefield(param_files)
    print("Force field loaded successfully")
    
    # Apply force field to structure
    print("\nApplying force field to structure...")
    structure.apply_forcefield(ff)
    print("Force field applied successfully")
    
    # Get nonbonded parameters for quick lookup
    nb_params = {param['atom_type']: param for param in ff.get_nonbonded_parameters()}
    
    # Print all atom information
    print("\nDetailed Atom Information:")
    print("-" * 120)
    print(f"{'Residue':10s} {'Seq':6s} {'Name':6s} {'Type':8s} {'TopoType':10s} {'Charge':10s} {'Mass':10s} {'Epsilon':10s} {'Rmin':10s}")
    print("-" * 120)
    atoms_data = structure.get_atoms_data()
    for atom in atoms_data:
        ff_params = nb_params.get(atom['topo_type'], {'epsilon': 0.0, 'rmin': 0.0})
        print(f"{atom['residue']:10s} {atom['sequence']:6d} {atom['name']:6s} {atom['type']:8s} "
              f"{atom['topo_type']:10s} {atom['topo_charge']:10.3f} {atom['topo_mass']:10.3f} "
              f"{ff_params['epsilon']:10.3f} {ff_params['rmin']:10.3f}")
    print("-" * 120)
    print()

def test_forcefield_data_access():
    """Test the new data access methods for force field information."""
    print("\n=== Testing Force Field Data Access ===")
    
    # Setup
    test_data_dir = get_test_data_dir()
    print("Creating project...")
    project = Project("test_project")
    print("Project created successfully")
    
    # Load structure and force field
    pdb_file = os.path.join(test_data_dir, "test.pdb")
    top_file = os.path.join(test_data_dir, "test.top")
    print(f"\nLoading structure from:")
    print(f"PDB: {pdb_file}")
    print(f"TOP: {top_file}")
    structure = project.load_structure(pdb_file, top_file)
    print("Structure loaded successfully")
    
    # Load force field
    print("\nLoading force field parameters:")
    param_files = [
        os.path.join(test_data_dir, "par_all36m_prot.prm"),
        os.path.join(test_data_dir, "par_all36_cgenff.prm"),
        os.path.join(test_data_dir, "toppar_water_ions.str"),
        os.path.join(test_data_dir, "silcs.str")
    ]
    for param_file in param_files:
        print(f"  {param_file}")
    ff = project.load_forcefield(param_files)
    print("Force field loaded successfully")
    
    # Apply force field to structure
    print("\nApplying force field to structure...")
    structure.apply_forcefield(ff)
    print("Force field applied successfully")
    
    # Print force field information
    print("\nForce Field Information:")
    print("------------------------")
    print("\nNonbonded Parameters:")
    print("-" * 80)
    print(f"{'Type':10s} {'Epsilon':15s} {'Rmin':15s}")
    print("-" * 80)
    nb_params = ff.get_nonbonded_parameters()
    for param in sorted(nb_params, key=lambda x: x['atom_type']):
        print(f"{param['atom_type']:10s} {param['epsilon']:15.3f} {param['rmin']:15.3f}")
    print("-" * 80)
    
    # Print NBFIX parameters
    print("\nNBFIX Parameters:")
    nbfix_params = ff.get_nbfix_parameters()
    print(f"Total NBFIX parameters: {len(nbfix_params)}")
    if nbfix_params:
        print("\nAll NBFIX parameters:")
        print("-" * 80)
        print(f"{'Type1':10s} {'Type2':10s} {'Epsilon':15s} {'Rmin':15s}")
        print("-" * 80)
        for param in sorted(nbfix_params, key=lambda x: (x['type1'], x['type2'])):
            print(f"{param['type1']:10s} {param['type2']:10s} {param['epsilon']:15.3f} {param['rmin']:15.3f}")
        print("-" * 80)
    
    # Print global parameters
    print("\nGlobal Parameters:")
    global_params = ff.get_global_parameters()
    for param, value in sorted(global_params.items()):
        print(f"  {param}: {value:.3f} Å")
    print()

def test_energy_data_access():
    """Test the new data access methods for energy calculations."""
    print("\n=== Testing Energy Data Access ===")
    
    # Setup
    test_data_dir = get_test_data_dir()
    print("Creating project...")
    project = Project("test_project")
    print("Project created successfully")
    
    # Load structure and force field
    print("\nLoading structure and force field...")
    pdb_file = os.path.join(test_data_dir, "test.pdb")
    top_file = os.path.join(test_data_dir, "test.top")
    structure = project.load_structure(pdb_file, top_file)
    
    param_files = [
        os.path.join(test_data_dir, "par_all36m_prot.prm"),
        os.path.join(test_data_dir, "par_all36_cgenff.prm"),
        os.path.join(test_data_dir, "toppar_water_ions.str"),
        os.path.join(test_data_dir, "silcs.str")
    ]
    ff = project.load_forcefield(param_files)
    print("Structure and force field loaded successfully")
    
    # Apply force field to structure
    print("\nApplying force field to structure...")
    structure.apply_forcefield(ff)
    print("Force field applied successfully")
    
    print()

if __name__ == '__main__':
    print("=== Starting main ===\n")
    test_structure_data_access()
    test_forcefield_data_access()
    test_energy_data_access()
    print("\nAll tests passed successfully!") 