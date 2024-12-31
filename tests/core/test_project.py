# tests/core/test_project.py

import os
import sys
from math import isnan

print("=== Starting imports ===")
print("Importing pygcmc...")
import pygcmc as gcmc
print("Successfully imported pygcmc")

def test_project():
    print("\n=== Starting test_project ===")
    try:
        # Get the directory containing test data
        print("Getting test data directory...")
        test_data_dir = os.path.join(os.path.dirname(__file__), "..", "data")
        print(f"Test data directory: {test_data_dir}")
        
        # Verify test files exist
        print("\nSetting up file paths...")
        pdb_file = os.path.join(test_data_dir, "test.pdb")
        top_file = os.path.join(test_data_dir, "test.top")
        param_file = os.path.join(test_data_dir, "silcs.str")
        
        print("\nChecking test files:")
        for file_path in [pdb_file, top_file, param_file]:
            print(f"Checking {file_path}...")
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"Required test file not found: {file_path}")
            print(f"Found: {file_path}")
        
        print("\n=== Creating project ===")
        print("About to create Project object...")
        project = gcmc.Project("test_project")
        print("Project object created successfully")
        print(f"Project name: {project.get_name()}")
        
        print("\n=== Loading structure ===")
        print(f"About to load structure from:")
        print(f"PDB: {pdb_file}")
        print(f"TOP: {top_file}")
        print("Calling load_structure...")
        structure = project.load_structure(pdb_file, top_file)
        print("Structure loaded successfully")
        print(f"Structure contains {len(structure)} atoms")
        
        print("\n=== Loading force field ===")
        print(f"About to load force field from: {param_file}")
        print("Creating parameter files list...")
        param_files = [param_file]
        print("Calling load_forcefield...")
        ff = project.load_forcefield(param_files)
        print("Force field loaded successfully")

        # Apply force field to structure
        print("\n=== Applying force field to structure ===")
        structure.apply_forcefield(ff)
        print("Force field applied successfully")
        
        # Print detailed information
        print("\n=== Force Field Parameters ===")
        print(f"Global parameters:")
        print(f"  Cutoff: {ff.cutoff}")
        print(f"  Switching: {ff.switching}")
        print(f"  Pairlist distance: {ff.pairlist_distance}")
        
        print("\n=== Nonbonded Parameters ===")
        nonbonded_params = ff.nonbonded_params
        print(f"Total nonbonded parameters: {len(nonbonded_params)}")
        print("\nSample nonbonded parameters:")
        for atom_type, params in list(nonbonded_params.items())[:5]:
            print(f"  {atom_type}: epsilon={params.epsilon}, rmin={params.rmin}")
        
        print("\n=== NBFIX Parameters ===")
        nbfix_params = ff.nbfix_params
        print(f"Total NBFIX parameters: {len(nbfix_params)}")
        print("\nAll NBFIX parameters:")
        for (type1, type2), params in nbfix_params.items():
            print(f"  {type1}-{type2}: epsilon={params.epsilon}, rmin={params.rmin}")
        
        print("\n=== Structure Information ===")
        try:
            # Get total number of atoms without storing the list
            print(f"Total atoms: {len(structure)}")
            
            print("\nSample atom information (first 5 atoms):")
            atom_count = 0
            for atom in structure.atoms:
                if atom_count >= 5:
                    break
                    
                print(f"\nAtom {atom_count + 1}:")
                try:
                    print(f"  Name: {atom.name}")
                    print(f"  Type: {atom.type}")
                    print(f"  Position: ({atom.x:.3f}, {atom.y:.3f}, {atom.z:.3f})")
                    
                    # Print topology information if available
                    if hasattr(atom, 'topo_type') and atom.topo_type:
                        print(f"  Topology type: {atom.topo_type}")
                        if hasattr(atom, 'topo_charge'):
                            print(f"  Charge: {atom.topo_charge:.3f}")
                        if hasattr(atom, 'topo_mass'):
                            print(f"  Mass: {atom.topo_mass:.3f}")
                    
                    # Print force field parameters if available
                    if (hasattr(atom, 'forcefield_epsilon') and 
                        not isnan(atom.forcefield_epsilon) and 
                        hasattr(atom, 'forcefield_rmin') and 
                        not isnan(atom.forcefield_rmin)):
                        print(f"  Force field parameters:")
                        print(f"    Epsilon: {atom.forcefield_epsilon:.3f}")
                        print(f"    Rmin: {atom.forcefield_rmin:.3f}")
                except Exception as e:
                    print(f"  Error printing atom details: {str(e)}")
                    continue
                
                atom_count += 1
                    
        except Exception as e:
            print(f"Error accessing structure information: {str(e)}")
        
        # Clear references to help with cleanup
        structure = None
        ff = None
        project = None
        
        print("\n=== Test completed successfully ===")
        return True
            
    except Exception as e:
        print(f"\nERROR during test execution: {str(e)}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return False

def main():
    print("\n=== Starting main ===")
    try:
        success = test_project()
        if success:
            print("\n=== Test completed successfully ===")
            return 0
        else:
            print("\n=== Test failed ===")
            return 1
    except Exception as e:
        print(f"\nFATAL ERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    print("=== Starting test script ===")
    result = main()
    print(f"\n=== Script finished with exit code {result} ===")
    sys.exit(result) 