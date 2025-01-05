import os
from pygcmc import Structure

def print_structure_info(pdb_file, psf_files=None, top_file=None):
    """Load structure and print all residue and atom information.
    
    Args:
        pdb_file (str): Path to PDB file
        psf_files (list[str], optional): List of PSF files
        top_file (str, optional): Path to TOP file
    """
    # Create and load structure
    if psf_files:
        # Load with PSF files
        print(f"Loading PDB file: {pdb_file}")
        for psf in psf_files:
            print(f"Loading PSF file: {psf}")
        structure = Structure(pdb=pdb_file, psf=psf_files)
    elif top_file:
        # Load with TOP file
        print(f"Loading PDB file: {pdb_file}")
        print(f"Loading TOP file: {top_file}")
        structure = Structure(pdb=pdb_file, top=top_file)
    else:
        raise ValueError("Must provide either psf_files or top_file")

    # Print basic structure info
    print(f"\nStructure Summary:")
    print(f"Total residues: {len(structure.residues)}")
    print(f"Total atoms: {len(structure.atoms)}")

    # Group residues by name for better organization
    residues_by_name = {}
    for residue in structure.residues:
        if residue.name not in residues_by_name:
            residues_by_name[residue.name] = []
        residues_by_name[residue.name].append(residue)

    # Print residue type summary
    print("\nResidue Type Summary:")
    print("-" * 80)
    for res_name, residues in sorted(residues_by_name.items()):
        print(f"{res_name}: {len(residues)} residues")
        # Print sequence numbers for this residue type
        seq_nums = [str(res.sequence_number) for res in residues]
        print(f"  Sequence numbers: {', '.join(seq_nums)}")

    # Print detailed residue information
    print("\nDetailed Residue Information:")
    print("=" * 80)
    
    for res_name, residues in sorted(residues_by_name.items()):
        print(f"\n{res_name} Residues:")
        print("-" * 80)
        
        for residue in sorted(residues, key=lambda x: x.sequence_number):
            print(f"\nResidue {residue.name} {residue.sequence_number}:")
            print(f"  Chain ID: '{residue.chain_id}'")  # Print chain ID with quotes to see if it's empty
            print(f"  Number of atoms: {len(residue.atoms)}")
            
            # Calculate and print center of mass
            com = residue.center_of_mass()
            print(f"  Center of mass: ({com[0]:.3f}, {com[1]:.3f}, {com[2]:.3f})")
            
            # Print atom information for this residue
            print("\n  Atom Information:")
            print("  " + "-" * 90)
            print("  {:<6} {:<6} {:<10} {:<12} {:<12} {:<30}".format(
                "Name", "Type", "Topo Type", "Charge", "Mass", "Coordinates"))
            print("  " + "-" * 90)
            
            for atom in sorted(residue.atoms, key=lambda x: x.name):
                # Format charge with sign
                charge = f"{atom.topo_charge:+.3f}" if atom.topo_charge else "N/A"
                # Format mass
                mass = f"{atom.topo_mass:.3f}" if atom.topo_mass else "N/A"
                # Format coordinates
                coords = f"({atom.x:8.3f}, {atom.y:8.3f}, {atom.z:8.3f})"
                
                print("  {:<6} {:<6} {:<10} {:<12} {:<12} {:<30}".format(
                    atom.name,
                    atom.type if atom.type else "N/A",
                    atom.topo_type if atom.topo_type else "N/A",
                    charge,
                    mass,
                    coords
                ))

if __name__ == "__main__":
    # Example usage with test data
    test_dir = os.path.join(os.path.dirname(__file__), "..", "tests", "data")
    pdb_file = os.path.join(test_dir, "test.pdb")
    
    # Option 1: Using PSF files
    psf_files = [
        os.path.join(test_dir, "test_proa.psf"),
        os.path.join(test_dir, "mols", "benx.psf"),
        os.path.join(test_dir, "mols", "prpx.psf"),
        os.path.join(test_dir, "mols", "sol.psf")
    ]
    print("\nLoading structure with PSF files:")
    print_structure_info(pdb_file, psf_files=psf_files)
    
    print("\n" + "=" * 80 + "\n")
    
    # Option 2: Using TOP file
    top_file = os.path.join(test_dir, "test.top")
    print("\nLoading structure with TOP file:")
    print_structure_info(pdb_file, top_file=top_file) 