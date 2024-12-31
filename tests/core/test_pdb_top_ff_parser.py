# tests/core/test_pdb_top_ff_parser.py

import os
import sys
import math  # 导入 math 模块以处理 'math.isnan'
import pyGCMC_bindings as core

def print_atom_info(atoms):
    print("\nDetailed Atom Information:")
    print("-" * 120)
    print(f"{'Residue':6} {'Seq':4} {'Name':4} {'Type':6} {'TopoType':8} {'Charge':8} {'Mass':8} {'Epsilon':8} {'Rmin':8}")
    print("-" * 120)
    
    for atom in atoms:
        epsilon = atom.forcefield_epsilon if not math.isnan(atom.forcefield_epsilon) else 'nan'
        rmin = atom.forcefield_rmin if not math.isnan(atom.forcefield_rmin) else 'nan'
        print(f"{atom.residue:6} {atom.sequence:4d} {atom.name:4} {atom.type:6} {atom.topo_type:8} "
              f"{atom.topo_charge:8.3f} {atom.topo_mass:8.3f} "
              f"{epsilon:8} {rmin:8}")

def debug_atom_info(atom):
    print(f"\nDebug info for atom {atom.name} in residue {atom.residue} {atom.sequence}:")
    print(f"  Basic info:")
    print(f"    serial: {atom.serial}")
    print(f"    name: {atom.name}")
    print(f"    residue: {atom.residue}")
    print(f"    sequence: {atom.sequence}")
    print(f"    position: ({atom.x}, {atom.y}, {atom.z})")
    print(f"  Topology info:")
    print(f"    type: {atom.type}")
    print(f"    topo_type: {atom.topo_type}")
    print(f"    topo_charge: {atom.topo_charge}")
    print(f"    topo_mass: {atom.topo_mass}")
    print(f"  Force field info:")
    print(f"    forcefield_epsilon: {atom.forcefield_epsilon}")
    print(f"    forcefield_rmin: {atom.forcefield_rmin}")
    print(f"  Status:")
    print(f"    has_topology_info: {atom.has_topology_info()}")
    print(f"    has_forcefield_info: {atom.has_forcefield_info()}")

def main():
    # Get the directory containing test data
    test_data_dir = os.path.join(os.path.dirname(__file__), "..", "data")
    
    # Initialize parsers
    pdb_parser = core.PDBParser()
    top_parser = core.TopParser()
    ff_parser_prot = core.FFParser()
    ff_parser_cgenff = core.FFParser()
    ff_parser_water = core.FFParser()
    ff_parser_silcs = core.FFParser()

    # Parse PDB file
    pdb_file = os.path.join(test_data_dir, "test.pdb")
    coords, residues = core.PDBParser.parse(pdb_file)
    if not residues:
        print("Error: PDB file contains no residues")
        return 1

    # Extract atoms from residues
    atoms = []
    for residue in residues:
        atoms.extend(residue.atoms)

    # Parse topology file with includes
    top_file = os.path.join(test_data_dir, "test.top")
    if not top_parser.parse_with_includes(top_file):
        print("Error: Failed to parse topology file")
        return 1

    # Parse force field files
    prot_prm_file = os.path.join(test_data_dir, "par_all36m_prot.prm")
    if not ff_parser_prot.parse(prot_prm_file):
        print("Error: Failed to parse protein force field file")
        return 1

    cgenff_prm_file = os.path.join(test_data_dir, "par_all36_cgenff.prm")
    if not ff_parser_cgenff.parse(cgenff_prm_file):
        print("Error: Failed to parse general force field file")
        return 1

    water_ions_str = os.path.join(test_data_dir, "toppar_water_ions.str")
    if not ff_parser_water.parse(water_ions_str):
        print("Error: Failed to parse water and ions force field file")
        return 1

    silcs_str = os.path.join(test_data_dir, "silcs.str")
    if not ff_parser_silcs.parse(silcs_str):
        print("Error: Failed to parse SILCS force field file")
        return 1

    print(f"Reading files from {test_data_dir}")
    print(f"PDB file: {pdb_file}")
    print(f"TOP file: {top_file}")
    print(f"Protein parameter file: {prot_prm_file}")
    print(f"CGenFF parameter file: {cgenff_prm_file}")
    print(f"Water parameter file: {water_ions_str}")
    print(f"SILCS parameter file: {silcs_str}")
    
    try:
        # 1. Parse PDB file
        coords, residues = core.PDBParser.parse(pdb_file)
        print(f"\nFound {len(residues)} residues in PDB file")
        
        # Extract atoms from residues
        atoms = []
        for residue in residues:
            atoms.extend(residue.atoms)
        print(f"Total number of atoms: {len(atoms)}")
        
        # Print debug info for first atom of each residue type
        seen_residues = set()
        for atom in atoms:
            if atom.residue not in seen_residues:
                debug_atom_info(atom)
                seen_residues.add(atom.residue)
        
        # 2. Parse topology file with includes
        if not top_parser.parse_with_includes(top_file):
            print("Failed to parse topology file")
            return 1
        print("\nSuccessfully parsed topology file")
        
        # 3. Parse force field files
        if not ff_parser_prot.parse(prot_prm_file):
            print("Failed to parse protein force field file")
            return 1
        if not ff_parser_cgenff.parse(cgenff_prm_file):
            print("Failed to parse CGenFF force field file")
            return 1
        if not ff_parser_water.parse(water_ions_str):
            print("Failed to parse water force field file")
            return 1
        if not ff_parser_silcs.parse(silcs_str):
            print("Failed to parse SILCS force field file")
            return 1
        print("Successfully parsed force field files")
        
        # Add debug print before topology update
        print("\nBefore topology update:")
        for atom in atoms[:5]:  # Print first 5 atoms for brevity
            print(f"Atom {atom.name} ({atom.residue} {atom.sequence}): type={atom.type}")
        
        # 4. Update atoms with topology information
        # 使用新的 update_pdb_atoms 方法
        topo_updated = top_parser.update_pdb_atoms(atoms)
        print(f"\nUpdated {topo_updated} atoms with topology information")
        
        # Add debug print after topology update
        print("\nAfter topology update:")
        for atom in atoms[:5]:
            print(f"Atom {atom.name} ({atom.residue} {atom.sequence}): "
                  f"type={atom.type}, topo_type={atom.topo_type}, "
                  f"topo_charge={atom.topo_charge}, topo_mass={atom.topo_mass}")
        
        # Add debug print for force field update
        print("\nBefore force field update:")
        for atom in atoms[:5]:
            print(f"Atom {atom.name} ({atom.residue} {atom.sequence}): "
                  f"type={atom.type}, topo_type={atom.topo_type}, "
                  f"epsilon={atom.forcefield_epsilon}, rmin={atom.forcefield_rmin}")
        
        # 5. Update atoms with force field parameters
        # 使用新的 update_pdb_atoms 方法
        ff_updated_prot = ff_parser_prot.update_pdb_atoms(atoms)
        ff_updated_gen = ff_parser_cgenff.update_pdb_atoms(atoms)
        ff_updated_water = ff_parser_water.update_pdb_atoms(atoms)
        ff_updated_silcs = ff_parser_silcs.update_pdb_atoms(atoms)
        print(f"Updated {ff_updated_prot} atoms with protein force field parameters")
        print(f"Updated {ff_updated_gen} atoms with CGenFF force field parameters")
        print(f"Updated {ff_updated_water} atoms with water force field parameters")
        print(f"Updated {ff_updated_silcs} atoms with SILCS force field parameters")
        
        # Add debug print after force field update
        print("\nAfter force field update:")
        for atom in atoms[:5]:
            epsilon = atom.forcefield_epsilon if not math.isnan(atom.forcefield_epsilon) else 'nan'
            rmin = atom.forcefield_rmin if not math.isnan(atom.forcefield_rmin) else 'nan'
            print(f"Atom {atom.name} ({atom.residue} {atom.sequence}): "
                  f"type={atom.type}, topo_type={atom.topo_type}, "
                  f"epsilon={epsilon}, rmin={rmin}")
        
        # Print detailed information for each atom
        print_atom_info(atoms)
        
    except Exception as e:
        print(f"Error: {e}")
        return 1
        
    return 0

if __name__ == "__main__":
    sys.exit(main()) 
