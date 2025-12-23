# tests/system/MonteCarloSystem/forcefield_tests.py

import pytest
from .helpers import *

def test_initialize_force_field(molecular_system, charmm_ff):
    print("\n=== Starting test_initialize_force_field ===")
    
    # Create Monte Carlo system
    mc_system = pygcmc.MonteCarloSystem()
    print("Created MonteCarloSystem")
    
    # Set reasonable max capacity
    info = pygcmc.MCInfo()
    info.max_residues = 1000
    info.max_atoms = 10000
    mc_system.initialize(info)
    print("Initialized MonteCarloSystem with max capacity")
    
    # Print molecular system info
    print(f"\nMolecular system info:")
    print(f"Number of residues: {molecular_system.get_num_residues()}")
    print(f"Number of atoms: {molecular_system.get_num_atoms()}")
    print(f"Box dimensions: {molecular_system.boxDimensions}")
    
    # Initialize with molecular system
    mc_system.initialize_from_molecular(molecular_system)
    print("\nInitialized from molecular system")
    
    # Add movement molecules
    # Load small molecules
    movement_molecules = []
    
    # Load benzene
    benx_pdb = os.path.join(TEST_DATA_DIR, "mols", "benx.pdb")
    benx_psf = os.path.join(TEST_DATA_DIR, "mols", "benx.psf")
    benx_structure = pygcmc.PDBParser.parse_file(benx_pdb)
    benx_topology = pygcmc.PSFParser.parse_file(benx_psf)
    benx_mol = pygcmc.MolecularSystem()
    benx = benx_mol.combine(benx_structure, benx_topology)
    movement_molecules.append(pygcmc.MovementMolecularInfo(benx, 10))  # 10 copies
    print("\nLoaded benzene as movement molecule")
    
    # Load imidazole
    imia_pdb = os.path.join(TEST_DATA_DIR, "mols", "imia.pdb")
    imia_psf = os.path.join(TEST_DATA_DIR, "mols", "imia.psf")
    imia_structure = pygcmc.PDBParser.parse_file(imia_pdb)
    imia_topology = pygcmc.PSFParser.parse_file(imia_psf)
    imia_mol = pygcmc.MolecularSystem()
    imia = imia_mol.combine(imia_structure, imia_topology)
    movement_molecules.append(pygcmc.MovementMolecularInfo(imia, 10))  # 10 copies
    print("Loaded imidazole as movement molecule")
    
    # Load water
    sol_pdb = os.path.join(TEST_DATA_DIR, "mols", "sol.pdb")
    sol_top = os.path.join(TEST_DATA_DIR, "mols", "sol.itp")
    sol_structure = pygcmc.PDBParser.parse_file(sol_pdb)
    sol_topology = pygcmc.TOPParser.parse_file(sol_top)
    sol_mol = pygcmc.MolecularSystem()
    sol = sol_mol.combine(sol_structure, sol_topology)
    movement_molecules.append(pygcmc.MovementMolecularInfo(sol, 20))  # 20 copies for water
    print("Loaded water as movement molecule")
    
    # Add movement molecules to system
    mc_system.add_movement_molecules(movement_molecules)
    print("\nAdded all movement molecules to system")
    
    # Get atom types before force field initialization
    atom_types = mc_system.get_type_maps()
    print(f"\nNumber of atom types: {len(atom_types.atomTypes)}")
    print("Atom types in the system:", atom_types.atomTypes)
    
    # Get state and check movement atom types
    state = mc_system.get_state()
    print("\nMovement atom types info:")
    print(f"Number of movement atom types: {state.numMovementAtomTypes}")
    print("Movement atom types indices:", state.movementAtomTypes)
    print("Movement residues info:")
    for info in state.movementResidues:
        print(f"  {info.resName}: start={info.startIndex}, active={info.activeCount}, total={info.totalCount}")
    
    # Print available LJ parameters in the force field
    print("\nAvailable LJ parameters in force field:")
    lj_params = charmm_ff.lj_params
    print(f"Number of LJ parameters: {len(lj_params)}")
    for atom_type, params in lj_params.items():
        print(f"{atom_type}: epsilon={params.epsilon:.4f} kcal/mol, rmin_half={params.rmin_half:.4f} Å")
    
    # Check if all atom types have LJ parameters
    missing_types = []
    for type_name in atom_types.atomTypes:
        try:
            params = charmm_ff.get_lj_params(type_name)
            print(f"Found LJ params for {type_name}: epsilon={params.epsilon:.4f} kcal/mol, rmin_half={params.rmin_half:.4f} Å")
        except Exception as e:
            print(f"Failed to get LJ params for {type_name}: {str(e)}")
            missing_types.append(type_name)
    
    if missing_types:
        print("\nMissing LJ parameters for atom types:", missing_types)
        pytest.fail(f"Missing LJ parameters for atom types: {missing_types}")
    
    print("\nAll atom types have LJ parameters, proceeding with force field initialization")
    
    # Initialize force field
    try:
        mc_system.initialize_force_field(charmm_ff)
        print("Force field initialization successful")
    except Exception as e:
        print(f"Force field initialization failed: {str(e)}")
        raise
    
    # Get state after force field initialization
    state = mc_system.get_state()
    print("\nGot state after force field initialization")
    
    # Check force field parameters
    print("\nForce field parameters:")
    print(f"numTotalTypes: {state.forcefield.numTotalTypes}")
    print(f"numMovementTypes: {state.forcefield.numMovementTypes}")
    print(f"ljSigma size: {len(state.forcefield.ljSigma)}")
    print(f"ljEps size: {len(state.forcefield.ljEps)}")
    
    # Basic assertions
    assert state.forcefield.numTotalTypes == len(atom_types.atomTypes)
    assert state.forcefield.numMovementTypes == state.numMovementAtomTypes
    # Modified assertion for array size, now should be numTotalTypes * numTotalTypes
    expected_size = state.forcefield.numTotalTypes * state.forcefield.numTotalTypes
    assert len(state.forcefield.ljSigma) == expected_size, \
        f"Expected ljSigma size {expected_size}, got {len(state.forcefield.ljSigma)}"
    assert len(state.forcefield.ljEps) == expected_size, \
        f"Expected ljEps size {expected_size}, got {len(state.forcefield.ljEps)}"
    
    # Unit conversion constants
    ANGSTROM_TO_NM = 0.1  # 1 Å = 0.1 nm
    KCAL_TO_KJ = 4.184    # 1 kcal/mol = 4.184 kJ/mol
    
    # Modified construction logic for type pairs, now checking all type pairs
    atom_type_pairs = []
    all_type_indices = range(len(atom_types.atomTypes))
    
    # Check all possible type pairs
    for type1_idx in all_type_indices:
        type1 = atom_types.atomTypes[type1_idx]
        for type2_idx in all_type_indices:
            type2 = atom_types.atomTypes[type2_idx]
            atom_type_pairs.append((type1, type2))
    
    print(f"\nConstructed {len(atom_type_pairs)} type pairs to check")
    print("First few pairs:", atom_type_pairs[:5])

    # Modified parameter index calculation
    for type1, type2 in atom_type_pairs:
        # Get type indices in the force field
        type1_idx = atom_types.get_or_add_type(type1)
        type2_idx = atom_types.get_or_add_type(type2)
        
        # Calculate pair index in the force field arrays
        pairIdx = type1_idx * state.forcefield.numTotalTypes + type2_idx
        
        # Get actual parameters from force field
        actual_eps = state.forcefield.ljEps[pairIdx]
        actual_sigma = state.forcefield.ljSigma[pairIdx]
        
        # Get NBFIX parameters if they exist
        nbfix_result = charmm_ff.get_nbfix(type1, type2)
        
        if nbfix_result[1]:  # If NBFIX exists
            # Convert NBFIX epsilon from kcal/mol to kJ/mol
            expected_eps = nbfix_result[0] * KCAL_TO_KJ
            # Convert NBFIX Rmin to sigma in nm
            expected_sigma = nbfix_result[1] / 2.0**(1.0/6.0) * ANGSTROM_TO_NM
            
            print(f"NBFIX {type1}-{type2}:")
            print(f"  epsilon: expected {expected_eps:.6f}, got {actual_eps:.6f} kJ/mol")
            print(f"  sigma: expected {expected_sigma:.6f}, got {actual_sigma:.6f} nm")
            
            assert abs(actual_eps - expected_eps) < 1e-6, \
                f"NBFIX epsilon mismatch for {type1}-{type2}: expected {expected_eps}, got {actual_eps} kJ/mol"
            assert abs(actual_sigma - expected_sigma) < 1e-6, \
                f"NBFIX sigma mismatch for {type1}-{type2}: expected {expected_sigma}, got {actual_sigma} nm"
        else:  # No NBFIX, use combination rules
            lj1 = charmm_ff.get_lj_params(type1)
            lj2 = charmm_ff.get_lj_params(type2)
            # Convert Rmin/2 from Å to nm and then to sigma
            sigma1 = 2 * (lj1.rmin_half / math.pow(2.0, 1.0/6.0)) * ANGSTROM_TO_NM
            sigma2 = 2 * (lj2.rmin_half / math.pow(2.0, 1.0/6.0)) * ANGSTROM_TO_NM
            
            # Calculate expected combined parameters
            # First combine in original units (kcal/mol), then convert to kJ/mol
            expected_eps = math.sqrt(lj1.epsilon * lj2.epsilon) * KCAL_TO_KJ
            # Combine sigma values (already converted to nm)
            expected_sigma = 0.5 * (sigma1 + sigma2)
            
            print(f"Combined LJ {type1}-{type2}:")
            print(f"  Expected: eps={expected_eps:.4f} kJ/mol, sigma={expected_sigma:.4f} nm")
            print(f"  Actual:   eps={actual_eps:.4f} kJ/mol, sigma={actual_sigma:.4f} nm")
            
            assert abs(actual_eps - expected_eps) < 1e-6, \
                f"Combined epsilon mismatch for {type1}-{type2}: expected {expected_eps}, got {actual_eps} kJ/mol"
            assert abs(actual_sigma - expected_sigma) < 1e-6, \
                f"Combined sigma mismatch for {type1}-{type2}: expected {expected_sigma}, got {actual_sigma} nm"
    
    print("\n=== test_initialize_force_field completed successfully ===")
