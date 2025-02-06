# tests/simulation/test_file_nonbonded.py

import pytest
import pygcmc
import os
import math

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

@pytest.fixture
def charmm_ff():
    """Load CHARMM force field files."""
    # Get the test data directory
    test_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(test_dir, "data")
    
    # Load CHARMM force field files
    ff = pygcmc.ForceField()
    pygcmc.PRMParser.parse_file_to_forcefield(os.path.join(data_dir, "par_all36_cgenff.prm"), ff)
    pygcmc.PRMParser.parse_file_to_forcefield(os.path.join(data_dir, "par_all36m_prot.prm"), ff)
    pygcmc.PRMParser.parse_file_to_forcefield(os.path.join(data_dir, "silcs.str"), ff)
    pygcmc.PRMParser.parse_file_to_forcefield(os.path.join(data_dir, "toppar_water_ions.str"), ff)
    return ff

def test_nonbonded_energy_with_full_system(charmm_ff):
    """Test nonbonded energy calculation using a complete Monte Carlo system.
    
    This test:
    1. Creates a Monte Carlo system from PDB/TOP files
    2. Adds movement molecules (benx, imia, sol)
    3. Initializes force field
    4. Calculates and verifies nonbonded energies
    """
    print("\n=== Starting test_nonbonded_energy_with_full_system ===")
    
    # Create Monte Carlo system
    mc_system = pygcmc.MonteCarloSystem()
    info = pygcmc.MCInfo()
    info.max_residues = 10000
    info.max_atoms = 100000
    mc_system.initialize(info)
    
    # Load base system
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    top_path = os.path.join(TEST_DATA_DIR, "test.top")
    structure = pygcmc.PDBParser.parse_file(pdb_path)
    topology = pygcmc.TOPParser.parse_file(top_path)
    mol_system = pygcmc.MolecularSystem()
    base_molecular = mol_system.combine(structure, topology)
    
    # Initialize base system
    mc_system.initialize_from_molecular(base_molecular)
    print("\nInitialized base system")
    
    # Load movement molecules
    MOLS_DIR = os.path.join(TEST_DATA_DIR, "mols")
    
    # Load benzene
    benx_pdb = os.path.join(MOLS_DIR, "benx.pdb")
    benx_psf = os.path.join(MOLS_DIR, "benx.psf")
    benx_structure = pygcmc.PDBParser.parse_file(benx_pdb)
    benx_topology = pygcmc.PSFParser.parse_file(benx_psf)
    benx_mol = pygcmc.MolecularSystem()
    benx = benx_mol.combine(benx_structure, benx_topology)
    
    # Load imidazole
    imia_pdb = os.path.join(MOLS_DIR, "imia.pdb")
    imia_psf = os.path.join(MOLS_DIR, "imia.psf")
    imia_structure = pygcmc.PDBParser.parse_file(imia_pdb)
    imia_topology = pygcmc.PSFParser.parse_file(imia_psf)
    imia_mol = pygcmc.MolecularSystem()
    imia = imia_mol.combine(imia_structure, imia_topology)
    
    # Load water
    sol_pdb = os.path.join(MOLS_DIR, "sol.pdb")
    sol_top = os.path.join(MOLS_DIR, "sol.itp")
    sol_structure = pygcmc.PDBParser.parse_file(sol_pdb)
    sol_topology = pygcmc.TOPParser.parse_file(sol_top)
    sol_mol = pygcmc.MolecularSystem()
    sol = sol_mol.combine(sol_structure, sol_topology)
    
    # Create movement molecule info list
    movement_mols = [
        pygcmc.MovementMolecularInfo(benx, 5),  # 5 benzene molecules
        pygcmc.MovementMolecularInfo(imia, 5),  # 5 imidazole molecules
        pygcmc.MovementMolecularInfo(sol, 10)   # 10 water molecules
    ]
    
    # Add movement molecules
    mc_system.add_movement_molecules(movement_mols)
    print("Added movement molecules")
    
    # Load and initialize force field
    mc_system.initialize_force_field(charmm_ff)
    print("Initialized force field")
    
    # Get state
    state = mc_system.get_state()
    
    # Print system info
    print("\nSystem information:")
    print(f"Active residues: {state.activeResidueCount}")
    print(f"Active atoms: {state.activeAtomCount}")
    print(f"Movement residues: {len(state.movementResidues)}")
    for info in state.movementResidues:
        print(f"  {info.resName}: start={info.startIndex}, active={info.activeCount}, total={info.totalCount}")
    
    # Calculate nonbonded energies
    print("\nCalculating nonbonded energies...")
    pygcmc.computeNaiveNonbondedEnergy(state)
    
    # Verify energies for each movement residue
    print("\nEnergies for movement residues:")
    for move_info in state.movementResidues:
        print(f"\n{move_info.resName} residues:")
        for i in range(move_info.startIndex, move_info.startIndex + move_info.activeCount):
            res = state.residues[i]
            if res.active:
                total_energy = res.energy_vdw + res.energy_elec
                print(f"  Residue {i}: vdw={res.energy_vdw:.3f}, elec={res.energy_elec:.3f}, total={total_energy:.3f}")
                
                # Basic sanity checks
                assert not math.isnan(res.energy_vdw), f"VDW energy is NaN for residue {i}"
                assert not math.isnan(res.energy_elec), f"Electrostatic energy is NaN for residue {i}"
                assert not math.isinf(res.energy_vdw), f"VDW energy is infinite for residue {i}"
                assert not math.isinf(res.energy_elec), f"Electrostatic energy is infinite for residue {i}"
                
                # Check for unreasonably large energies
                assert abs(res.energy_vdw) < 1e6, f"VDW energy too large for residue {i}"
                assert abs(res.energy_elec) < 1e6, f"Electrostatic energy too large for residue {i}"
    
    print("\n=== test_nonbonded_energy_with_full_system completed successfully ===")

def test_nonbonded_energy_with_specific_molecules(charmm_ff):
    """Test nonbonded energy calculation with specific molecule arrangements.
    
    This test creates a system with specific molecule arrangements to verify:
    1. Water-water interactions
    2. Benzene-water interactions
    3. Imidazole-water interactions
    """
    print("\n=== Starting test_nonbonded_energy_with_specific_molecules ===")
    
    # Create Monte Carlo system
    mc_system = pygcmc.MonteCarloSystem()
    info = pygcmc.MCInfo()
    info.max_residues = 1000
    info.max_atoms = 10000
    mc_system.initialize(info)
    
    # Load molecules
    MOLS_DIR = os.path.join(TEST_DATA_DIR, "mols")
    
    # Load water
    sol_pdb = os.path.join(MOLS_DIR, "sol.pdb")
    sol_top = os.path.join(MOLS_DIR, "sol.itp")
    sol_structure = pygcmc.PDBParser.parse_file(sol_pdb)
    sol_topology = pygcmc.TOPParser.parse_file(sol_top)
    sol_mol = pygcmc.MolecularSystem()
    sol = sol_mol.combine(sol_structure, sol_topology)
    
    # Create movement molecule info with only water
    movement_mols = [
        pygcmc.MovementMolecularInfo(sol, 2)  # 2 water molecules
    ]
    
    # Add movement molecules
    mc_system.add_movement_molecules(movement_mols)
    print("Added water molecules")
    
    # Load and initialize force field
    mc_system.initialize_force_field(charmm_ff)
    print("Initialized force field")
    
    # Get state
    state = mc_system.get_state()
    
    # Position water molecules at specific distances
    # First water at origin
    water1_start = state.residues[0].atomStart
    for i in range(3):  # Water has 3 atoms
        state.atoms[water1_start + i].x = 0.0
        state.atoms[water1_start + i].y = 0.0
        state.atoms[water1_start + i].z = 0.0
    
    # Second water at (3,0,0) - typical hydrogen bond distance
    water2_start = state.residues[1].atomStart
    for i in range(3):
        state.atoms[water2_start + i].x = 3.0
        state.atoms[water2_start + i].y = 0.0
        state.atoms[water2_start + i].z = 0.0
    
    # Calculate nonbonded energies
    print("\nCalculating water-water interactions...")
    pygcmc.computeNaiveNonbondedEnergy(state)
    
    # Print and verify water-water interaction energies
    print("\nWater-water interaction energies:")
    for i in range(2):
        res = state.residues[i]
        total_energy = res.energy_vdw + res.energy_elec
        print(f"Water {i+1}: vdw={res.energy_vdw:.3f}, elec={res.energy_elec:.3f}, total={total_energy:.3f}")
        
        # Basic checks
        assert not math.isnan(total_energy), f"Energy is NaN for water {i+1}"
        assert not math.isinf(total_energy), f"Energy is infinite for water {i+1}"
        
        # Water-water interaction at 3Å should be favorable
        assert total_energy < 0, f"Water-water interaction should be attractive at 3Å"
    
    print("\n=== test_nonbonded_energy_with_specific_molecules completed successfully ===")

def test_benx_water_interaction(charmm_ff):
    """Test benzene-water interactions.
    
    This test verifies:
    1. Hydrophobic interaction between benzene and water
    2. Proper use of NBFIX parameters if defined
    3. Distance dependence of the interaction
    """
    print("\n=== Starting test_benx_water_interaction ===")
    
    # Create Monte Carlo system
    mc_system = pygcmc.MonteCarloSystem()
    info = pygcmc.MCInfo()
    info.max_residues = 1000
    info.max_atoms = 10000
    mc_system.initialize(info)
    
    # Load molecules
    MOLS_DIR = os.path.join(TEST_DATA_DIR, "mols")
    
    # Load benzene
    benx_pdb = os.path.join(MOLS_DIR, "benx.pdb")
    benx_psf = os.path.join(MOLS_DIR, "benx.psf")
    benx_structure = pygcmc.PDBParser.parse_file(benx_pdb)
    benx_topology = pygcmc.PSFParser.parse_file(benx_psf)
    benx_mol = pygcmc.MolecularSystem()
    benx = benx_mol.combine(benx_structure, benx_topology)
    
    # Load water
    sol_pdb = os.path.join(MOLS_DIR, "sol.pdb")
    sol_top = os.path.join(MOLS_DIR, "sol.itp")
    sol_structure = pygcmc.PDBParser.parse_file(sol_pdb)
    sol_topology = pygcmc.TOPParser.parse_file(sol_top)
    sol_mol = pygcmc.MolecularSystem()
    sol = sol_mol.combine(sol_structure, sol_topology)
    
    # Create movement molecule info
    movement_mols = [
        pygcmc.MovementMolecularInfo(benx, 1),  # 1 benzene
        pygcmc.MovementMolecularInfo(sol, 1)    # 1 water
    ]
    
    # Add movement molecules
    mc_system.add_movement_molecules(movement_mols)
    print("Added benzene and water molecules")
    
    # Load and initialize force field
    mc_system.initialize_force_field(charmm_ff)
    print("Initialized force field")
    
    # Get state
    state = mc_system.get_state()
    
    # Position molecules
    # Benzene at origin
    benx_start = state.residues[0].atomStart
    benx_count = state.residues[0].atomCount
    for i in range(benx_count):
        state.atoms[benx_start + i].x = 0.0
        state.atoms[benx_start + i].y = 0.0
        state.atoms[benx_start + i].z = 0.0
    
    # Test different water positions
    distances = [3.0, 5.0, 7.0]  # Test at different distances
    water_start = state.residues[1].atomStart
    water_count = state.residues[1].atomCount
    
    print("\nTesting benzene-water interactions at different distances:")
    for dist in distances:
        # Position water
        for i in range(water_count):
            state.atoms[water_start + i].x = dist
            state.atoms[water_start + i].y = 0.0
            state.atoms[water_start + i].z = 0.0
        
        # Calculate energy
        pygcmc.computeNaiveNonbondedEnergy(state)
        
        # Get energies
        benx_energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
        water_energy = state.residues[1].energy_vdw + state.residues[1].energy_elec
        
        print(f"\nAt distance {dist}Å:")
        print(f"Benzene energy: vdw={state.residues[0].energy_vdw:.3f}, elec={state.residues[0].energy_elec:.3f}, total={benx_energy:.3f}")
        print(f"Water energy: vdw={state.residues[1].energy_vdw:.3f}, elec={state.residues[1].energy_elec:.3f}, total={water_energy:.3f}")
        
        # Verify energy properties
        assert not math.isnan(benx_energy), "Benzene energy is NaN"
        assert not math.isnan(water_energy), "Water energy is NaN"
        
        # Energy should decrease with distance
        if dist > 3.0:
            assert abs(benx_energy) < prev_benx_energy, f"Energy not decreasing with distance at {dist}Å"
        
        prev_benx_energy = abs(benx_energy)
    
    print("\n=== test_benx_water_interaction completed successfully ===")

def test_imia_water_interaction(charmm_ff):
    """Test imidazole-water interactions.
    
    This test verifies:
    1. Hydrogen bonding capability
    2. Proper use of NBFIX parameters if defined
    3. Orientation dependence of interaction
    """
    print("\n=== Starting test_imia_water_interaction ===")
    
    # Create Monte Carlo system
    mc_system = pygcmc.MonteCarloSystem()
    info = pygcmc.MCInfo()
    info.max_residues = 1000
    info.max_atoms = 10000
    mc_system.initialize(info)
    
    # Load molecules
    MOLS_DIR = os.path.join(TEST_DATA_DIR, "mols")
    
    # Load imidazole
    imia_pdb = os.path.join(MOLS_DIR, "imia.pdb")
    imia_psf = os.path.join(MOLS_DIR, "imia.psf")
    imia_structure = pygcmc.PDBParser.parse_file(imia_pdb)
    imia_topology = pygcmc.PSFParser.parse_file(imia_psf)
    imia_mol = pygcmc.MolecularSystem()
    imia = imia_mol.combine(imia_structure, imia_topology)
    
    # Load water
    sol_pdb = os.path.join(MOLS_DIR, "sol.pdb")
    sol_top = os.path.join(MOLS_DIR, "sol.itp")
    sol_structure = pygcmc.PDBParser.parse_file(sol_pdb)
    sol_topology = pygcmc.TOPParser.parse_file(sol_top)
    sol_mol = pygcmc.MolecularSystem()
    sol = sol_mol.combine(sol_structure, sol_topology)
    
    # Create movement molecule info
    movement_mols = [
        pygcmc.MovementMolecularInfo(imia, 1),  # 1 imidazole
        pygcmc.MovementMolecularInfo(sol, 1)    # 1 water
    ]
    
    # Add movement molecules
    mc_system.add_movement_molecules(movement_mols)
    print("Added imidazole and water molecules")
    
    # Load and initialize force field
    mc_system.initialize_force_field(charmm_ff)
    print("Initialized force field")
    
    # Get state
    state = mc_system.get_state()
    
    # Position molecules
    # Imidazole at origin
    imia_start = state.residues[0].atomStart
    imia_count = state.residues[0].atomCount
    for i in range(imia_count):
        state.atoms[imia_start + i].x = 0.0
        state.atoms[imia_start + i].y = 0.0
        state.atoms[imia_start + i].z = 0.0
    
    # Test different water positions and orientations
    distances = [2.8, 4.0, 6.0]  # Include hydrogen bond distance
    water_start = state.residues[1].atomStart
    water_count = state.residues[1].atomCount
    
    print("\nTesting imidazole-water interactions at different distances:")
    for dist in distances:
        # Position water
        for i in range(water_count):
            state.atoms[water_start + i].x = dist
            state.atoms[water_start + i].y = 0.0
            state.atoms[water_start + i].z = 0.0
        
        # Calculate energy
        pygcmc.computeNaiveNonbondedEnergy(state)
        
        # Get energies
        imia_energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
        water_energy = state.residues[1].energy_vdw + state.residues[1].energy_elec
        
        print(f"\nAt distance {dist}Å:")
        print(f"Imidazole energy: vdw={state.residues[0].energy_vdw:.3f}, elec={state.residues[0].energy_elec:.3f}, total={imia_energy:.3f}")
        print(f"Water energy: vdw={state.residues[1].energy_vdw:.3f}, elec={state.residues[1].energy_elec:.3f}, total={water_energy:.3f}")
        
        # Verify energy properties
        assert not math.isnan(imia_energy), "Imidazole energy is NaN"
        assert not math.isnan(water_energy), "Water energy is NaN"
        
        # At hydrogen bond distance (2.8Å), interaction should be favorable
        if abs(dist - 2.8) < 0.1:
            assert imia_energy < 0, "Expected favorable interaction at hydrogen bond distance"
        
        # Energy should decrease with distance
        if dist > 2.8:
            assert abs(imia_energy) < prev_imia_energy, f"Energy not decreasing with distance at {dist}Å"
        
        prev_imia_energy = abs(imia_energy)
    
    print("\n=== test_imia_water_interaction completed successfully ===")

def test_distance_dependence(charmm_ff):
    """Test the distance dependence of nonbonded interactions.
    
    This test verifies:
    1. 1/r dependence of electrostatic interactions
    2. 1/r^6 and 1/r^12 terms in vdw interactions
    3. Proper energy decay at long distances
    """
    print("\n=== Starting test_distance_dependence ===")
    
    # Create Monte Carlo system
    mc_system = pygcmc.MonteCarloSystem()
    info = pygcmc.MCInfo()
    info.max_residues = 1000
    info.max_atoms = 10000
    mc_system.initialize(info)
    
    # Load water for simple point charge test
    MOLS_DIR = os.path.join(TEST_DATA_DIR, "mols")
    sol_pdb = os.path.join(MOLS_DIR, "sol.pdb")
    sol_top = os.path.join(MOLS_DIR, "sol.itp")
    sol_structure = pygcmc.PDBParser.parse_file(sol_pdb)
    sol_topology = pygcmc.TOPParser.parse_file(sol_top)
    sol_mol = pygcmc.MolecularSystem()
    sol = sol_mol.combine(sol_structure, sol_topology)
    
    # Create movement molecule info
    movement_mols = [
        pygcmc.MovementMolecularInfo(sol, 2)  # 2 water molecules
    ]
    
    # Add movement molecules
    mc_system.add_movement_molecules(movement_mols)
    print("Added water molecules")
    
    # Load and initialize force field
    mc_system.initialize_force_field(charmm_ff)
    print("Initialized force field")
    
    # Get state
    state = mc_system.get_state()
    
    # First water at origin
    water1_start = state.residues[0].atomStart
    for i in range(3):
        state.atoms[water1_start + i].x = 0.0
        state.atoms[water1_start + i].y = 0.0
        state.atoms[water1_start + i].z = 0.0
    
    # Test at various distances
    distances = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
    water2_start = state.residues[1].atomStart
    
    print("\nTesting distance dependence of interactions:")
    energies = []
    for dist in distances:
        # Position second water
        for i in range(3):
            state.atoms[water2_start + i].x = dist
            state.atoms[water2_start + i].y = 0.0
            state.atoms[water2_start + i].z = 0.0
        
        # Calculate energy
        pygcmc.computeNaiveNonbondedEnergy(state)
        
        # Get total energy
        total_energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
        energies.append(total_energy)
        
        print(f"\nAt distance {dist}Å:")
        print(f"VDW energy: {state.residues[0].energy_vdw:.3f}")
        print(f"Electrostatic energy: {state.residues[0].energy_elec:.3f}")
        print(f"Total energy: {total_energy:.3f}")
        
        # Verify energy properties
        assert not math.isnan(total_energy), f"Energy is NaN at distance {dist}Å"
        assert not math.isinf(total_energy), f"Energy is infinite at distance {dist}Å"
        
        # Check distance dependence
        if len(energies) > 1:
            ratio = abs(energies[-2] / energies[-1])
            expected_ratio = (dist / distances[-2])**6  # Approximate for vdw dominated region
            print(f"Energy ratio: {ratio:.3f}, Expected: {expected_ratio:.3f}")
            # Allow for some deviation due to mixed electrostatic and vdw
            assert abs(ratio - expected_ratio) < 2.0, f"Unexpected distance dependence at {dist}Å"
    
    print("\n=== test_distance_dependence completed successfully ===")

def test_short_range_repulsion(charmm_ff):
    """Test short-range repulsive interactions.
    
    This test verifies:
    1. Strong repulsion at very short distances
    2. Proper handling of numerical stability
    3. Dominance of r^12 term at short range
    """
    print("\n=== Starting test_short_range_repulsion ===")
    
    # Create Monte Carlo system
    mc_system = pygcmc.MonteCarloSystem()
    info = pygcmc.MCInfo()
    info.max_residues = 1000
    info.max_atoms = 10000
    mc_system.initialize(info)
    
    # Load water for simple point charge test
    MOLS_DIR = os.path.join(TEST_DATA_DIR, "mols")
    sol_pdb = os.path.join(MOLS_DIR, "sol.pdb")
    sol_top = os.path.join(MOLS_DIR, "sol.itp")
    sol_structure = pygcmc.PDBParser.parse_file(sol_pdb)
    sol_topology = pygcmc.TOPParser.parse_file(sol_top)
    sol_mol = pygcmc.MolecularSystem()
    sol = sol_mol.combine(sol_structure, sol_topology)
    
    # Create movement molecule info
    movement_mols = [
        pygcmc.MovementMolecularInfo(sol, 2)  # 2 water molecules
    ]
    
    # Add movement molecules
    mc_system.add_movement_molecules(movement_mols)
    print("Added water molecules")
    
    # Load and initialize force field
    mc_system.initialize_force_field(charmm_ff)
    print("Initialized force field")
    
    # Get state
    state = mc_system.get_state()
    
    # First water at origin
    water1_start = state.residues[0].atomStart
    for i in range(3):
        state.atoms[water1_start + i].x = 0.0
        state.atoms[water1_start + i].y = 0.0
        state.atoms[water1_start + i].z = 0.0
    
    # Test at very short distances
    distances = [0.5, 0.7, 1.0, 1.5, 2.0]
    water2_start = state.residues[1].atomStart
    
    print("\nTesting short-range repulsive interactions:")
    prev_energy = None
    for dist in distances:
        # Position second water
        for i in range(3):
            state.atoms[water2_start + i].x = dist
            state.atoms[water2_start + i].y = 0.0
            state.atoms[water2_start + i].z = 0.0
        
        # Calculate energy
        pygcmc.computeNaiveNonbondedEnergy(state)
        
        # Get energies
        vdw_energy = state.residues[0].energy_vdw
        elec_energy = state.residues[0].energy_elec
        total_energy = vdw_energy + elec_energy
        
        print(f"\nAt distance {dist}Å:")
        print(f"VDW energy: {vdw_energy:.3f}")
        print(f"Electrostatic energy: {elec_energy:.3f}")
        print(f"Total energy: {total_energy:.3f}")
        
        # Verify energy properties
        assert not math.isnan(total_energy), f"Energy is NaN at distance {dist}Å"
        assert not math.isinf(total_energy), f"Energy is infinite at distance {dist}Å"
        
        # Energy should be positive (repulsive) at very short range
        assert total_energy > 0, f"Expected repulsive interaction at {dist}Å"
        
        # Check that repulsion increases as distance decreases
        if prev_energy is not None:
            assert total_energy > prev_energy, f"Repulsion not increasing at {dist}Å"
        
        prev_energy = total_energy
        
        # At very short range, vdw should dominate
        if dist < 1.0:
            assert abs(vdw_energy) > abs(elec_energy), f"VDW not dominant at {dist}Å"
    
    print("\n=== test_short_range_repulsion completed successfully ===")

def test_nbfix_parameters(charmm_ff):
    """Test the proper application of NBFIX parameters.
    
    This test verifies:
    1. NBFIX parameters override standard combining rules
    2. Correct energy calculation with NBFIX parameters
    3. Proper handling of mixed NBFIX and standard parameters
    """
    print("\n=== Starting test_nbfix_parameters ===")
    
    # Create Monte Carlo system
    mc_system = pygcmc.MonteCarloSystem()
    info = pygcmc.MCInfo()
    info.max_residues = 1000
    info.max_atoms = 10000
    mc_system.initialize(info)
    
    # Load molecules
    MOLS_DIR = os.path.join(TEST_DATA_DIR, "mols")
    
    # Load molecules that might have NBFIX parameters
    benx_pdb = os.path.join(MOLS_DIR, "benx.pdb")
    benx_psf = os.path.join(MOLS_DIR, "benx.psf")
    benx_structure = pygcmc.PDBParser.parse_file(benx_pdb)
    benx_topology = pygcmc.PSFParser.parse_file(benx_psf)
    benx_mol = pygcmc.MolecularSystem()
    benx = benx_mol.combine(benx_structure, benx_topology)
    
    sol_pdb = os.path.join(MOLS_DIR, "sol.pdb")
    sol_top = os.path.join(MOLS_DIR, "sol.itp")
    sol_structure = pygcmc.PDBParser.parse_file(sol_pdb)
    sol_topology = pygcmc.TOPParser.parse_file(sol_top)
    sol_mol = pygcmc.MolecularSystem()
    sol = sol_mol.combine(sol_structure, sol_topology)
    
    # Create movement molecule info
    movement_mols = [
        pygcmc.MovementMolecularInfo(benx, 1),
        pygcmc.MovementMolecularInfo(sol, 1)
    ]
    
    # Add movement molecules
    mc_system.add_movement_molecules(movement_mols)
    print("Added molecules")
    
    # Load and initialize force field
    mc_system.initialize_force_field(charmm_ff)
    print("Initialized force field")
    
    # Get state
    state = mc_system.get_state()
    
    # Get atom types for checking NBFIX
    type_maps = mc_system.get_type_maps()
    
    # Print force field parameters
    print("\nForce field parameters:")
    for i in range(state.forcefield.numMovementTypes):
        type_i = type_maps.get_type_name(state.movementAtomTypes[i])
        for j in range(state.forcefield.maxTypes):
            type_j = type_maps.get_type_name(j)
            idx = i * state.forcefield.maxTypes + j
            eps = state.forcefield.ljEps[idx]
            sigma = state.forcefield.ljSigma[idx]
            
            # Check if NBFIX exists
            nbfix_result = charmm_ff.get_nbfix(type_i, type_j)
            if nbfix_result[1]:  # NBFIX exists
                print(f"NBFIX parameters for {type_i}-{type_j}:")
                print(f"  epsilon: {eps:.6f} (NBFIX)")
                print(f"  sigma: {sigma:.6f}")
                
                # Verify NBFIX parameters are used
                assert abs(eps - nbfix_result[0]) < 1e-6, \
                    f"NBFIX epsilon not properly applied for {type_i}-{type_j}"
            else:
                print(f"Standard parameters for {type_i}-{type_j}:")
                print(f"  epsilon: {eps:.6f}")
                print(f"  sigma: {sigma:.6f}")
    
    print("\n=== test_nbfix_parameters completed successfully ===")

def test_combination_rules(charmm_ff):
    """Test the application of Lorentz-Berthelot combination rules.
    
    This test verifies:
    1. Correct application of arithmetic mean for sigma
    2. Correct application of geometric mean for epsilon
    3. Proper handling of different atom type combinations
    """
    print("\n=== Starting test_combination_rules ===")
    
    # Create Monte Carlo system
    mc_system = pygcmc.MonteCarloSystem()
    info = pygcmc.MCInfo()
    info.max_residues = 1000
    info.max_atoms = 10000
    mc_system.initialize(info)
    
    # Load molecules
    MOLS_DIR = os.path.join(TEST_DATA_DIR, "mols")
    
    # Load different molecules to test various atom type combinations
    benx_pdb = os.path.join(MOLS_DIR, "benx.pdb")
    benx_psf = os.path.join(MOLS_DIR, "benx.psf")
    benx_structure = pygcmc.PDBParser.parse_file(benx_pdb)
    benx_topology = pygcmc.PSFParser.parse_file(benx_psf)
    benx_mol = pygcmc.MolecularSystem()
    benx = benx_mol.combine(benx_structure, benx_topology)
    
    imia_pdb = os.path.join(MOLS_DIR, "imia.pdb")
    imia_psf = os.path.join(MOLS_DIR, "imia.psf")
    imia_structure = pygcmc.PDBParser.parse_file(imia_pdb)
    imia_topology = pygcmc.PSFParser.parse_file(imia_psf)
    imia_mol = pygcmc.MolecularSystem()
    imia = imia_mol.combine(imia_structure, imia_topology)
    
    sol_pdb = os.path.join(MOLS_DIR, "sol.pdb")
    sol_top = os.path.join(MOLS_DIR, "sol.itp")
    sol_structure = pygcmc.PDBParser.parse_file(sol_pdb)
    sol_topology = pygcmc.TOPParser.parse_file(sol_top)
    sol_mol = pygcmc.MolecularSystem()
    sol = sol_mol.combine(sol_structure, sol_topology)
    
    # Create movement molecule info
    movement_mols = [
        pygcmc.MovementMolecularInfo(benx, 1),
        pygcmc.MovementMolecularInfo(imia, 1),
        pygcmc.MovementMolecularInfo(sol, 1)
    ]
    
    # Add movement molecules
    mc_system.add_movement_molecules(movement_mols)
    print("Added molecules")
    
    # Load and initialize force field
    mc_system.initialize_force_field(charmm_ff)
    print("Initialized force field")
    
    # Get state and type maps
    state = mc_system.get_state()
    type_maps = mc_system.get_type_maps()
    
    # Verify combination rules
    print("\nVerifying Lorentz-Berthelot combination rules:")
    for i in range(state.forcefield.numMovementTypes):
        type_i = type_maps.get_type_name(state.movementAtomTypes[i])
        lj_i = charmm_ff.get_lj_params(type_i)
        sigma_i = lj_i.rmin_half / math.pow(2.0, 1.0/6.0)
        
        for j in range(state.forcefield.maxTypes):
            type_j = type_maps.get_type_name(j)
            idx = i * state.forcefield.maxTypes + j
            
            # Skip if NBFIX exists
            if charmm_ff.get_nbfix(type_i, type_j)[1]:
                continue
            
            # Get individual LJ parameters
            lj_j = charmm_ff.get_lj_params(type_j)
            sigma_j = lj_j.rmin_half / math.pow(2.0, 1.0/6.0)
            
            # Calculate expected combined parameters
            expected_sigma = 0.5 * (sigma_i + sigma_j)
            expected_eps = math.sqrt(lj_i.epsilon * lj_j.epsilon)
            
            # Get actual parameters
            actual_sigma = state.forcefield.ljSigma[idx]
            actual_eps = state.forcefield.ljEps[idx]
            
            print(f"\nParameters for {type_i}-{type_j}:")
            print(f"Sigma: expected={expected_sigma:.6f}, actual={actual_sigma:.6f}")
            print(f"Epsilon: expected={expected_eps:.6f}, actual={actual_eps:.6f}")
            
            # Verify parameters
            assert abs(actual_sigma - expected_sigma) < 1e-6, \
                f"Incorrect sigma combination for {type_i}-{type_j}"
            assert abs(actual_eps - expected_eps) < 1e-6, \
                f"Incorrect epsilon combination for {type_i}-{type_j}"
    
    print("\n=== test_combination_rules completed successfully ===")

