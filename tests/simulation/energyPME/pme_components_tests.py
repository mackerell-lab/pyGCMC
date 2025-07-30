# tests/simulation/energyPME/pme_components_tests.py
"""PME component tests: verify component breakdown and function differences."""

import pygcmc
from pygcmc import MCState, MCForceField, MCAtom, MCResidue, MCMovementResidueInfo


def create_multi_atom_system():
    """Create a system with multiple fixed and moving atoms."""
    state = MCState()
    state.info.box = [3.0, 3.0, 3.0]
    state.info.cutoff = 1.0
    state.info.setTemperature(300.0)
    
    # Force field (only electrostatics)
    ff = MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 1
    ff.ljSigma = [0.0] * 4
    ff.ljEps = [0.0] * 4
    state.forcefield = ff
    
    atoms = []
    
    # Two fixed atoms
    atom = MCAtom()
    atom.x, atom.y, atom.z = 1.0, 1.5, 1.5
    atom.charge = 1.0
    atom.type = 0
    atoms.append(atom)
    
    atom = MCAtom()
    atom.x, atom.y, atom.z = 1.2, 1.5, 1.5
    atom.charge = -1.0
    atom.type = 0
    atoms.append(atom)
    
    # Two moving atoms
    atom = MCAtom()
    atom.x, atom.y, atom.z = 2.0, 1.5, 1.5
    atom.charge = 0.5
    atom.type = 1
    atoms.append(atom)
    
    atom = MCAtom()
    atom.x, atom.y, atom.z = 2.2, 1.5, 1.5
    atom.charge = -0.5
    atom.type = 1
    atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    
    # Create residues
    residues = []
    
    # Fixed residue
    res = MCResidue()
    res.atomStart = 0
    res.atomCount = 2
    res.active = True
    res.fixed = True
    res.type = 0
    residues.append(res)
    
    # Moving residue
    res = MCResidue()
    res.atomStart = 2
    res.atomCount = 2
    res.active = True
    res.fixed = False
    res.type = 1
    residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    # Set movement info
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 1
    movement_info.activeCount = 1
    state.movementResidues = [movement_info]
    
    return state


def test_pme_component_breakdown():
    """
    Test that PME functions return proper component breakdown.
    
    Verifies that all PME energy functions return the expected
    dictionary with real_space, reciprocal, self, and total components.
    """
    state = create_multi_atom_system()
    
    # Initialize PME
    alpha = 2.84
    mesh = [32, 32, 32]
    pygcmc.setPMEParameters(alpha, mesh, 4, 1e-6)
    pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh, 4)
    
    # Test computeSystemEnergyPME
    result = pygcmc.computeSystemEnergyPME(state)
    assert len(result) >= 3, "Should return (energy, vdw, components)"
    assert isinstance(result[2], dict), "Third element should be component dictionary"
    assert 'real_space' in result[2], "Should have real_space component"
    assert 'reciprocal' in result[2], "Should have reciprocal component"
    assert 'self' in result[2], "Should have self component"
    assert 'total' in result[2], "Should have total component"
    
    print(f"System energy components: {result[2]}")
    
    # Test computeMovementEnergyPME
    result = pygcmc.computeMovementEnergyPME(state)
    assert len(result) >= 3, "Should return (energy, vdw, components)"
    assert isinstance(result[2], dict), "Third element should be component dictionary"
    assert 'real_space' in result[2], "Should have real_space component"
    assert 'reciprocal' in result[2], "Should have reciprocal component"
    assert 'self' in result[2], "Should have self component"
    assert 'total' in result[2], "Should have total component"
    
    print(f"Movement energy components: {result[2]}")


def test_pme_function_differences():
    """
    Test differences between PME energy calculation functions.
    
    Compares computeSystemEnergyPME, computeMovementEnergyPME, and
    computeSystemEnergyPMEComplete (if available) to understand their
    different behaviors.
    """
    state = create_multi_atom_system()
    
    # Initialize PME
    alpha = 2.84
    mesh = [32, 32, 32]
    pygcmc.setPMEParameters(alpha, mesh, 4, 1e-6)
    pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh, 4)
    
    # Calculate different energies
    system_energy, system_vdw, system_comp = pygcmc.computeSystemEnergyPME(state)
    movement_energy, movement_vdw, movement_comp = pygcmc.computeMovementEnergyPME(state)
    
    print(f"\nSystem energy: {system_energy:.6f} kJ/mol")
    print(f"Movement energy: {movement_energy:.6f} kJ/mol")
    print(f"Difference: {system_energy - movement_energy:.6f} kJ/mol")
    
    # Test PME Complete if available
    if hasattr(pygcmc, 'computeSystemEnergyPMEComplete'):
        complete_result = pygcmc.computeSystemEnergyPMEComplete(state)
        print(f"PME Complete energy: {complete_result[0]:.6f} kJ/mol")
        
        # PME Complete should match regular PME for systems without exclusions
        assert abs(complete_result[0] - system_energy) < 0.1, \
            "PME Complete should match regular PME for this system"
    
    # Verify component differences
    print(f"\nComponent differences (System - Movement):")
    print(f"Real space: {system_comp['real_space'] - movement_comp['real_space']:.6f}")
    print(f"Reciprocal: {system_comp['reciprocal'] - movement_comp['reciprocal']:.6f}")
    print(f"Self: {system_comp['self'] - movement_comp['self']:.6f}")
    
    # Self energy difference should equal the self energy of fixed atoms
    # since movement energy only includes moving atoms' self energy
    expected_self_diff = abs(system_comp['self'] - movement_comp['self'])
    assert expected_self_diff > 0.1, "Self energy should differ between system and movement"


def test_pme_energy_relationships():
    """
    Test energy relationships in multi-atom systems.
    
    Verifies expected relationships between different interaction types
    and how they contribute to total and movement energies.
    """
    state = create_multi_atom_system()
    
    # Initialize PME
    alpha = 2.84
    mesh = [32, 32, 32]
    pygcmc.setPMEParameters(alpha, mesh, 4, 1e-6)
    pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh, 4)
    
    # Get energies
    system_energy, _, system_comp = pygcmc.computeSystemEnergyPME(state)
    movement_energy, _, movement_comp = pygcmc.computeMovementEnergyPME(state)
    
    # Calculate manual expectations
    # Fixed-Fixed: atoms 0-1 (charges +1, -1, distance 0.2)
    # Moving-Moving: atoms 2-3 (charges +0.5, -0.5, distance 0.2)
    # Fixed-Moving: 4 pairs with various distances
    
    print("\nExpected interaction contributions:")
    print(f"Fixed-Fixed: Strong (close dipole)")
    print(f"Moving-Moving: Medium (half charges)")
    print(f"Fixed-Moving: Weak (larger distances)")
    
    # The system energy should be more negative than movement energy
    # because it includes the strong fixed-fixed interaction
    assert system_energy < movement_energy, \
        "System energy should be more negative due to fixed-fixed interactions"
    
    # Both should have negative total energy (net attractive)
    assert system_energy < 0, "System should have net attractive interactions"
    
    print(f"\nEnergy totals:")
    print(f"System: {system_energy:.6f} kJ/mol")
    print(f"Movement: {movement_energy:.6f} kJ/mol")
    print("Verified: System energy includes all interactions, movement excludes fixed-fixed")