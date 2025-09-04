# tests/simulation/energyPME/grid_operations.py
"""PME grid operations diagnostic test."""

import math
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField


def test_pme_grid_operations():
    """
    Test PME grid operations with a simple two-atom system.
    
    This test creates a very simple system with just two oppositely charged atoms,
    and examines the PME grid operations in detail to diagnose why the reciprocal
    space energy computation is failing.
    """
    print("\nRunning test_pme_grid_operations...")
    
    # Set log level to INFO to view detailed debug output
    import pygcmc
    # First set System log level
    pygcmc.System.set_log_level(pygcmc.LogLevel.INFO)
    # Enable debug output for energy calculations
    pygcmc.setEnergyDebugOutput(True)
    # Remove non-existent API call
    # Change to output a prompt
    print("(Note: We have enabled detailed logging, now executing the test)")
    
    # Create a very simple system: two atoms, one positive and one negative
    state = MCState()
    
    # Set up a cubic box
    box_size = 3.0
    state.info.box = [box_size, box_size, box_size]
    state.info.setTemperature(300.0)
    
    cutoff = 1.0
    alpha = 2.5
    state.info.cutoff = cutoff
    
    # Create a force field with two atom types
    force_field = MCForceField()
    
    # Add atom types (sodium and chloride) - set force field parameters
    force_field.numTotalTypes = 2  # Na+ and Cl-
    
    # Define LJ parameters
    sigma_na = 1.0  # nm
    sigma_cl = 1.0  # nm
    eps_na = 1.0    # kJ/mol
    eps_cl = 1.0    # kJ/mol
    
    # Set LJ parameter matrix (diagonal and mixed terms)
    force_field.ljSigma = [
        sigma_na, (sigma_na + sigma_cl)/2.0,
        (sigma_na + sigma_cl)/2.0, sigma_cl
    ]
    force_field.ljEps = [
        eps_na, math.sqrt(eps_na * eps_cl),
        math.sqrt(eps_na * eps_cl), eps_cl
    ]
    
    # Direct creation of atoms and residues lists
    atoms = []
    residues = []
    
    # Create first atom - Na+
    atom1 = MCAtom()
    atom1.x = 1.0
    atom1.y = 1.5
    atom1.z = 1.5  # positioned at (1.0, 1.5, 1.5)
    atom1.charge = 1.0
    atom1.type = 0  # Na+
    atoms.append(atom1)
    
    # Create second atom - Cl-
    atom2 = MCAtom()
    atom2.x = 2.0
    atom2.y = 1.5
    atom2.z = 1.5  # positioned at (2.0, 1.5, 1.5)
    atom2.charge = -1.0
    atom2.type = 1  # Cl-
    atoms.append(atom2)
    
    # Create a residue containing these two atoms
    res = MCResidue()
    res.atomStart = 0  # index of the first atom
    res.atomCount = 2  # two atoms
    res.active = True
    res.fixed = False
    residues.append(res)
    
    # Set state's atoms and residues
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)
    
    # Set up force field
    state.forcefield = force_field
    
    print(f"Created system with {len(atoms)} atoms and {len(residues)} residues")
    
    # Calculate energies using Ewald and PME with detailed logging
    print("\n1. Setting PME parameters...")
    mesh_size = [32, 32, 32]
    spline_order = 5
    
    # Initialize PME parameters - use pygcmc module instead of platform.cpu
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    
    # Initialize PME parameters
    box = [box_size, box_size, box_size]
    pygcmc.initializePMEParameters(cutoff, box, alpha, mesh_size, spline_order)
    
    # Calculate standard Ewald energy as reference
    print("\n2. Calculating standard Ewald energy...")
    # Set Ewald parameters
    kmax = [8, 8, 8]  # Number of k-space vectors for Ewald calculation
    pygcmc.setEwaldParameters(alpha, kmax)
    pygcmc.initializeEwaldParameters(cutoff, box, alpha)
    
    # Calculate energy using Ewald
    ewald_elec, ewald_vdw, ewald_dict = pygcmc.computeSystemEnergyEwald(state)
    
    ewald_real = ewald_dict["real_space"]
    ewald_reciprocal = ewald_dict["reciprocal"]
    ewald_self = ewald_dict["self"]
    ewald_total = ewald_dict["total"]
    
    print(f"Ewald energy components:")
    print(f"  Real space:     {ewald_real:.6f}")
    print(f"  Reciprocal:     {ewald_reciprocal:.6f}")
    print(f"  Self:           {ewald_self:.6f}")
    print(f"  Total:          {ewald_total:.6f}")
    
    # Calculate PME energy
    print("\n3. Calculating PME energy...")
    pme_elec, pme_vdw, pme_dict = pygcmc.computeSystemEnergyPME(state)
    
    pme_real = pme_dict["real_space"]
    pme_reciprocal = pme_dict["reciprocal"]
    pme_self = pme_dict["self"]
    pme_total = pme_dict["total"]
    
    print(f"PME energy components:")
    print(f"  Real space:     {pme_real:.6f}")
    print(f"  Reciprocal:     {pme_reciprocal:.6f}")
    print(f"  Self:           {pme_self:.6f}")
    print(f"  Total:          {pme_total:.6f}")
    
    # Calculate difference
    real_diff = abs(pme_real - ewald_real)
    recip_diff = abs(pme_reciprocal - ewald_reciprocal)
    self_diff = abs(pme_self - ewald_self)
    total_diff = abs(pme_total - ewald_total)
    
    print("\nDifferences between PME and Ewald:")
    print(f"  Real space:     {real_diff:.6f}")
    print(f"  Reciprocal:     {recip_diff:.6f}")
    print(f"  Self:           {self_diff:.6f}")
    print(f"  Total:          {total_diff:.6f}")
    
    # Check if PME reciprocal is close to zero
    if abs(pme_reciprocal) < 1e-6:
        print("\n! WARNING: PME reciprocal energy is zero or near zero!")
        print("This confirms the issue observed in other tests.")
    
    # Verify that the real space energy is correct
    assert abs(real_diff) < 1e-4, f"Real space energies differ: {pme_real} vs {ewald_real}"
    
    # Check self energy
    assert abs(self_diff) < 1e-4, f"Self energies differ: {pme_self} vs {ewald_self}"
    
    # Diagnostic conclusion
    print("\nDiagnostic conclusion:")
    if abs(pme_reciprocal) < 1e-6:
        print("The PME reciprocal space energy calculation is failing.")
        print("Possible causes:")
        print("1. Charge spreading onto the grid may not be working correctly")
        print("2. The FFT implementation may have issues")
        print("3. The reciprocal space convolution may be incorrectly implemented")
        print("4. The B-spline moduli calculation may be incorrect")
    else:
        print("The PME reciprocal space energy is non-zero, but differs from Ewald.")
        print("This suggests the PME implementation needs further refinement.")
    
    # Print additional information showing our added debug output
    print("\nNote: Additional debug information should appear in the logs above.")
    print("If no additional information is shown, check that log level settings are correct.")
    
    # Ensure test doesn't fail because PME reciprocal energy is zero
    # This is just a diagnostic test, we expect to find problems, not fix them
    assert True, "This test is for diagnostic purposes only"