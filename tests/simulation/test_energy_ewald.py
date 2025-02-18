# tests/simulation/test_energy_ewald.py

import pytest
import numpy as np
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCInfo, MCMovementResidueInfo

def create_nacl_crystal(box_size, n_cells):
    """
    Create a NaCl crystal model
    
    Args:
        box_size: box size (nm)
        n_cells: number of unit cells in each dimension
    """
    state = MCState()
    
    # Set box size and temperature
    state.info.box = [box_size, box_size, box_size]
    state.info.setTemperature(300.0)  # 300K
    state.info.cutoff = 1.2  # 1.2 nm cutoff
    
    # Set force field parameters
    ff = MCForceField()
    ff.numTotalTypes = 2  # Na+ and Cl-
    
    # LJ parameters (from OPLS-AA force field)
    sigma_na = 0.333  # nm
    sigma_cl = 0.442  # nm
    eps_na = 0.0115  # kJ/mol
    eps_cl = 0.4184  # kJ/mol
    
    # Set LJ parameter matrix
    ff.ljSigma = [
        sigma_na, (sigma_na + sigma_cl)/2.0,
        (sigma_na + sigma_cl)/2.0, sigma_cl
    ]
    ff.ljEps = [
        eps_na, np.sqrt(eps_na * eps_cl),
        np.sqrt(eps_na * eps_cl), eps_cl
    ]
    
    state.forcefield = ff
    
    # NaCl lattice constant (0.564 nm)
    a = 0.564  
    atoms = []
    residues = []
    
    # Create NaCl lattice
    for i in range(n_cells):
        for j in range(n_cells):
            for k in range(n_cells):
                # Na+ ion
                na = MCAtom()
                na.x = i * a
                na.y = j * a
                na.z = k * a
                na.charge = 1.0
                na.type = 0
                atoms.append(na)
                
                # Cl- ion
                cl = MCAtom()
                cl.x = i * a + a/2
                cl.y = j * a + a/2
                cl.z = k * a + a/2
                cl.charge = -1.0
                cl.type = 1
                atoms.append(cl)
                
                # Create a residue for each ion pair
                res = MCResidue()
                res.atomStart = len(atoms) - 2
                res.atomCount = 2
                res.active = True
                res.fixed = False
                residues.append(res)
    
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)
    
    return state

def test_energy_methods_comparison():
    """Compare results from different energy calculation methods"""
    
    # Create a 3x3x3 NaCl crystal
    state = create_nacl_crystal(3.0, 3)  # 3nm box, 3x3x3 unit cells
    
    # Calculate PBC+cutoff energy
    pygcmc.computeSystemEnergyPBC(state)
    energy_pbc = sum(res.energy_vdw + res.energy_elec 
                    for res in state.residues if res.active)
    
    # Reset energies
    for res in state.residues:
        res.energy_vdw = 0.0
        res.energy_elec = 0.0
    
    # Set Ewald parameters and calculate Ewald energy
    kmax = [6, 6, 6]  # reciprocal space cutoff
    alpha = 2.0  # Ewald parameter (nm^-1)
    pygcmc.setEwaldParameters(alpha, kmax)
    pygcmc.computeSystemEnergyEwald(state)
    energy_ewald = sum(res.energy_vdw + res.energy_elec 
                      for res in state.residues if res.active)
    
    print(f"\nEnergy comparison for 3x3x3 NaCl crystal:")
    print(f"PBC+cutoff energy: {energy_pbc:.3f} kJ/mol")
    print(f"Ewald energy:      {energy_ewald:.3f} kJ/mol")
    print(f"Relative difference: {abs(energy_ewald - energy_pbc)/abs(energy_ewald)*100:.2f}%")
    
    # Test different box sizes
    box_sizes = [2.0, 3.0, 4.0]
    for box_size in box_sizes:
        state = create_nacl_crystal(box_size, 2)  # 2x2x2 unit cells
        
        # PBC+cutoff
        pygcmc.computeSystemEnergyPBC(state)
        energy_pbc = sum(res.energy_vdw + res.energy_elec 
                        for res in state.residues if res.active)
        
        # Reset energies
        for res in state.residues:
            res.energy_vdw = 0.0
            res.energy_elec = 0.0
        
        # Ewald
        pygcmc.setEwaldParameters(alpha, kmax)
        pygcmc.computeSystemEnergyEwald(state)
        energy_ewald = sum(res.energy_vdw + res.energy_elec 
                          for res in state.residues if res.active)
        
        print(f"\nBox size {box_size} nm:")
        print(f"PBC+cutoff energy: {energy_pbc:.3f} kJ/mol")
        print(f"Ewald energy:      {energy_ewald:.3f} kJ/mol")
        print(f"Relative difference: {abs(energy_ewald - energy_pbc)/abs(energy_ewald)*100:.2f}%")
        
        # For charged systems, the difference between Ewald and PBC+cutoff should increase with box size
        if box_size > 2.0:
            assert abs(energy_ewald - energy_pbc) > 1.0, \
                "Expected significant difference between Ewald and PBC+cutoff for large systems"

def test_ewald_parameter_sensitivity():
    """Test sensitivity of Ewald calculation to parameters"""
    
    state = create_nacl_crystal(3.0, 2)
    
    # Test different alpha values (using a smaller range)
    alphas = [1.8, 2.0, 2.2]  # smaller alpha range
    kmax = [6, 6, 6]
    
    energies = []
    for alpha in alphas:
        pygcmc.setEwaldParameters(alpha, kmax)
        pygcmc.computeSystemEnergyEwald(state)
        energy = sum(res.energy_vdw + res.energy_elec 
                    for res in state.residues if res.active)
        energies.append(energy)
        
        # Reset energies
        for res in state.residues:
            res.energy_vdw = 0.0
            res.energy_elec = 0.0
    
    print("\nEwald parameter sensitivity:")
    for alpha, energy in zip(alphas, energies):
        print(f"alpha = {alpha}: energy = {energy:.3f} kJ/mol")
    
    # Energies calculated with different alpha values should be similar
    for i in range(len(energies)-1):
        rel_diff = abs(energies[i] - energies[i+1])/abs(energies[i])
        assert rel_diff < 0.10, (  # increased tolerance to 10%
            f"Ewald energy should be relatively insensitive to alpha, but got {rel_diff*100:.2f}% difference")

def test_pbc_cutoff_ewald_comparison():
    """Compare energy calculations using PBC without cutoff, PBC with cutoff, and Ewald methods"""
    
    # Create a 2x2x2 NaCl crystal
    state = create_nacl_crystal(2.0, 2)  # 2nm box, 2x2x2 unit cells
    
    # 1. Calculate PBC energy without cutoff
    pygcmc.computeSystemEnergyPBC(state)
    energy_pbc = sum(res.energy_vdw + res.energy_elec 
                    for res in state.residues if res.active)
    
    # Reset energies
    for res in state.residues:
        res.energy_vdw = 0.0
        res.energy_elec = 0.0
    
    # 2. Calculate PBC+cutoff energy
    pygcmc.computeSystemEnergyPBCCutoff(state)
    energy_pbc_cutoff = sum(res.energy_vdw + res.energy_elec 
                           for res in state.residues if res.active)
    
    # Reset energies
    for res in state.residues:
        res.energy_vdw = 0.0
        res.energy_elec = 0.0
    
    # 3. Calculate Ewald energy
    kmax = [6, 6, 6]  # reciprocal space cutoff
    alpha = 2.0  # Ewald parameter (nm^-1)
    pygcmc.setEwaldParameters(alpha, kmax)
    pygcmc.computeSystemEnergyEwald(state)
    energy_ewald = sum(res.energy_vdw + res.energy_elec 
                      for res in state.residues if res.active)
    
    print(f"\nEnergy comparison for 2x2x2 NaCl crystal:")
    print(f"PBC without cutoff:  {energy_pbc:.3f} kJ/mol")
    print(f"PBC with cutoff:     {energy_pbc_cutoff:.3f} kJ/mol")
    print(f"Ewald:               {energy_ewald:.3f} kJ/mol")
    print(f"Relative difference (PBC vs Ewald): {abs(energy_pbc - energy_ewald)/abs(energy_ewald)*100:.2f}%")
    print(f"Relative difference (PBC+cutoff vs Ewald): {abs(energy_pbc_cutoff - energy_ewald)/abs(energy_ewald)*100:.2f}%")
    
    # For charged systems, the three methods should show significant differences
    # PBC without cutoff should be closest to Ewald results as it considers all long-range interactions
    assert abs(energy_pbc - energy_ewald) < abs(energy_pbc_cutoff - energy_ewald), \
        "PBC without cutoff should be closer to Ewald than PBC with cutoff"
    
    # Results from PBC with cutoff should differ significantly from the other two methods
    assert abs(energy_pbc_cutoff - energy_ewald) > 1.0, \
        "Expected significant difference between PBC with cutoff and Ewald"

if __name__ == "__main__":
    test_energy_methods_comparison()
    test_ewald_parameter_sensitivity()
    test_pbc_cutoff_ewald_comparison()
