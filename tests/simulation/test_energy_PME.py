import pytest
import numpy as np
import math
import pygcmc
import os
from pygcmc import MCState, MCInfo, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo

# 直接复制create_nacl_crystal函数代码
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
        eps_na, math.sqrt(eps_na * eps_cl),
        math.sqrt(eps_na * eps_cl), eps_cl
    ]
    
    state.forcefield = ff
    
    # NaCl lattice constant (0.564 nm)
    a = 0.564  
    atoms = []
    residues = []
    
    # Create NaCl lattice
    print(f"\nCreating {n_cells}x{n_cells}x{n_cells} NaCl crystal...")
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
                
    print(f"Creation complete, added a total of {len(atoms)} atoms and {len(residues)} residues.")
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)
    
    return state

def test_pme_initialization():
    """
    Test PME parameter initialization
    """
    # Create a basic system
    box_size = 4.0
    n_cells = 2
    state = create_nacl_crystal(box_size, n_cells)
    
    # 设置系统的box尺寸
    box = [box_size, box_size, box_size]
    cutoff = box_size / 2.0
    
    # 明确设置PME参数，不使用自动调整
    alpha = 0.3
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    # 确保使用带有所有参数的显式调用
    print(f"Initializing PME with parameters: alpha={alpha}, mesh_size={mesh_size}, cutoff={cutoff}")
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    pygcmc.initializePMEParameters(cutoff, box, alpha, mesh_size, spline_order)
    
    # 输出参数便于调试
    print(f"PME parameters set. Ready to compute energy.")
    
    # Calculate energy
    elec_energy, vdw_energy, ewald_dict = pygcmc.computeSystemEnergyPME(state)
    
    # Check that energy components are reasonable
    assert ewald_dict["real_space"] != 0.0
    assert ewald_dict["reciprocal"] != 0.0
    assert ewald_dict["self"] != 0.0
    assert ewald_dict["total"] != 0.0
    
    print(f"PME energy components: real_space={ewald_dict['real_space']:.6f}, "
          f"reciprocal={ewald_dict['reciprocal']:.6f}, "
          f"self={ewald_dict['self']:.6f}, "
          f"total={ewald_dict['total']:.6f}")
    
    # Check that total energy is the sum of components
    assert abs(ewald_dict["total"] - (ewald_dict["real_space"] + 
                                    ewald_dict["reciprocal"] + 
                                    ewald_dict["self"] + vdw_energy)) < 1e-6

def test_pme_vs_ewald():
    """
    Compare PME and Ewald methods for a simple system
    """
    # Create a basic system
    box_size = 4.0
    n_cells = 2
    state = create_nacl_crystal(box_size, n_cells)
    
    # Set parameters for both methods
    alpha = 0.3
    kmax = [5, 5, 5]
    mesh_size = [16, 16, 16]
    
    box = [box_size, box_size, box_size]
    cutoff = box_size / 2.0
    
    # First calculate with Ewald
    pygcmc.setEwaldParameters(alpha, kmax)
    pygcmc.initializeEwaldParameters(cutoff, box, alpha)
    ewald_elec, ewald_vdw, ewald_dict = pygcmc.computeSystemEnergyEwald(state)
    
    ewald_real = ewald_dict["real_space"]
    ewald_recip = ewald_dict["reciprocal"]
    ewald_self = ewald_dict["self"]
    ewald_total = ewald_dict["total"]
    
    # Then calculate with PME
    pygcmc.setPMEParameters(alpha, mesh_size)
    pygcmc.initializePMEParameters(cutoff, box, alpha)
    pme_elec, pme_vdw, pme_dict = pygcmc.computeSystemEnergyPME(state)
    
    pme_real = pme_dict["real_space"]
    pme_recip = pme_dict["reciprocal"]
    pme_self = pme_dict["self"]
    pme_total = pme_dict["total"]
    
    # Print comparison
    print("Comparison of Ewald and PME energies:")
    print(f"Ewald: real={ewald_real:.6f}, recip={ewald_recip:.6f}, self={ewald_self:.6f}, total={ewald_total:.6f}")
    print(f"PME:   real={pme_real:.6f}, recip={pme_recip:.6f}, self={pme_self:.6f}, total={pme_total:.6f}")
    
    # Check that real-space and self energies are nearly identical
    assert abs(ewald_real - pme_real) < 1e-3 * abs(ewald_real)
    assert abs(ewald_self - pme_self) < 1e-6 * abs(ewald_self)
    
    # Reciprocal space may differ by a small percentage due to PME approximation
    assert abs(ewald_recip - pme_recip) < 0.05 * abs(ewald_recip)
    
    # Total energy should be close
    assert abs(ewald_total - pme_total) < 0.05 * abs(ewald_total)

def test_pme_spline_order():
    """
    Test the effect of different B-spline orders in PME
    """
    # Create a basic system
    box_size = 4.0
    n_cells = 2
    state = create_nacl_crystal(box_size, n_cells)
    
    box = [box_size, box_size, box_size]
    cutoff = box_size / 2.0
    alpha = 0.3
    mesh_size = [32, 32, 32]
    
    # Test different spline orders
    energies = []
    for spline_order in [2, 4, 6, 8]:
        # Initialize PME with specified spline order
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        pygcmc.initializePMEParameters(cutoff, box, alpha, mesh_size, spline_order)
        
        # Calculate energy
        _, _, pme_dict = pygcmc.computeSystemEnergyPME(state)
        total_energy = pme_dict["total"]
        energies.append(total_energy)
        
        print(f"Spline order {spline_order}: total energy = {total_energy:.6f}")
    
    # Higher order splines should converge to a stable value
    energy_diffs = [abs(energies[i] - energies[i-1]) for i in range(1, len(energies))]
    print(f"Energy differences between orders: {energy_diffs}")
    
    # Check that differences decrease with higher order
    for i in range(1, len(energy_diffs)):
        assert energy_diffs[i] <= energy_diffs[i-1] * 1.2  # Allow some numerical variance

def test_pme_error_tolerance():
    """
    Test PME error tolerance control
    """
    # Create a basic system
    box_size = 4.0
    n_cells = 2
    state = create_nacl_crystal(box_size, n_cells)
    
    box = [box_size, box_size, box_size]
    cutoff = box_size / 2.0
    
    # Test different error tolerances
    tolerances = [1e-3, 1e-4, 1e-5, 1e-6]
    energies = []
    
    for tol in tolerances:
        # Initialize with auto parameters and specified tolerance
        pygcmc.initializePMEParameters(cutoff, box, 0.0, None, 4, tol)
        
        # Calculate energy
        _, _, pme_dict = pygcmc.computeSystemEnergyPME(state)
        total_energy = pme_dict["total"]
        energies.append(total_energy)
        
        print(f"Error tolerance {tol}: total energy = {total_energy:.6f}")
    
    # Energies should converge with decreasing tolerance
    energy_diffs = [abs(energies[i] - energies[i-1]) for i in range(1, len(energies))]
    print(f"Energy differences between tolerances: {energy_diffs}")
    
    # Check that differences decrease with stricter tolerance
    for i in range(1, len(energy_diffs)):
        assert energy_diffs[i] <= energy_diffs[i-1] * 1.5  # Allow some numerical variance

def test_pme_movement_energy():
    """
    Test movement energy calculation with PME
    """
    # Create a basic system
    box_size = 4.0
    n_cells = 2
    state = create_nacl_crystal(box_size, n_cells)
    
    # Set PME parameters
    box = [box_size, box_size, box_size]
    cutoff = box_size / 2.0
    alpha = 0.3
    mesh_size = [32, 32, 32]
    
    pygcmc.setPMEParameters(alpha, mesh_size)
    pygcmc.initializePMEParameters(cutoff, box, alpha)
    
    # Set up movement residues
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 0
    movement_info.activeCount = 2  # Move 2 residues
    state.movementResidues.append(movement_info)
    
    # First calculate full system energy
    _, _, system_pme_dict = pygcmc.computeSystemEnergyPME(state)
    full_energy = system_pme_dict["total"]
    
    # Then calculate movement energy
    _, _, movement_pme_dict = pygcmc.computeMovementEnergyPME(state)
    movement_energy = movement_pme_dict["total"]
    
    print(f"Full system energy: {full_energy:.6f}")
    print(f"Movement energy: {movement_energy:.6f}")
    
    # Movement energy should be a subset of the full energy
    assert abs(movement_energy) <= abs(full_energy) * 1.1  # Allow some numerical variance

if __name__ == "__main__":
    test_pme_initialization()
    test_pme_vs_ewald()
    test_pme_spline_order()
    test_pme_error_tolerance()
    test_pme_movement_energy()
