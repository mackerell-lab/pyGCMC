"""
Test PME with correctly structured NaCl crystal where each ion is its own residue
"""

import pygcmc
from pygcmc import MCState, MCInfo, MCAtom, MCResidue, MCForceField
import math


def create_nacl_crystal_correct(box_size, n_cells):
    """
    Create a NaCl crystal where each ion is its own residue
    
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
    print(f"\nCreating {n_cells}x{n_cells}x{n_cells} NaCl crystal (correct structure)...")
    atom_count = 0
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
                
                # Create a residue for Na+
                res_na = MCResidue()
                res_na.atomStart = atom_count
                res_na.atomCount = 1
                res_na.active = True
                res_na.fixed = False
                res_na.type = 0
                residues.append(res_na)
                atom_count += 1
                
                # Cl- ion
                cl = MCAtom()
                cl.x = i * a + a/2
                cl.y = j * a + a/2
                cl.z = k * a + a/2
                cl.charge = -1.0
                cl.type = 1
                atoms.append(cl)
                
                # Create a residue for Cl-
                res_cl = MCResidue()
                res_cl.atomStart = atom_count
                res_cl.atomCount = 1
                res_cl.active = True
                res_cl.fixed = False
                res_cl.type = 1
                residues.append(res_cl)
                atom_count += 1
                
    print(f"Creation complete, added {len(atoms)} atoms and {len(residues)} residues.")
    print(f"Each ion is in its own residue (correct structure)")
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)
    
    return state


def test_pme_vs_ewald_correct():
    """Test PME vs Ewald with correctly structured NaCl crystal"""
    
    # Create a basic system
    box_size = 4.0
    n_cells = 2
    state = create_nacl_crystal_correct(box_size, n_cells)
    
    # Set parameters for both methods
    alpha = 0.3
    kmax = [5, 5, 5]
    mesh_size = [32, 32, 32]
    spline_order = 4
    
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
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    pygcmc.initializePMEParameters(cutoff, box, alpha, mesh_size, spline_order)
    
    pme_elec, pme_vdw, pme_dict = pygcmc.computeSystemEnergyPME(state)
    
    pme_real = pme_dict["real_space"]
    pme_recip = pme_dict["reciprocal"]
    pme_self = pme_dict["self"]
    pme_total = pme_dict["total"]
    
    # Calculate relative errors
    real_rel_error = abs(ewald_real - pme_real) / abs(ewald_real) if abs(ewald_real) > 1e-6 else 0
    recip_rel_error = abs(ewald_recip - pme_recip) / abs(ewald_recip) if abs(ewald_recip) > 1e-6 else 0
    self_rel_error = abs(ewald_self - pme_self) / abs(ewald_self) if abs(ewald_self) > 1e-6 else 0
    total_rel_error = abs(ewald_total - pme_total) / abs(ewald_total) if abs(ewald_total) > 1e-6 else 0
    
    # Print comparison
    print("\nEnergy comparison (correctly structured NaCl):")
    print(f"              Real Space      Reciprocal     Self           Total")
    print(f"Ewald:        {ewald_real:10.6f}  {ewald_recip:10.6f}  {ewald_self:10.6f}  {ewald_total:10.6f}")
    print(f"PME:          {pme_real:10.6f}  {pme_recip:10.6f}  {pme_self:10.6f}  {pme_total:10.6f}")
    print(f"Rel. Error:   {real_rel_error:10.6f}  {recip_rel_error:10.6f}  {self_rel_error:10.6f}  {total_rel_error:10.6f}")
    
    # Analysis
    print("\nAnalysis:")
    print(f"Real space difference: {abs(ewald_real - pme_real):.6f} kJ/mol")
    print(f"Total energy difference: {abs(ewald_total - pme_total):.6f} kJ/mol")
    
    # Check if PME and Ewald agree when structure is correct
    assert total_rel_error < 1e-5, f"PME and Ewald should agree for correctly structured system, but error is {total_rel_error:.6%}"
    print("\n✓ PME and Ewald agree for correctly structured NaCl crystal!")


if __name__ == "__main__":
    test_pme_vs_ewald_correct()