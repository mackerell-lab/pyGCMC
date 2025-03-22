import unittest
import numpy as np
import math
import pygcmc
from pygcmc import MCState, MCInfo, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo

# Direct copy of create_nacl_crystal function from test_energy_PME.py
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

class TestEnergyPGP(unittest.TestCase):
    """
    Test the PGP (Pair-Grid PME) implementation
    """
    
    def setUp(self):
        """
        Set up a NaCl crystal system for testing.
        """
        # Create a basic NaCl crystal for testing
        self.box_size = 2.82  # nm, approximately 28.2 Å
        self.n_cells = 2      # 2x2x2 supercell
        self.system = create_nacl_crystal(self.box_size, self.n_cells)
        
        # Set cutoff within 1/2 box size to avoid problems with PBCs
        self.cutoff = 1.0   # nm
        self.pair_cutoff = 0.5  # nm
        
        # Temperature in K
        self.temp = 300.0
        self.system.set_temperature(self.temp)
        
    def test_pgp_initialization(self):
        """
        Test initialization of PGP parameters
        """
        # Initialize with default parameters
        alpha = 0.29  # 1/nm
        mesh_size = [32, 32, 32]
        pair_grid_size = [16, 16, 16]
        spline_order = 4
        tolerance = 1e-5
        
        # Initialize PGP parameters
        pygcmc.initializePGPParameters(
            cutoff=self.cutoff,
            pair_cutoff=self.pair_cutoff,
            box=[self.box_size, self.box_size, self.box_size],
            alpha=alpha,
            meshSize=mesh_size,
            pairGridSize=pair_grid_size,
            splineOrder=spline_order,
            tolerance=tolerance
        )
        
        # Calculate energy with PGP
        elec, vdw, energy_dict = pygcmc.computeSystemEnergyPGP(self.system.state)
        
        # Basic validation - energy should be finite and reasonable
        self.assertTrue(np.isfinite(elec))
        self.assertTrue(np.isfinite(vdw))
        self.assertTrue(np.isfinite(energy_dict["total"]))
        
        # Print the energy components for debugging
        print(f"PGP Energy: elec={elec:.6f}, vdw={vdw:.6f}, total={energy_dict['total']:.6f}")
        print(f"  Components: real={energy_dict['real_space']:.6f}, recip={energy_dict['reciprocal']:.6f}, self={energy_dict['self']:.6f}")
        
    def test_pgp_vs_pme(self):
        """
        Compare PGP results with PME for validation
        """
        # Parameters
        alpha = 0.29  # 1/nm
        mesh_size = [32, 32, 32]
        pair_grid_size = [16, 16, 16]
        spline_order = 4
        tolerance = 1e-5
        
        # Initialize PME parameters
        pygcmc.initializePMEParameters(
            cutoff=self.cutoff,
            box=[self.box_size, self.box_size, self.box_size],
            alpha=alpha,
            meshSize=mesh_size,
            splineOrder=spline_order,
            tolerance=tolerance
        )
        
        # Calculate energy with PME
        elec_pme, vdw_pme, energy_dict_pme = pygcmc.computeSystemEnergyPME(self.system.state)
        
        # Reset system energy
        self.system.reset_energy()
        
        # Initialize PGP parameters with same alpha and mesh
        pygcmc.initializePGPParameters(
            cutoff=self.cutoff,
            pair_cutoff=self.pair_cutoff,
            box=[self.box_size, self.box_size, self.box_size],
            alpha=alpha,
            meshSize=mesh_size,
            pairGridSize=pair_grid_size,
            splineOrder=spline_order,
            tolerance=tolerance
        )
        
        # Calculate energy with PGP
        elec_pgp, vdw_pgp, energy_dict_pgp = pygcmc.computeSystemEnergyPGP(self.system.state)
        
        # Compare results - for now, PGP should give similar results to PME
        # since the pair-grid part is just a placeholder
        # In future, this test will be updated once the pair-grid is fully implemented
        print(f"PME Energy: elec={elec_pme:.6f}, vdw={vdw_pme:.6f}, total={energy_dict_pme['total']:.6f}")
        print(f"PGP Energy: elec={elec_pgp:.6f}, vdw={vdw_pgp:.6f}, total={energy_dict_pgp['total']:.6f}")
        
        # For now, values should be similar since we're using PME as base implementation
        self.assertAlmostEqual(elec_pme, elec_pgp, delta=abs(elec_pme)*0.01)  # 1% tolerance
        self.assertAlmostEqual(vdw_pme, vdw_pgp, delta=abs(vdw_pme)*0.01)     # 1% tolerance
        
    def test_pgp_movement_energy(self):
        """
        Test PGP energy calculation for movement residues
        """
        # Create a separate test system
        system = create_nacl_crystal(self.box_size, self.n_cells)
        system.set_temperature(self.temp)
        
        # Initialize PGP parameters
        pygcmc.initializePGPParameters(
            cutoff=self.cutoff,
            pair_cutoff=self.pair_cutoff,
            box=[self.box_size, self.box_size, self.box_size],
            alpha=0.29,
            meshSize=[32, 32, 32],
            pairGridSize=[16, 16, 16],
            splineOrder=4,
            tolerance=1e-5
        )
        
        # Get total system energy
        elec_total, vdw_total, energy_dict_total = pygcmc.computeSystemEnergyPGP(system.state)
        
        # Create movement info for a subset of atoms (1/8 of the system)
        n_residues = len(system.state.residues)
        movement_residues = list(range(0, n_residues, 8))
        system.set_movement_residues(movement_residues)
        
        # Calculate movement energy
        elec_move, vdw_move, energy_dict_move = pygcmc.computeMovementEnergyPGP(system.state)
        
        # Validate - movement energy should be smaller than total energy
        print(f"Total Energy:   elec={elec_total:.6f}, vdw={vdw_total:.6f}, total={energy_dict_total['total']:.6f}")
        print(f"Movement Energy: elec={elec_move:.6f}, vdw={vdw_move:.6f}, total={energy_dict_move['total']:.6f}")
        
        # Movement energy should be non-zero but less than total energy
        self.assertGreater(abs(elec_move), 0.0)
        self.assertLess(abs(elec_move), abs(elec_total) * 1.001)  # With a small margin for numerical errors
        
    def test_pgp_pair_grid_parameters(self):
        """
        Test different pair grid parameters in PGP
        """
        # Try different pair grid sizes
        pair_grid_sizes = [
            [8, 8, 8],
            [16, 16, 16],
            [32, 32, 32]
        ]
        
        # Common parameters
        alpha = 0.29
        mesh_size = [32, 32, 32]
        spline_order = 4
        
        # Test different pair grid sizes
        results = []
        for pg_size in pair_grid_sizes:
            # Initialize PGP
            pygcmc.initializePGPParameters(
                cutoff=self.cutoff,
                pair_cutoff=self.pair_cutoff,
                box=[self.box_size, self.box_size, self.box_size],
                alpha=alpha,
                meshSize=mesh_size,
                pairGridSize=pg_size,
                splineOrder=spline_order,
                tolerance=1e-5
            )
            
            # Calculate energy
            elec, vdw, energy_dict = pygcmc.computeSystemEnergyPGP(self.system.state)
            
            # Store results
            results.append({
                'pair_grid_size': pg_size,
                'elec': elec,
                'vdw': vdw,
                'total': energy_dict['total'],
                'real_space': energy_dict['real_space'],
                'reciprocal': energy_dict['reciprocal'],
                'self': energy_dict['self']
            })
            
            # Results should be finite
            self.assertTrue(np.isfinite(elec))
            self.assertTrue(np.isfinite(vdw))
            self.assertTrue(np.isfinite(energy_dict['total']))
            
        # Print results for comparison
        for res in results:
            print(f"Pair Grid {res['pair_grid_size']}: elec={res['elec']:.6f}, total={res['total']:.6f}")
            
        # Results should be similar until pair grid component is fully implemented
        # In future, we'll test the convergence with increasing grid size

if __name__ == '__main__':
    unittest.main() 