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
    
    def test_pgp_parameter_setting(self):
        """
        Test that PGP parameters can be set properly
        """
        # Just test that the parameter setting doesn't throw an exception
        alpha = 0.29  # 1/nm
        mesh_size = [32, 32, 32]
        pair_grid_size = [16, 16, 16]
        spline_order = 4
        tolerance = 1e-5
        
        # Set the parameters
        pygcmc.setPGPParameters(
            alpha=alpha,
            meshSize=mesh_size,
            pair_cutoff=self.pair_cutoff,
            pairGridSize=pair_grid_size,
            splineOrder=spline_order,
            tolerance=tolerance
        )
        
        # This test passes if setPGPParameters doesn't throw an exception
        self.assertTrue(True)
        
    def test_precompute_grid_potential(self):
        """
        Test precomputing the grid potential
        """
        # 设置PGP参数
        alpha = 0.29  # 1/nm
        mesh_size = [32, 32, 32]
        pair_grid_size = [16, 16, 16]
        spline_order = 4
        tolerance = 1e-5
        
        # 初始化参数
        pygcmc.setPGPParameters(
            alpha=alpha,
            meshSize=mesh_size,
            pair_cutoff=self.pair_cutoff,
            pairGridSize=pair_grid_size,
            splineOrder=spline_order,
            tolerance=tolerance
        )
        
        # 创建一个包含固定和移动部分的模型
        system = create_nacl_crystal(self.box_size, self.n_cells)
        
        # 将一半的残基标记为固定
        n_residues = len(system.state.residues)
        for i in range(0, n_residues, 2):
            system.state.residues[i].fixed = True
        
        # 预计算固定部分的网格电势
        pygcmc.precomputeGridPotential(system.state, fixed_only=True)
        
        # 测试成功执行而不崩溃
        self.assertTrue(True)
        print("Grid potential precomputation succeeded")
        
    def test_interpolate_molecule_energy(self):
        """
        Test interpolating molecule energy from the precomputed grid
        """
        # 设置PGP参数
        alpha = 0.29  # 1/nm
        mesh_size = [32, 32, 32]
        pair_grid_size = [16, 16, 16]
        spline_order = 4
        tolerance = 1e-5
        
        # 初始化参数
        pygcmc.setPGPParameters(
            alpha=alpha,
            meshSize=mesh_size,
            pair_cutoff=self.pair_cutoff,
            pairGridSize=pair_grid_size,
            splineOrder=spline_order,
            tolerance=tolerance
        )
        
        # 创建一个包含固定和移动部分的模型
        system = create_nacl_crystal(self.box_size, self.n_cells)
        
        # 将一半的残基标记为固定，一半为移动
        n_residues = len(system.state.residues)
        for i in range(0, n_residues, 2):
            system.state.residues[i].fixed = True
        
        # 设置移动残基
        movement_residues = []
        for i in range(1, n_residues, 2):  # 选择非固定残基
            movement_residues.append(i)
        system.set_movement_residues(movement_residues)
        
        # 预计算固定部分的网格电势
        pygcmc.precomputeGridPotential(system.state, fixed_only=True)
        
        # 计算插值能量 - 使用新的函数名
        energy = pygcmc.calculateMoleculeEnergy(system.state)
        
        # 同时测试两个等效函数
        energy2 = pygcmc.interpolateMoleculeEnergy(system.state)
        
        # 验证两个函数返回相同结果
        self.assertEqual(energy, energy2)
        
        # 验证结果
        self.assertTrue(np.isfinite(energy))
        print(f"Interpolated energy: {energy}")

if __name__ == '__main__':
    unittest.main() 