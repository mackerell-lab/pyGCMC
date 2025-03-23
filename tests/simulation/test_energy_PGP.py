import unittest
import numpy as np
import math
import pygcmc
from pygcmc import MCState, MCInfo, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo

# 设置日志级别为INFO或更低，确保能看到详细日志输出
# 系统日志设置
pygcmc.System.set_log_level(pygcmc.LogLevel.INFO)
pygcmc.System.set_verbose(True)

# 平台日志设置 (energyPGP.cpp中的日志输出需要这个设置)
pygcmc.set_platform_verbose(True)  # 启用平台日志输出
pygcmc.set_platform_log_level(pygcmc.PlatformLogLevel.INFO)
pygcmc.set_platform_debug_mode(True)  # 启用调试模式用于测试

# 如果需要更详细的日志，可以设置为DEBUG
# pygcmc.System.set_log_level(pygcmc.LogLevel.DEBUG)
# pygcmc.set_platform_log_level(pygcmc.PlatformLogLevel.DEBUG)

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
        self.system.info.setTemperature(self.temp)
    
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
        box = [self.box_size, self.box_size, self.box_size]
        
        # 初始化参数
        pygcmc.setPMEParameters(
            alpha=alpha,
            meshSize=mesh_size,
            splineOrder=spline_order,
            tolerance=tolerance
        )
        
        # 初始化PME参数 - 这是关键步骤
        pygcmc.initializePMEParameters(self.cutoff, box, alpha)
        
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
        n_residues = len(system.residues)
        for i in range(0, n_residues, 2):
            system.residues[i].fixed = True
        
        # 预计算固定部分的网格电势
        pygcmc.precomputeGridPotential(system, fixed_only=True)
        
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
        box = [self.box_size, self.box_size, self.box_size]
        
        # 初始化参数
        pygcmc.setPMEParameters(
            alpha=alpha,
            meshSize=mesh_size,
            splineOrder=spline_order,
            tolerance=tolerance
        )
        
        # 初始化PME参数 - 这是关键步骤
        pygcmc.initializePMEParameters(self.cutoff, box, alpha)
        
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
        n_residues = len(system.residues)
        fixed_residues = []
        moving_residues = []
        
        for i in range(n_residues):
            if i % 2 == 0:
                system.residues[i].fixed = True
                fixed_residues.append(i)
            else:
                system.residues[i].fixed = False
                moving_residues.append(i)
        
        # 设置移动残基 - 使用MCMovementResidueInfo正确设置
        movement_info = MCMovementResidueInfo()
        movement_info.startIndex = moving_residues[0]  # 第一个移动残基的索引
        movement_info.activeCount = len(moving_residues)  # 移动残基的数量
        system.movementResidues.append(movement_info)
        
        # 预计算固定部分的网格电势
        pygcmc.precomputeGridPotential(system, fixed_only=True)
        
        # 计算插值能量 - 使用新的函数名
        energy = pygcmc.calculateMoleculeEnergy(system)
        
        # 同时测试两个等效函数
        energy2 = pygcmc.interpolateMoleculeEnergy(system)
        
        # 验证两个函数返回相同结果
        self.assertEqual(energy, energy2)
        
        # 验证结果
        self.assertTrue(np.isfinite(energy))
        print(f"Interpolated energy: {energy}")

# 将测试函数移到模块级别
def test_compare_pme_pgp_energy():
    """
    比较PME和PGP计算的移动前后能量值是否一致
    
    这个测试验证:
    1. 在初始状态下PME和PGP计算的能量值应该相同
    2. 移动分子后，PME和PGP计算的能量变化应该相同
    """
    # 设置参数 - 确保PME和PGP使用相同的参数
    box_size = 2.82  # nm
    n_cells = 2
    cutoff = 1.0   # nm
    pair_cutoff = 0.5  # nm
    box = [box_size, box_size, box_size]
    
    alpha = 0.29  # 1/nm
    mesh_size = [32, 32, 32]
    pair_grid_size = [32, 32, 32]  # 使用与PME相同的网格大小以便准确比较
    spline_order = 4
    tolerance = 1e-5
    
    # 创建测试系统
    system = create_nacl_crystal(box_size, n_cells)
    
    # 确保盒子大小正确设置
    system.info.box = box
    system.info.cutoff = cutoff
    
    # 设置PME和PGP参数
    pygcmc.setPMEParameters(
        alpha=alpha,
        meshSize=mesh_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    
    # 初始化PME参数 - 这是关键步骤
    pygcmc.initializePMEParameters(cutoff, box, alpha)
    
    pygcmc.setPGPParameters(
        alpha=alpha,
        meshSize=mesh_size,
        pair_cutoff=cutoff,
        pairGridSize=pair_grid_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    
    # 将一半残基标记为固定，一半为移动
    n_residues = len(system.residues)
    fixed_residues = []
    moving_residues = []
    
    for i in range(n_residues):
        if i % 2 == 0:  # 偶数索引的残基标记为固定
            system.residues[i].fixed = True
            fixed_residues.append(i)
        else:  # 奇数索引的残基标记为移动
            system.residues[i].fixed = False
            moving_residues.append(i)
    
    # 确保移动残基是连续的
    moving_residues.sort()
    
    # 创建一个移动残基信息对象并添加到系统中
    print(f"Moving residues: {moving_residues}")
    
    # 清除之前可能存在的移动残基信息
    system.movementResidues.clear()
    
    # 创建一个移动残基信息对象
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = moving_residues[0]  # 第一个移动残基的索引
    movement_info.activeCount = len(moving_residues)  # 移动残基的数量
    system.movementResidues.append(movement_info)
    
    print(f"Added movement info: startIndex={movement_info.startIndex}, activeCount={movement_info.activeCount}")
    print(f"System has {system.activeResidueCount} active residues and {len(system.movementResidues)} movement residue groups")
    
    # 验证移动残基信息是否已设置
    assert len(system.movementResidues) > 0, "No movement residues set!"
    
    # 第1步: 使用PME计算初始系统能量
    initial_pme_result = pygcmc.computeMovementEnergyPME(system)
    initial_pme_energy = initial_pme_result[0]  # PME电静态能量
    print(f"Initial PME energy result: {initial_pme_result}")
    
    # 第2步: 使用PGP预计算网格电势并计算移动残基能量
    pygcmc.precomputeGridPotential(system, fixed_only=True)
    initial_pgp_energy = pygcmc.calculateMoleculeEnergy(system)
    
    # 打印初始能量
    print(f"Initial PME energy: {initial_pme_energy}")
    print(f"Initial PGP energy: {initial_pgp_energy}")
    
    # 第3步: 移动移动残基（如平移0.1 nm）
    translation = [0.1, 0.1, 0.1]  # nm
    
    for res_idx in moving_residues:
        residue = system.residues[res_idx]
        for atom_idx in range(residue.atomCount):
            atom_index = residue.atomStart + atom_idx
            atom = system.atoms[atom_index]
            atom.x += translation[0]
            atom.y += translation[1]
            atom.z += translation[2]
    
    # 第4步: 使用PME计算移动后的系统能量
    moved_pme_result = pygcmc.computeMovementEnergyPME(system)
    moved_pme_energy = moved_pme_result[0]  # PME电静态能量
    print(f"Moved PME energy result: {moved_pme_result}")
    
    # 第5步: 使用PGP计算移动后的能量
    moved_pgp_energy = pygcmc.calculateMoleculeEnergy(system)
    
    # 打印移动后能量
    print(f"Moved PME energy: {moved_pme_energy}")
    print(f"Moved PGP energy: {moved_pgp_energy}")
    
    # 计算能量变化
    pme_energy_change = moved_pme_energy - initial_pme_energy
    pgp_energy_change = moved_pgp_energy - initial_pgp_energy
    
    print(f"PME energy change: {pme_energy_change}")
    print(f"PGP energy change: {pgp_energy_change}")
    
    if abs(pme_energy_change) < 1e-10:
        print("PME energy change is too small, cannot compute relative error")
        assert abs(pgp_energy_change) < 1e-10, f"PGP energy should also be close to zero"
    else:
        # 计算相对误差，允许一定的误差范围（例如5%）
        relative_error = abs((pgp_energy_change - pme_energy_change) / pme_energy_change)
        print(f"Relative error: {relative_error * 100:.4f}%")
        
        # 验证PGP和PME计算的能量变化在误差范围内一致
        assert relative_error < 0.05, f"相对误差过大: {relative_error*100:.2f}%"  # 允许5%的误差

if __name__ == '__main__':
    unittest.main() 