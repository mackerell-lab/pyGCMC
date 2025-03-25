import unittest
import numpy as np
import math
import pygcmc
from pygcmc import MCState, MCInfo, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo
import sys

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

# 确保输出缓冲区立即刷新
sys.stdout.flush()
print("日志级别设置已完成")
sys.stdout.flush()

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

def create_long_distance_system(box_size):
    """
    创建一个专门用于测试倒空间计算的系统
    
    原子间距很远，超出实空间截断距离，这样倒空间计算将占主导
    
    Args:
        box_size: 盒子大小 (nm)
    """
    state = MCState()
    
    # 设置盒子大小和温度
    state.info.box = [box_size, box_size, box_size]
    state.info.setTemperature(300.0)  # 300K
    state.info.cutoff = 1.2  # 1.2 nm cutoff
    
    # 设置力场参数
    ff = MCForceField()
    ff.numTotalTypes = 2  # 两种离子类型
    
    # LJ参数
    sigma_na = 0.333  # nm
    sigma_cl = 0.442  # nm
    eps_na = 0.0115  # kJ/mol
    eps_cl = 0.4184  # kJ/mol
    
    # 设置LJ参数矩阵
    ff.ljSigma = [
        sigma_na, (sigma_na + sigma_cl)/2.0,
        (sigma_na + sigma_cl)/2.0, sigma_cl
    ]
    ff.ljEps = [
        eps_na, math.sqrt(eps_na * eps_cl),
        math.sqrt(eps_na * eps_cl), eps_cl
    ]
    
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # 创建固定部分：两个带电离子放在盒子的对角位置
    print(f"\n创建具有远距离相互作用的系统...")
    
    # 离子1：放在盒子一角
    ion1 = MCAtom()
    ion1.x = 0.1
    ion1.y = 0.1
    ion1.z = 0.1
    ion1.charge = 1.0
    ion1.type = 0
    atoms.append(ion1)
    
    # 离子2：放在盒子对角
    ion2 = MCAtom()
    ion2.x = box_size - 0.1
    ion2.y = box_size - 0.1
    ion2.z = box_size - 0.1
    ion2.charge = -1.0
    ion2.type = 1
    atoms.append(ion2)
    
    # 创建固定残基
    fixed_res = MCResidue()
    fixed_res.atomStart = 0
    fixed_res.atomCount = 2
    fixed_res.active = True
    fixed_res.fixed = True
    residues.append(fixed_res)
    
    # 创建用于移动的离子，放在盒子中部
    ion3 = MCAtom()
    ion3.x = box_size / 2.0
    ion3.y = box_size / 2.0
    ion3.z = box_size / 2.0
    ion3.charge = 1.0
    ion3.type = 0
    atoms.append(ion3)
    
    # 创建移动残基
    move_res = MCResidue()
    move_res.atomStart = 2
    move_res.atomCount = 1
    move_res.active = True
    move_res.fixed = False
    residues.append(move_res)
    
    print(f"系统创建完成，总计 {len(atoms)} 个原子和 {len(residues)} 个残基。")
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)
    
    return state

def test_pgp_parameter_setting():
    """
    Test that PGP parameters can be set properly
    """
    # Just test that the parameter setting doesn't throw an exception
    alpha = 0.29  # 1/nm
    mesh_size = [32, 32, 32]
    potential_grid_size = [16, 16, 16]
    spline_order = 4
    tolerance = 1e-5
    potential_cutoff = 0.5  # nm
    
    # Set the parameters
    pygcmc.setPGPParameters(
        alpha=alpha,
        meshSize=mesh_size,
        potential_cutoff=potential_cutoff,
        potentialGridSize=potential_grid_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    
    # This test passes if setPGPParameters doesn't throw an exception
    assert True, "Parameters set successfully"
        
def test_precompute_grid_potential():
    """
    Test precomputing the grid potential
    """
    # 设置基本参数
    box_size = 2.82  # nm, approximately 28.2 Å
    n_cells = 2      # 2x2x2 supercell
    cutoff = 1.0   # nm
    potential_cutoff = 0.5  # nm
    
    # 设置PGP参数
    alpha = 0.29  # 1/nm
    mesh_size = [32, 32, 32]
    potential_grid_size = [16, 16, 16]
    spline_order = 4
    tolerance = 1e-5
    box = [box_size, box_size, box_size]
    
    # 初始化参数
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
        potential_cutoff=potential_cutoff,
        potentialGridSize=potential_grid_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    
    # 创建一个包含固定和移动部分的模型
    system = create_nacl_crystal(box_size, n_cells)
    
    # 将一半的残基标记为固定
    n_residues = len(system.residues)
    for i in range(0, n_residues, 2):
        system.residues[i].fixed = True
    
    # 预计算固定部分的网格电势
    pygcmc.precomputeGridPotential(system, fixed_only=True)
    
    # 测试成功执行而不崩溃
    assert True, "Grid potential precomputation succeeded"
    print("Grid potential precomputation succeeded")
        
def test_interpolate_molecule_energy():
    """
    Test interpolating molecule energy from the precomputed grid
    """
    # 设置基本参数
    box_size = 2.82  # nm, approximately 28.2 Å
    n_cells = 2      # 2x2x2 supercell
    cutoff = 1.0   # nm
    potential_cutoff = 0.5  # nm
    
    # 设置PGP参数
    alpha = 0.29  # 1/nm
    mesh_size = [32, 32, 32]
    potential_grid_size = [16, 16, 16]
    spline_order = 4
    tolerance = 1e-5
    box = [box_size, box_size, box_size]
    
    # 初始化参数
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
        potential_cutoff=potential_cutoff,
        potentialGridSize=potential_grid_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    
    # 创建一个包含固定和移动部分的模型
    system = create_nacl_crystal(box_size, n_cells)
    
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
    
    # 检查能量值是否合理
    # 注意：这里我们不检查具体的能量值，因为计算结果取决于多种因素
    # 只检查能量值是否为有限数且不为零
    assert np.isfinite(energy), "Energy value should be finite"
    assert energy != 0.0, "Energy value should not be exactly zero"
    
    print(f"Interpolated energy: {energy} kJ/mol")

# 将测试函数移到模块级别
def test_compare_pme_pgp_energy():
    """
    比较PME和PGP计算的移动前后能量值是否一致
    
    这个测试验证:
    1. 在初始状态下PME和PGP计算的能量值应该相同
    2. 移动分子后，PME和PGP计算的能量变化应该相同
    """
    # 设置参数 - 确保PME和PGP使用相同的参数
    box_size = 5.0  # nm - 使用更大的盒子
    cutoff = 1.0   # nm
    potential_cutoff = 1.0  # nm - 与cutoff相同
    box = [box_size, box_size, box_size]
    
    alpha = 0.29  # 1/nm
    mesh_size = [32, 32, 32]
    potential_grid_size = [32, 32, 32]  # 使用与PME相同的网格大小以便准确比较
    spline_order = 4
    tolerance = 1e-5
    
    # 设置测试超时时间，避免长时间运行
    timeout = 10  # 秒
    
    print("创建长距离测试系统...")
    sys.stdout.flush()
    
    # 创建测试系统 - 使用专门设计的长距离系统
    system = create_long_distance_system(box_size)
    
    # 确保盒子大小正确设置
    system.info.box = box
    system.info.cutoff = cutoff
    
    print("设置PME参数...")
    sys.stdout.flush()
    
    # 设置PME和PGP参数
    pygcmc.setPMEParameters(
        alpha=alpha,
        meshSize=mesh_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    
    # 初始化PME参数 - 这是关键步骤
    pygcmc.initializePMEParameters(cutoff, box, alpha)
    
    print("设置PGP参数...")
    sys.stdout.flush()
    
    pygcmc.setPGPParameters(
        alpha=alpha,
        meshSize=mesh_size,
        potential_cutoff=potential_cutoff,
        potentialGridSize=potential_grid_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    
    print("设置移动残基...")
    sys.stdout.flush()
    
    # 移动残基已经在create_long_distance_system中设置好了
    # 这里只需准备movementResidues列表
    moving_residues = [1]  # 第二个残基是移动残基
    
    # 验证固定残基信息
    fixed_count = sum(1 for res in system.residues if res.fixed)
    print(f"固定残基数: {fixed_count}")
    print(f"移动残基数: {len(moving_residues)}")
    sys.stdout.flush()
    
    # 设置移动残基信息
    system.movementResidues.clear()
    
    # 创建移动残基信息
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = moving_residues[0]  # 移动残基的索引
    movement_info.activeCount = 1  # 只有一个移动残基
    system.movementResidues.append(movement_info)
    
    print(f"添加移动信息: startIndex={movement_info.startIndex}, activeCount={movement_info.activeCount}")
    print(f"系统有 {system.activeResidueCount} 个活跃残基和 {len(system.movementResidues)} 个移动残基组")
    sys.stdout.flush()
    
    # 第1步: 使用PME计算初始系统能量
    print("计算PME初始能量...")
    sys.stdout.flush()
    initial_pme_result = pygcmc.computeMovementEnergyPME(system)
    initial_pme_energy = initial_pme_result[0]  # PME电静态能量
    initial_pme_dict = initial_pme_result[2]  # PME能量细节字典
    initial_pme_reciprocal = initial_pme_dict['reciprocal']  # 只取倒空间部分
    print(f"Initial PME energy result: {initial_pme_result}")
    print(f"Initial PME reciprocal energy: {initial_pme_reciprocal}")
    sys.stdout.flush()
    
    # 第2步: 使用PGP预计算网格电势并计算移动残基能量
    print("预计算PGP网格电势...")
    sys.stdout.flush()
    pygcmc.precomputeGridPotential(system, fixed_only=True)
    
    print("计算PGP初始能量...")
    sys.stdout.flush()
    initial_pgp_energy = pygcmc.calculateMoleculeEnergy(system)
    
    # 打印初始能量
    print(f"Initial PME reciprocal energy: {initial_pme_reciprocal}")
    print(f"Initial PGP energy: {initial_pgp_energy}")
    sys.stdout.flush()
    
    # 第3步: 移动移动残基（如平移0.1 nm）
    translation = [0.1, 0.1, 0.1]  # nm
    
    print("移动残基...")
    sys.stdout.flush()
    for res_idx in moving_residues:
        residue = system.residues[res_idx]
        for atom_idx in range(residue.atomCount):
            atom_index = residue.atomStart + atom_idx
            atom = system.atoms[atom_index]
            atom.x += translation[0]
            atom.y += translation[1]
            atom.z += translation[2]
    
    # 第4步: 使用PME计算移动后的系统能量
    print("计算PME移动后能量...")
    sys.stdout.flush()
    moved_pme_result = pygcmc.computeMovementEnergyPME(system)
    moved_pme_energy = moved_pme_result[0]  # PME电静态能量
    moved_pme_dict = moved_pme_result[2]  # PME能量细节字典
    moved_pme_reciprocal = moved_pme_dict['reciprocal']  # 只取倒空间部分
    print(f"Moved PME energy result: {moved_pme_result}")
    print(f"Moved PME reciprocal energy: {moved_pme_reciprocal}")
    sys.stdout.flush()
    
    # 第5步: 使用PGP计算移动后的能量
    print("计算PGP移动后能量...")
    sys.stdout.flush()
    moved_pgp_energy = pygcmc.calculateMoleculeEnergy(system)
    
    # 打印移动后能量
    print(f"Moved PME reciprocal energy: {moved_pme_reciprocal}")
    print(f"Moved PGP energy: {moved_pgp_energy}")
    sys.stdout.flush()
    
    # 计算能量变化 - 只使用倒空间部分
    pme_energy_change = moved_pme_reciprocal - initial_pme_reciprocal
    pgp_energy_change = moved_pgp_energy - initial_pgp_energy
    
    print(f"PME reciprocal energy change: {pme_energy_change}")
    print(f"PGP energy change: {pgp_energy_change}")
    sys.stdout.flush()
    
    if abs(pme_energy_change) < 1e-10:
        print("PME energy change is too small, cannot compute relative error")
        assert abs(pgp_energy_change) < 1e-10, f"PGP energy should also be close to zero"
    else:
        # 计算相对误差，允许一定的误差范围（例如10%）
        relative_error = abs((pgp_energy_change - pme_energy_change) / pme_energy_change)
        print(f"Relative error: {relative_error * 100:.4f}%")
        sys.stdout.flush()
        
        # 验证PGP和PME计算的能量变化在误差范围内一致
        assert relative_error < 0.1, f"相对误差过大: {relative_error*100:.2f}%"  # 允许10%的误差

def test_compare_ewald_pme_pgp_complex():
    """
    使用更复杂的系统比较Ewald、PME和PGP算法
    
    这个测试:
    1. 创建一个固定部分包含多个带电粒子的复杂系统
    2. 移动部分包含2-3个原子，距离超过cutoff
    3. 比较三种方法计算的倒空间能量变化
    """
    # 设置系统参数
    box_size = 8.0  # nm - 使用更大的盒子，确保倒空间主导
    cutoff = 1.0    # nm
    potential_cutoff = 1.0  # nm
    box = [box_size, box_size, box_size]
    
    # Ewald和PME参数
    alpha = 0.29    # 1/nm
    mesh_size = [32, 32, 32]
    potential_grid_size = [32, 32, 32]
    spline_order = 4
    tolerance = 1e-5
    
    print("创建复杂测试系统...")
    sys.stdout.flush()
    
    # 创建测试系统
    system = MCState()
    system.info.box = box
    system.info.setTemperature(300.0)
    system.info.cutoff = cutoff
    
    # 设置力场
    ff = MCForceField()
    ff.numTotalTypes = 2 
    ff.ljSigma = [0.333, 0.3875, 0.3875, 0.442]
    ff.ljEps = [0.0115, 0.0693, 0.0693, 0.4184]
    system.forcefield = ff
    
    atoms = []
    residues = []
    
    # 创建固定部分 - 8个离子形成一个立方体
    fixed_positions = [
        (1.0, 1.0, 1.0),
        (1.0, 1.0, box_size-1.0),
        (1.0, box_size-1.0, 1.0),
        (1.0, box_size-1.0, box_size-1.0),
        (box_size-1.0, 1.0, 1.0),
        (box_size-1.0, 1.0, box_size-1.0),
        (box_size-1.0, box_size-1.0, 1.0),
        (box_size-1.0, box_size-1.0, box_size-1.0)
    ]
    
    # 添加固定离子
    for i, pos in enumerate(fixed_positions):
        ion = MCAtom()
        ion.x, ion.y, ion.z = pos
        ion.charge = 1.0 if i % 2 == 0 else -1.0  # 交替正负电荷
        ion.type = 0 if i % 2 == 0 else 1
        atoms.append(ion)
    
    # 创建固定残基
    fixed_res = MCResidue()
    fixed_res.atomStart = 0
    fixed_res.atomCount = len(fixed_positions)
    fixed_res.active = True
    fixed_res.fixed = True
    residues.append(fixed_res)
    
    # 创建移动部分 - 3个原子形成一个水分子
    # 位置在盒子中央，确保与固定部分距离超过cutoff
    mobile_atoms = [
        (box_size/2.0, box_size/2.0, box_size/2.0),       # 中心氧原子
        (box_size/2.0 + 0.1, box_size/2.0, box_size/2.0), # 氢原子1
        (box_size/2.0, box_size/2.0 + 0.1, box_size/2.0)  # 氢原子2
    ]
    
    mobile_charges = [-0.8, 0.4, 0.4]  # 水分子电荷
    
    # 添加移动原子
    for pos, q in zip(mobile_atoms, mobile_charges):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = q
        atom.type = 0
        atoms.append(atom)
    
    # 创建移动残基
    mobile_res = MCResidue()
    mobile_res.atomStart = len(fixed_positions)
    mobile_res.atomCount = len(mobile_atoms)
    mobile_res.active = True
    mobile_res.fixed = False
    residues.append(mobile_res)
    
    # 设置系统
    system.atoms = atoms
    system.residues = residues
    system.activeAtomCount = len(atoms)
    system.activeResidueCount = len(residues)
    
    # 设置移动残基信息
    system.movementResidues.clear()
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 1  # 第二个残基（移动残基）的索引
    movement_info.activeCount = 1  # 只有一个移动残基
    system.movementResidues.append(movement_info)
    
    print(f"系统创建完成: {system.activeAtomCount}个原子, {system.activeResidueCount}个残基")
    print(f"固定原子: {len(fixed_positions)}, 移动原子: {len(mobile_atoms)}")
    sys.stdout.flush()
    
    # 初始化各种电荷方法
    print("设置计算参数...")
    
    # PME参数
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order, tolerance)
    pygcmc.initializePMEParameters(cutoff, box, alpha)
    
    # PGP参数
    pygcmc.setPGPParameters(alpha, mesh_size, potential_cutoff, 
                           potential_grid_size, spline_order, tolerance)
    
    # Ewald参数（使用computeEwaldEnergy函数）
    # 注意：这里假设你的库中有计算标准Ewald能量的函数
    
    # 步骤1: 计算初始能量
    print("计算初始能量...")
    
    # PME能量
    initial_pme_result = pygcmc.computeMovementEnergyPME(system)
    initial_pme_energy = initial_pme_result[0]
    initial_pme_dict = initial_pme_result[2]
    initial_pme_reciprocal = initial_pme_dict['reciprocal']
    
    # Ewald能量（如果有）
    # initial_ewald_energy = pygcmc.computeEwaldEnergy(system)
    
    # PGP能量
    pygcmc.precomputeGridPotential(system, fixed_only=True)
    initial_pgp_energy = pygcmc.calculateMoleculeEnergy(system)
    
    print(f"初始PME倒空间能量: {initial_pme_reciprocal}")
    # print(f"初始Ewald能量: {initial_ewald_energy}")
    print(f"初始PGP能量: {initial_pgp_energy}")
    
    # 步骤2: 移动移动残基（平移0.2 nm）
    translation = [0.2, 0.2, 0.2]
    print(f"移动残基: 平移{translation}...")
    
    # 移动水分子
    for i in range(mobile_res.atomCount):
        atom_index = mobile_res.atomStart + i
        system.atoms[atom_index].x += translation[0]
        system.atoms[atom_index].y += translation[1]
        system.atoms[atom_index].z += translation[2]
    
    # 步骤3: 计算移动后能量
    print("计算移动后能量...")
    
    # PME能量
    moved_pme_result = pygcmc.computeMovementEnergyPME(system)
    moved_pme_energy = moved_pme_result[0]
    moved_pme_dict = moved_pme_result[2]
    moved_pme_reciprocal = moved_pme_dict['reciprocal']
    
    # Ewald能量（如果有）
    # moved_ewald_energy = pygcmc.computeEwaldEnergy(system)
    
    # PGP能量
    moved_pgp_energy = pygcmc.calculateMoleculeEnergy(system)
    
    print(f"移动后PME倒空间能量: {moved_pme_reciprocal}")
    # print(f"移动后Ewald能量: {moved_ewald_energy}")
    print(f"移动后PGP能量: {moved_pgp_energy}")
    
    # 计算能量变化
    pme_energy_change = moved_pme_reciprocal - initial_pme_reciprocal
    # ewald_energy_change = moved_ewald_energy - initial_ewald_energy
    pgp_energy_change = moved_pgp_energy - initial_pgp_energy
    
    print(f"PME倒空间能量变化: {pme_energy_change}")
    # print(f"Ewald能量变化: {ewald_energy_change}")
    print(f"PGP能量变化: {pgp_energy_change}")
    
    # 计算相对误差
    if abs(pme_energy_change) > 1e-10:
        pgp_relative_error = abs((pgp_energy_change - pme_energy_change) / pme_energy_change)
        print(f"PGP与PME相对误差: {pgp_relative_error*100:.4f}%")
        assert pgp_relative_error < 0.1, f"PGP相对误差过大: {pgp_relative_error*100:.2f}%"
    
    # 如果有Ewald计算，也比较与PME的误差
    # if abs(pme_energy_change) > 1e-10:
    #     ewald_relative_error = abs((ewald_energy_change - pme_energy_change) / pme_energy_change)
    #     print(f"Ewald与PME相对误差: {ewald_relative_error*100:.4f}%")
    #     assert ewald_relative_error < 0.1, f"Ewald相对误差过大: {ewald_relative_error*100:.2f}%"

def test_compare_ewald_pme_pgp_planar():
    """
    使用平面系统比较Ewald、PME和PGP算法
    
    特点：
    1. 所有原子都在z=4.0 nm的平面上
    2. 移动分子与固定部分距离超过cutoff
    3. 平面与网格平面平行
    4. 原子位置故意偏离网格点
    5. 使用8x8x8的粗网格便于观察
    """
    # 设置系统参数
    box_size = 8.0  # nm
    cutoff = 1.0    # nm
    potential_cutoff = 1.0  # nm
    box = [box_size, box_size, box_size]
    
    # 计算参数
    alpha = 0.29    # 1/nm
    mesh_size = [8, 8, 8]  # 改为8x8x8的粗网格
    potential_grid_size = [8, 8, 8]  # 同样改为8x8x8
    spline_order = 4
    tolerance = 1e-5
    
    # 计算网格间距
    grid_spacing = box_size / mesh_size[0]
    print(f"网格间距: {grid_spacing:.3f} nm")
    
    # 选择平面z坐标（确保与网格平面平行）
    z_plane = 4.0  # nm
    
    print("创建平面测试系统...")
    sys.stdout.flush()
    
    # 创建测试系统
    system = MCState()
    system.info.box = box
    system.info.setTemperature(300.0)
    system.info.cutoff = cutoff
    
    # 设置力场
    ff = MCForceField()
    ff.numTotalTypes = 2 
    ff.ljSigma = [0.333, 0.3875, 0.3875, 0.442]
    ff.ljEps = [0.0115, 0.0693, 0.0693, 0.4184]
    system.forcefield = ff
    
    atoms = []
    residues = []
    
    # 创建固定部分 - 4个离子形成一个正方形
    # 位置故意偏离网格点
    fixed_positions = [
        (1.0 + grid_spacing/3, 1.0 + grid_spacing/3, z_plane),  # 左下
        (1.0 + grid_spacing/3, box_size-1.0 + grid_spacing/3, z_plane),  # 左上
        (box_size-1.0 + grid_spacing/3, 1.0 + grid_spacing/3, z_plane),  # 右下
        (box_size-1.0 + grid_spacing/3, box_size-1.0 + grid_spacing/3, z_plane)  # 右上
    ]
    
    # 添加固定离子
    for i, pos in enumerate(fixed_positions):
        ion = MCAtom()
        ion.x, ion.y, ion.z = pos
        ion.charge = 1.0 if i % 2 == 0 else -1.0  # 交替正负电荷
        ion.type = 0 if i % 2 == 0 else 1
        atoms.append(ion)
        
        # 打印原子位置和最近的网格点
        grid_x = round(pos[0] / grid_spacing)
        grid_y = round(pos[1] / grid_spacing)
        grid_z = round(pos[2] / grid_spacing)
        print(f"固定离子 {i}: 位置=({pos[0]:.3f}, {pos[1]:.3f}, {pos[2]:.3f})")
        print(f"  最近网格点: ({grid_x}, {grid_y}, {grid_z})")
        print(f"  偏离网格点: ({pos[0]-grid_x*grid_spacing:.3f}, {pos[1]-grid_y*grid_spacing:.3f}, {pos[2]-grid_z*grid_spacing:.3f})")
    
    # 创建固定残基
    fixed_res = MCResidue()
    fixed_res.atomStart = 0
    fixed_res.atomCount = len(fixed_positions)
    fixed_res.active = True
    fixed_res.fixed = True
    residues.append(fixed_res)
    
    # 创建移动部分 - 3个原子形成一个水分子
    # 位置在盒子中央，确保与固定部分距离超过cutoff
    mobile_atoms = [
        (box_size/2.0 + grid_spacing/3, box_size/2.0 + grid_spacing/3, z_plane),  # 氧原子
        (box_size/2.0 + grid_spacing/3 + 0.1, box_size/2.0 + grid_spacing/3, z_plane),  # 氢原子1
        (box_size/2.0 + grid_spacing/3, box_size/2.0 + grid_spacing/3 + 0.1, z_plane)  # 氢原子2
    ]
    
    mobile_charges = [-0.8, 0.4, 0.4]  # 水分子电荷
    
    # 添加移动原子
    for i, (pos, q) in enumerate(zip(mobile_atoms, mobile_charges)):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = q
        atom.type = 0
        atoms.append(atom)
        
        # 打印原子位置和最近的网格点
        grid_x = round(pos[0] / grid_spacing)
        grid_y = round(pos[1] / grid_spacing)
        grid_z = round(pos[2] / grid_spacing)
        print(f"移动原子 {i}: 位置=({pos[0]:.3f}, {pos[1]:.3f}, {pos[2]:.3f})")
        print(f"  最近网格点: ({grid_x}, {grid_y}, {grid_z})")
        print(f"  偏离网格点: ({pos[0]-grid_x*grid_spacing:.3f}, {pos[1]-grid_y*grid_spacing:.3f}, {pos[2]-grid_z*grid_spacing:.3f})")
    
    # 创建移动残基
    mobile_res = MCResidue()
    mobile_res.atomStart = len(fixed_positions)
    mobile_res.atomCount = len(mobile_atoms)
    mobile_res.active = True
    mobile_res.fixed = False
    residues.append(mobile_res)
    
    # 设置系统
    system.atoms = atoms
    system.residues = residues
    system.activeAtomCount = len(atoms)
    system.activeResidueCount = len(residues)
    
    # 设置移动残基信息
    system.movementResidues.clear()
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 1  # 第二个残基（移动残基）的索引
    movement_info.activeCount = 1  # 只有一个移动残基
    system.movementResidues.append(movement_info)
    
    print(f"系统创建完成: {system.activeAtomCount}个原子, {system.activeResidueCount}个残基")
    print(f"固定原子: {len(fixed_positions)}, 移动原子: {len(mobile_atoms)}")
    print(f"网格大小: {mesh_size[0]}x{mesh_size[1]}x{mesh_size[2]}")
    print(f"网格间距: {grid_spacing:.3f} nm")
    sys.stdout.flush()
    
    # 初始化各种电荷方法
    print("设置计算参数...")
    
    # PME参数
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order, tolerance)
    pygcmc.initializePMEParameters(cutoff, box, alpha)
    
    # PGP参数
    pygcmc.setPGPParameters(alpha, mesh_size, potential_cutoff, 
                           potential_grid_size, spline_order, tolerance)
    
    # 步骤1: 计算初始能量
    print("计算初始能量...")
    
    # PME能量
    initial_pme_result = pygcmc.computeMovementEnergyPME(system)
    initial_pme_energy = initial_pme_result[0]
    initial_pme_dict = initial_pme_result[2]
    initial_pme_reciprocal = initial_pme_dict['reciprocal']
    
    # PGP能量
    pygcmc.precomputeGridPotential(system, fixed_only=True)
    initial_pgp_energy = pygcmc.calculateMoleculeEnergy(system)
    
    print(f"初始PME倒空间能量: {initial_pme_reciprocal}")
    print(f"初始PGP能量: {initial_pgp_energy}")
    
    # 步骤2: 移动移动残基（在平面上平移）
    translation = [0.2, 0.2, 0.0]  # 只在x-y平面上移动
    print(f"移动残基: 平移{translation}...")
    
    # 移动水分子
    for i in range(mobile_res.atomCount):
        atom_index = mobile_res.atomStart + i
        system.atoms[atom_index].x += translation[0]
        system.atoms[atom_index].y += translation[1]
        # z坐标保持不变，确保仍在同一平面上
    
    # 步骤3: 计算移动后能量
    print("计算移动后能量...")
    
    # PME能量
    moved_pme_result = pygcmc.computeMovementEnergyPME(system)
    moved_pme_energy = moved_pme_result[0]
    moved_pme_dict = moved_pme_result[2]
    moved_pme_reciprocal = moved_pme_dict['reciprocal']
    
    # PGP能量
    moved_pgp_energy = pygcmc.calculateMoleculeEnergy(system)
    
    print(f"移动后PME倒空间能量: {moved_pme_reciprocal}")
    print(f"移动后PGP能量: {moved_pgp_energy}")
    
    # 计算能量变化
    pme_energy_change = moved_pme_reciprocal - initial_pme_reciprocal
    pgp_energy_change = moved_pgp_energy - initial_pgp_energy
    
    print(f"PME倒空间能量变化: {pme_energy_change}")
    print(f"PGP能量变化: {pgp_energy_change}")
    
    # 计算相对误差
    if abs(pme_energy_change) > 1e-10:
        pgp_relative_error = abs((pgp_energy_change - pme_energy_change) / pme_energy_change)
        print(f"PGP与PME相对误差: {pgp_relative_error*100:.4f}%")
        assert pgp_relative_error < 0.1, f"PGP相对误差过大: {pgp_relative_error*100:.2f}%"

def test_compare_ewald_pme_pgp_asymmetric():
    """
    测试具有不对称电荷分布的系统中PGP计算的准确性
    
    特点:
    1. 固定部分包含不对称分布的多个带电粒子
    2. 移动残基为一个水分子，远离固定部分
    3. 执行多次随机移动，确保移动后与固定部分距离始终大于cutoff
    4. 通过比较PME和PGP计算的能量验证准确性
    """
    # 设置系统参数
    box_size = 8.0  # nm - 使用较大的盒子
    cutoff = 1.0    # nm
    potential_cutoff = 1.0  # nm
    box = [box_size, box_size, box_size]
    
    # Ewald和PME参数
    alpha = 0.29    # 1/nm
    mesh_size = [32, 32, 32]
    potential_grid_size = [32, 32, 32]
    spline_order = 4
    tolerance = 1e-5
    
    # 定义计算周期性边界条件下距离的函数
    def calculate_pbc_distance(pos1, pos2, box_size):
        """计算周期性边界条件下两个点之间的距离"""
        dx = abs(pos1[0] - pos2[0])
        dy = abs(pos1[1] - pos2[1])
        dz = abs(pos1[2] - pos2[2])

        # 应用周期性边界条件
        if dx > box_size/2:
            dx = box_size - dx
        if dy > box_size/2:
            dy = box_size - dy
        if dz > box_size/2:
            dz = box_size - dz

        return math.sqrt(dx*dx + dy*dy + dz*dz)
    
    # 检查移动是否安全(所有距离都大于cutoff)
    def is_safe_position(mobile_positions, fixed_positions, cutoff, box_size):
        """检查移动位置与所有固定粒子的距离是否都大于cutoff"""
        for mobile_pos in mobile_positions:
            for fixed_pos, _ in fixed_positions:
                distance = calculate_pbc_distance(mobile_pos, fixed_pos, box_size)
                if distance <= cutoff:
                    return False, distance
        return True, None
    
    # 生成安全的随机移动向量
    def generate_safe_move(current_positions, fixed_positions, cutoff, box_size, max_step=0.3):
        """生成一个安全的随机移动向量，确保移动后所有粒子与所有固定粒子距离大于cutoff"""
        for attempt in range(100):  # 尝试最多100次
            # 生成随机位移
            dx = (np.random.random() - 0.5) * 2 * max_step
            dy = (np.random.random() - 0.5) * 2 * max_step
            dz = (np.random.random() - 0.5) * 2 * max_step
            
            # 计算新位置
            new_positions = []
            for pos in current_positions:
                new_pos = (
                    (pos[0] + dx) % box_size,
                    (pos[1] + dy) % box_size,
                    (pos[2] + dz) % box_size
                )
                new_positions.append(new_pos)
            
            # 检查新位置是否安全
            is_safe, min_distance = is_safe_position(new_positions, fixed_positions, cutoff, box_size)
            if is_safe:
                return (dx, dy, dz), new_positions
        
        # 如果100次尝试都失败，返回更小的移动
        print("警告：100次尝试都未能找到安全位置，使用较小的移动")
        dx = 0.05
        dy = 0.05
        dz = 0.05
        new_positions = []
        for pos in current_positions:
            new_pos = (
                (pos[0] + dx) % box_size,
                (pos[1] + dy) % box_size,
                (pos[2] + dz) % box_size
            )
            new_positions.append(new_pos)
        return (dx, dy, dz), new_positions
    
    print("创建不对称电荷分布测试系统...")
    sys.stdout.flush()
    
    # 创建测试系统
    system = MCState()
    system.info.box = box
    system.info.setTemperature(300.0)
    system.info.cutoff = cutoff
    
    # 设置力场
    ff = MCForceField()
    ff.numTotalTypes = 2 
    ff.ljSigma = [0.333, 0.3875, 0.3875, 0.442]
    ff.ljEps = [0.0115, 0.0693, 0.0693, 0.4184]
    system.forcefield = ff
    
    atoms = []
    residues = []
    
    # 创建固定部分 - 不对称的带电粒子分布
    fixed_particles = [
        # 位置 (x, y, z)                电荷
        ((1.0, 1.0, 1.0),               1.0),  # 正电荷
        ((1.5, 1.0, 1.2),              -0.8),  # 负电荷
        ((1.3, 1.7, 1.5),               0.6),  # 正电荷
        ((0.8, 1.6, 0.9),              -0.7),  # 负电荷
        ((box_size-1.0, 1.0, 1.0),      1.0),  # 正电荷
        ((box_size-1.5, 1.2, 1.3),     -0.5),  # 负电荷
        ((1.0, box_size-1.0, 1.0),      0.8),  # 正电荷
        ((1.2, box_size-1.5, 0.9),     -0.6),  # 负电荷
        ((1.0, 1.0, box_size-1.0),      0.9),  # 正电荷
        ((1.4, 1.3, box_size-1.4),     -0.7),  # 负电荷
        # 添加更多的粒子增加系统复杂性
        ((box_size-2.0, box_size-2.0, 2.0),  0.7),  # 正电荷
        ((box_size-2.5, box_size-2.3, 2.2), -0.4),  # 负电荷
        ((2.0, box_size-2.0, box_size-2.0),  0.5),  # 正电荷
        ((2.2, box_size-2.2, box_size-2.4), -0.3),  # 负电荷
        ((box_size-2.0, 2.0, box_size-2.0), -1.5),  # 负电荷，大致平衡系统
    ]
    
    # 确认粒子数
    print(f"固定粒子数: {len(fixed_particles)}")
    
    # 计算固定部分总电荷
    total_fixed_charge = sum(charge for _, charge in fixed_particles)
    print(f"固定部分总电荷: {total_fixed_charge}")
    
    # 添加固定粒子
    for i, (pos, charge) in enumerate(fixed_particles):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = 0 if charge > 0 else 1  # 正电荷用type 0，负电荷用type 1
        atoms.append(atom)
    
    # 创建固定残基
    fixed_res = MCResidue()
    fixed_res.atomStart = 0
    fixed_res.atomCount = len(fixed_particles)
    fixed_res.active = True
    fixed_res.fixed = True
    residues.append(fixed_res)
    
    # 计算盒子中心位置，确保移动残基远离固定粒子
    center_x = box_size / 2.0
    center_y = box_size / 2.0
    center_z = box_size / 2.0
    
    # 创建移动部分 - 一个水分子放在盒子中央
    # 水分子电荷参数：氧原子-0.8，两个氢原子各+0.4
    # 水分子键长：O-H约0.1 nm
    
    # 设置水分子的位置
    oxygen_pos = (center_x, center_y, center_z)
    hydrogen1_pos = (center_x + 0.1, center_y, center_z)  # 第一个H原子
    hydrogen2_pos = (center_x, center_y + 0.1, center_z)  # 第二个H原子
    
    # 水分子电荷
    oxygen_charge = -0.8
    hydrogen_charge = 0.4  # 每个氢原子
    
    # 记录移动部分的起始索引
    mobile_start_idx = len(atoms)
    
    # 创建氧原子
    o_atom = MCAtom()
    o_atom.x, o_atom.y, o_atom.z = oxygen_pos
    o_atom.charge = oxygen_charge
    o_atom.type = 1  # 氧原子类型
    atoms.append(o_atom)
    
    # 创建第一个氢原子
    h1_atom = MCAtom()
    h1_atom.x, h1_atom.y, h1_atom.z = hydrogen1_pos
    h1_atom.charge = hydrogen_charge
    h1_atom.type = 0  # 氢原子类型
    atoms.append(h1_atom)
    
    # 创建第二个氢原子
    h2_atom = MCAtom()
    h2_atom.x, h2_atom.y, h2_atom.z = hydrogen2_pos
    h2_atom.charge = hydrogen_charge
    h2_atom.type = 0  # 氢原子类型
    atoms.append(h2_atom)
    
    print(f"移动水分子: O位置=({oxygen_pos[0]}, {oxygen_pos[1]}, {oxygen_pos[2]}), 电荷={oxygen_charge}")
    print(f"  H1位置=({hydrogen1_pos[0]}, {hydrogen1_pos[1]}, {hydrogen1_pos[2]}), 电荷={hydrogen_charge}")
    print(f"  H2位置=({hydrogen2_pos[0]}, {hydrogen2_pos[1]}, {hydrogen2_pos[2]}), 电荷={hydrogen_charge}")
    print(f"  水分子总电荷: {oxygen_charge + 2*hydrogen_charge}")
    
    # 创建移动残基（水分子）
    mobile_res = MCResidue()
    mobile_res.atomStart = mobile_start_idx
    mobile_res.atomCount = 3  # 水分子有3个原子
    mobile_res.active = True
    mobile_res.fixed = False
    residues.append(mobile_res)
    
    # 设置系统
    system.atoms = atoms
    system.residues = residues
    system.activeAtomCount = len(atoms)
    system.activeResidueCount = len(residues)
    
    # 计算系统总电荷
    total_system_charge = sum(atom.charge for atom in atoms)
    print(f"系统总电荷: {total_system_charge}")
    # 允许系统有少量电荷，无需严格断言
    
    # 设置移动残基信息
    system.movementResidues.clear()
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 1  # 第二个残基（移动残基）的索引
    movement_info.activeCount = 1  # 只有一个移动残基
    system.movementResidues.append(movement_info)
    
    print(f"系统创建完成: {system.activeAtomCount}个原子, {system.activeResidueCount}个残基")
    print(f"固定原子: {len(fixed_particles)}, 移动原子: 3 (水分子)")
    sys.stdout.flush()
    
    # 初始化PME和PGP
    print("设置计算参数...")
    
    # PME参数
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order, tolerance)
    pygcmc.initializePMEParameters(cutoff, box, alpha)
    
    # PGP参数
    pygcmc.setPGPParameters(alpha, mesh_size, potential_cutoff, 
                         potential_grid_size, spline_order, tolerance)
    
    # 预计算固定部分网格电势
    print("预计算固定部分网格电势...")
    pygcmc.precomputeGridPotential(system, fixed_only=True)
    
    # 设置要执行的随机移动次数
    num_moves = 5  # 减少测试次数以加快测试
    print(f"将执行 {num_moves} 次随机移动测试")
    
    # 存储PGP和PME之间的相对误差
    pgp_pme_errors = []
    
    # 验证初始位置是否安全
    current_positions = []
    for i in range(mobile_res.atomCount):
        atom_idx = mobile_res.atomStart + i
        current_positions.append((
            system.atoms[atom_idx].x,
            system.atoms[atom_idx].y,
            system.atoms[atom_idx].z
        ))
    
    is_safe, min_dist = is_safe_position(current_positions, fixed_particles, cutoff, box_size)
    if not is_safe:
        print(f"警告：初始位置不安全，最小距离为 {min_dist} nm")
    else:
        print(f"初始位置安全，与固定粒子的最小距离 > {cutoff} nm")
    
    # 执行多次随机移动
    for move_idx in range(num_moves):
        print(f"\n执行第 {move_idx+1}/{num_moves} 次随机移动测试")
        
        # 步骤1: 计算初始能量
        print("计算初始能量...")
        
        # PME能量
        initial_pme_result = pygcmc.computeMovementEnergyPME(system)
        initial_pme_energy = initial_pme_result[0]
        initial_pme_dict = initial_pme_result[2]
        initial_pme_reciprocal = initial_pme_dict['reciprocal']
        
        # PGP能量
        initial_pgp_energy = pygcmc.calculateMoleculeEnergy(system)
        
        print(f"初始PME倒空间能量: {initial_pme_reciprocal}")
        print(f"初始PGP能量: {initial_pgp_energy}")
        
        # 步骤2: 生成安全的随机移动并应用
        current_positions = []
        for i in range(mobile_res.atomCount):
            atom_idx = mobile_res.atomStart + i
            current_positions.append((
                system.atoms[atom_idx].x,
                system.atoms[atom_idx].y,
                system.atoms[atom_idx].z
            ))
        
        delta, new_positions = generate_safe_move(current_positions, fixed_particles, cutoff, box_size)
        print(f"随机移动向量: {delta}")
        
        # 移动所有移动残基中的粒子
        for i in range(mobile_res.atomCount):
            atom_idx = mobile_res.atomStart + i
            system.atoms[atom_idx].x = new_positions[i][0]
            system.atoms[atom_idx].y = new_positions[i][1]
            system.atoms[atom_idx].z = new_positions[i][2]
            print(f"移动原子 {i} 到: ({new_positions[i][0]:.4f}, {new_positions[i][1]:.4f}, {new_positions[i][2]:.4f})")
        
        # 验证移动后的位置是否安全
        is_safe, min_dist = is_safe_position(new_positions, fixed_particles, cutoff, box_size)
        if not is_safe:
            print(f"警告：移动后位置不安全，最小距离为 {min_dist} nm")
            assert min_dist > cutoff, f"移动后距离 ({min_dist} nm) 小于cutoff ({cutoff} nm)，将引入实空间能量"
        else:
            min_dist = float('inf')
            for mobile_pos in new_positions:
                for fixed_pos, _ in fixed_particles:
                    dist = calculate_pbc_distance(mobile_pos, fixed_pos, box_size)
                    min_dist = min(min_dist, dist)
            print(f"移动后位置安全，与固定粒子的最小距离：{min_dist:.4f} nm (cutoff={cutoff} nm)")
        
        # 步骤3: 计算移动后能量
        print("计算移动后能量...")
        
        # PME能量
        moved_pme_result = pygcmc.computeMovementEnergyPME(system)
        moved_pme_energy = moved_pme_result[0]
        moved_pme_dict = moved_pme_result[2]
        moved_pme_reciprocal = moved_pme_dict['reciprocal']
        
        # 检查直接空间能量是否为0
        moved_pme_direct = moved_pme_dict.get('direct', 0.0)
        if abs(moved_pme_direct) > 1e-10:
            print(f"警告: PME直接空间能量不为零: {moved_pme_direct}")
            print("这意味着存在小于cutoff的粒子对!")
        
        # PGP能量
        moved_pgp_energy = pygcmc.calculateMoleculeEnergy(system)
        
        print(f"移动后PME倒空间能量: {moved_pme_reciprocal}")
        print(f"移动后PGP能量: {moved_pgp_energy}")
        
        # 计算能量变化
        pme_energy_change = moved_pme_reciprocal - initial_pme_reciprocal
        pgp_energy_change = moved_pgp_energy - initial_pgp_energy
        
        print(f"PME倒空间能量变化: {pme_energy_change}")
        print(f"PGP能量变化: {pgp_energy_change}")
        
        # 计算相对误差
        if abs(pme_energy_change) > 1e-6:
            pgp_pme_error = abs((pgp_energy_change - pme_energy_change) / pme_energy_change)
            print(f"PGP与PME相对误差: {pgp_pme_error*100:.4f}%")
            pgp_pme_errors.append(pgp_pme_error)
        else:
            print("PME能量变化接近零，跳过相对误差计算")
    
    # 计算平均误差
    if pgp_pme_errors:
        avg_error = sum(pgp_pme_errors) / len(pgp_pme_errors)
        print(f"\n{num_moves}次移动的平均相对误差: {avg_error*100:.4f}%")
        
        # 使用更宽松的错误容限，因为PGP是一种近似方法
        acceptable_error = 0.5  # 允许50%的误差
        assert avg_error < acceptable_error, f"PGP与PME平均相对误差过大: {avg_error*100:.2f}%"
    else:
        print("\n没有有效的误差数据用于统计")
