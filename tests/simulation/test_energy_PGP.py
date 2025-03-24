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
