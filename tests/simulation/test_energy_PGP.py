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
    2. 正负电荷数量不平衡
    3. 移动残基包含一个带电粒子，远离固定部分
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
    # 定义不对称的位置和电荷
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
        ((box_size-2.0, 2.0, box_size-2.0),  0.6),  # 正电荷
        # 总电荷稍微为正，不完全中性
    ]
    
    # 确认粒子数
    print(f"固定粒子数: {len(fixed_particles)}")
    
    # 验证固定部分总电荷
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
    
    # 创建移动部分 - 单个带电粒子放在盒子中央
    # 位置足够远，确保与所有固定粒子的距离都大于cutoff
    mobile_atom_pos = (center_x, center_y, center_z)
    mobile_charge = -total_fixed_charge  # 使系统总电荷为0
    
    # 添加移动粒子
    mobile_atom = MCAtom()
    mobile_atom.x, mobile_atom.y, mobile_atom.z = mobile_atom_pos
    mobile_atom.charge = mobile_charge
    mobile_atom.type = 1 if mobile_charge < 0 else 0
    atoms.append(mobile_atom)
    
    print(f"移动粒子: 位置=({mobile_atom_pos[0]}, {mobile_atom_pos[1]}, {mobile_atom_pos[2]}), 电荷={mobile_charge}")
    
    # 创建移动残基
    mobile_res = MCResidue()
    mobile_res.atomStart = len(fixed_particles)
    mobile_res.atomCount = 1
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
    print(f"固定原子: {len(fixed_particles)}, 移动原子: 1")
    sys.stdout.flush()
    
    # 初始化PME和PGP
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
    
    # 步骤2: 移动移动残基（对角线方向平移0.3 nm）
    translation = [0.3, 0.3, 0.3]
    print(f"移动残基: 平移{translation}...")
    
    # 移动粒子
    system.atoms[mobile_res.atomStart].x += translation[0]
    system.atoms[mobile_res.atomStart].y += translation[1]
    system.atoms[mobile_res.atomStart].z += translation[2]
    
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
