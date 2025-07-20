#!/usr/bin/env python3
"""
详细测试Drude粒子的力平衡
计算每个Drude粒子受到的净力
"""

import pygcmc
import numpy as np
import pickle
import os

def calculate_drude_forces(state, force):
    """
    计算所有Drude粒子受到的力
    """
    n_waters = state.activeResidueCount
    
    # 首先运行SCF优化Drude位置
    print("运行SCF优化Drude位置...")
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 10.0  # 较严格的容差
    params.maxIterations = 200
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    # 复制状态以保持原始位置
    optimized_state = state.copy()
    
    try:
        energy = force.calculateEnergySCF(optimized_state)
        print(f"SCF收敛，能量: {energy:.2f} kJ/mol")
        print(f"能量/水: {energy/n_waters:.2f} kJ/mol")
    except Exception as e:
        print(f"SCF未收敛: {e}")
        return None
    
    # 分析每个Drude粒子的力平衡
    print("\n分析Drude粒子力平衡...")
    
    drude_forces = []
    drude_displacements = []
    parent_child_forces = []
    
    # 常数
    k_drude = 418400.0  # kJ/mol/nm² (Drude弹簧常数)
    ONE_4PI_EPS0 = 138.935456  # kJ*nm/mol/e²
    
    for i in range(n_waters):
        o_idx = i * 5      # O原子
        d_idx = i * 5 + 1  # Drude粒子
        
        # 获取位置
        o_pos = np.array([optimized_state.atoms[o_idx].x, 
                         optimized_state.atoms[o_idx].y, 
                         optimized_state.atoms[o_idx].z])
        d_pos = np.array([optimized_state.atoms[d_idx].x, 
                         optimized_state.atoms[d_idx].y, 
                         optimized_state.atoms[d_idx].z])
        
        # 计算位移
        disp_vec = d_pos - o_pos
        disp = np.linalg.norm(disp_vec)
        drude_displacements.append(disp * 1000)  # nm -> pm
        
        # 1. 计算弹簧回复力 (指向parent)
        f_spring = -k_drude * disp_vec
        f_spring_mag = np.linalg.norm(f_spring)
        parent_child_forces.append(f_spring_mag)
        
        # 2. 计算Drude受到的总电场力
        # 需要考虑来自所有其他原子的库仑力
        f_electric = np.zeros(3)
        
        # Drude电荷
        q_drude = optimized_state.atoms[d_idx].charge
        
        # 遍历所有原子计算电场
        for j in range(optimized_state.activeAtomCount):
            if j == d_idx:  # 跳过自己
                continue
                
            # 同一残基内的原子不产生电场力（分子内排除）
            j_res = j // 5
            if j_res == i:
                continue
            
            # 其他原子的位置和电荷
            other_pos = np.array([optimized_state.atoms[j].x,
                                 optimized_state.atoms[j].y,
                                 optimized_state.atoms[j].z])
            q_other = optimized_state.atoms[j].charge
            
            # 计算距离向量（考虑PBC）
            delta = d_pos - other_pos
            box = state.info.box[0]  # 假设立方盒子
            delta = delta - box * np.round(delta / box)
            
            r2 = np.dot(delta, delta)
            if r2 < 0.01:  # 避免太近的距离
                continue
                
            r = np.sqrt(r2)
            
            # 库仑力: F = k * q1 * q2 * r_vec / r³
            f_coulomb = ONE_4PI_EPS0 * q_drude * q_other * delta / (r2 * r)
            f_electric += f_coulomb
        
        f_electric_mag = np.linalg.norm(f_electric)
        
        # 3. 计算净力
        f_net = f_spring + f_electric
        f_net_mag = np.linalg.norm(f_net)
        
        drude_forces.append({
            'spring_force': f_spring_mag,
            'electric_force': f_electric_mag,
            'net_force': f_net_mag,
            'displacement': disp * 1000  # pm
        })
    
    return drude_forces, drude_displacements

def analyze_force_balance(filename, description):
    """
    分析系统的Drude力平衡
    """
    print(f"\n{'='*70}")
    print(f"{description}")
    print(f"{'='*70}")
    
    if not os.path.exists(filename):
        print(f"文件不存在: {filename}")
        return None
        
    with open(filename, 'rb') as f:
        data = pickle.load(f)
    
    n_waters = data['n_waters']
    positions = data['positions']
    box_length = data['box_length']
    
    print(f"系统信息:")
    print(f"  水分子数: {n_waters}")
    print(f"  盒子长度: {box_length:.3f} nm")
    print(f"  密度: {data.get('density', 'N/A')} g/cm³")
    
    # 创建PyGCMC状态
    state = pygcmc.MCState()
    atoms = []
    residues = []
    
    charges = data['charges']
    atom_types = [0, 1, 2, 2, 3]
    
    for i in range(n_waters * 5):
        atom = pygcmc.MCAtom()
        atom.x = positions[i][0]
        atom.y = positions[i][1]
        atom.z = positions[i][2]
        atom.charge = charges[i % 5]
        atom.type = atom_types[i % 5]
        atoms.append(atom)
    
    for i in range(n_waters):
        res = pygcmc.MCResidue()
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = n_waters
    
    state.info.box = np.array([box_length, box_length, box_length])
    state.info.cutoff = min(0.9, box_length / 2 - 0.01)
    
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    # 创建DrudeForce
    force = pygcmc.DrudeForce()
    
    # 添加Drude粒子
    for i in range(n_waters):
        force.addParticle(
            drudeIndex=5*i+1,
            parentIndex=5*i,
            aniso1Index=-1,
            aniso2Index=-1,
            aniso3Index=-1,
            aniso4Index=-1,
            charge=-1.71636,
            polarizability=0.0009782237,
            aniso12=1.0,
            aniso34=1.0
        )
    
    # 添加适量Thole对
    print("\n添加Thole对...")
    n_thole_pairs = 0
    thole_cutoff = 0.8
    
    # 为了准确性，添加更多Thole对
    for i in range(min(100, n_waters)):
        o1_idx = i * 5
        for j in range(i+1, min(i+50, n_waters)):
            o2_idx = j * 5
            
            dx = state.atoms[o1_idx].x - state.atoms[o2_idx].x
            dy = state.atoms[o1_idx].y - state.atoms[o2_idx].y
            dz = state.atoms[o1_idx].z - state.atoms[o2_idx].z
            
            dx -= box_length * round(dx / box_length)
            dy -= box_length * round(dy / box_length)
            dz -= box_length * round(dz / box_length)
            
            dist = np.sqrt(dx*dx + dy*dy + dz*dz)
            
            if dist < thole_cutoff:
                force.addScreenedPair(i, j, 1.3)
                n_thole_pairs += 1
    
    print(f"  添加了 {n_thole_pairs} 个Thole对")
    
    # 计算力平衡
    result = calculate_drude_forces(state, force)
    
    if result is None:
        return None
    
    drude_forces, displacements = result
    
    # 统计分析
    print("\n力平衡统计:")
    print("-"*60)
    
    # 分析前100个Drude粒子
    sample_size = min(100, len(drude_forces))
    
    net_forces = [f['net_force'] for f in drude_forces[:sample_size]]
    spring_forces = [f['spring_force'] for f in drude_forces[:sample_size]]
    electric_forces = [f['electric_force'] for f in drude_forces[:sample_size]]
    
    print(f"\n净力统计 (kJ/mol/nm):")
    print(f"  平均: {np.mean(net_forces):.2f}")
    print(f"  标准差: {np.std(net_forces):.2f}")
    print(f"  最大: {np.max(net_forces):.2f}")
    print(f"  最小: {np.min(net_forces):.2f}")
    
    print(f"\n弹簧力统计 (kJ/mol/nm):")
    print(f"  平均: {np.mean(spring_forces):.2f}")
    print(f"  最大: {np.max(spring_forces):.2f}")
    
    print(f"\n电场力统计 (kJ/mol/nm):")
    print(f"  平均: {np.mean(electric_forces):.2f}")
    print(f"  最大: {np.max(electric_forces):.2f}")
    
    print(f"\nDrude位移统计 (pm):")
    print(f"  平均: {np.mean(displacements):.2f}")
    print(f"  标准差: {np.std(displacements):.2f}")
    print(f"  最大: {np.max(displacements):.2f}")
    print(f"  最小: {np.min(displacements):.2f}")
    
    # 力平衡质量评估
    print(f"\n力平衡质量评估:")
    well_balanced = sum(1 for f in net_forces if f < 10.0)
    moderately_balanced = sum(1 for f in net_forces if 10.0 <= f < 50.0)
    poorly_balanced = sum(1 for f in net_forces if f >= 50.0)
    
    print(f"  优秀 (净力<10): {well_balanced} ({well_balanced/sample_size*100:.1f}%)")
    print(f"  中等 (10≤净力<50): {moderately_balanced} ({moderately_balanced/sample_size*100:.1f}%)")
    print(f"  较差 (净力≥50): {poorly_balanced} ({poorly_balanced/sample_size*100:.1f}%)")
    
    # 判断整体平衡情况
    avg_net_force = np.mean(net_forces)
    if avg_net_force < 20.0:
        print(f"\n✓ Drude粒子整体力平衡良好 (平均净力 {avg_net_force:.1f} kJ/mol/nm)")
    else:
        print(f"\n✗ Drude粒子力平衡需要改进 (平均净力 {avg_net_force:.1f} kJ/mol/nm)")
    
    return {
        'avg_net_force': avg_net_force,
        'avg_displacement': np.mean(displacements),
        'well_balanced_ratio': well_balanced/sample_size
    }

def main():
    """
    主函数
    """
    print("Drude粒子力平衡详细测试")
    print("="*70)
    
    # 测试系统
    test_systems = [
        ('../tests/performance/large_water_systems/water_256.pkl', '256水 - 未优化'),
        ('../tests/performance/optimized_water_systems/water_256_nvt_simple.pkl', '256水 - NVT优化'),
    ]
    
    results = {}
    
    for filename, description in test_systems:
        result = analyze_force_balance(filename, description)
        if result:
            results[description] = result
    
    # 对比总结
    print(f"\n\n{'='*70}")
    print("力平衡对比总结")
    print("="*70)
    
    if len(results) == 2:
        print(f"\n{'系统':^20} {'平均净力(kJ/mol/nm)':^20} {'平均位移(pm)':^15} {'优秀比例':^12}")
        print("-"*70)
        
        for desc, res in results.items():
            print(f"{desc:20} {res['avg_net_force']:^20.1f} {res['avg_displacement']:^15.1f} "
                  f"{res['well_balanced_ratio']*100:^12.1f}%")
        
        print(f"\n结论:")
        
        # 判断哪个系统更好
        unopt = results.get('256水 - 未优化')
        opt = results.get('256水 - NVT优化')
        
        if unopt and opt:
            if opt['avg_net_force'] < unopt['avg_net_force']:
                print("1. ✓ NVT优化后的系统力平衡更好")
            else:
                print("1. ✗ 未优化系统的力平衡更好")
                
            if opt['avg_displacement'] < 20.0 and unopt['avg_displacement'] < 20.0:
                print("2. ✓ 两个系统的Drude位移都在合理范围内(<20pm)")
            else:
                print("2. ⚠ 某些系统的Drude位移偏大")
                
            if opt['well_balanced_ratio'] > 0.8 and unopt['well_balanced_ratio'] > 0.8:
                print("3. ✓ 两个系统都有超过80%的Drude粒子达到良好平衡")
            else:
                print("3. ⚠ 部分Drude粒子的力平衡需要改进")

if __name__ == "__main__":
    main()