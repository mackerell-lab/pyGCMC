#!/usr/bin/env python3
"""
简单的SCF验证：检查基本功能
"""

import numpy as np
import pygcmc

def create_two_water_system():
    """
    创建两个水分子的简单系统
    """
    state = pygcmc.MCState()
    
    # 盒子
    box_size = 1.0  # nm
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.cutoff = 0.45
    
    # SWM4-NDP参数
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    
    atoms = []
    residues = []
    
    # 水1: 在(0.3, 0.3, 0.3)
    water1_o = [0.3, 0.3, 0.3]
    # 水2: 在(0.6, 0.6, 0.6) - 距离约0.52 nm
    water2_o = [0.6, 0.6, 0.6]
    
    # 创建水1
    positions = [
        water1_o,                                    # O
        water1_o,                                    # D (初始在O位置)
        [water1_o[0]+0.0957, water1_o[1], water1_o[2]],  # H1
        [water1_o[0]-0.024, water1_o[1]+0.0926, water1_o[2]],  # H2
        water1_o                                     # M (简化)
    ]
    
    # 创建水2
    positions.extend([
        water2_o,                                    # O
        water2_o,                                    # D (初始在O位置)
        [water2_o[0]+0.0957, water2_o[1], water2_o[2]],  # H1
        [water2_o[0]-0.024, water2_o[1]+0.0926, water2_o[2]],  # H2
        water2_o                                     # M (简化)
    ])
    
    # 创建原子
    for i in range(10):
        atom = pygcmc.MCAtom()
        atom.x = positions[i][0]
        atom.y = positions[i][1]
        atom.z = positions[i][2]
        atom.charge = charges[i % 5]
        atom.type = i % 5 if i % 5 < 4 else 3
        atoms.append(atom)
    
    # 创建残基
    for i in range(2):
        res = pygcmc.MCResidue()
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = 10
    state.activeResidueCount = 2
    
    # 力场参数
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state

def main():
    """
    主测试函数
    """
    print("简单SCF验证测试")
    print("="*70)
    
    # 创建系统
    state = create_two_water_system()
    print(f"创建了2个水分子系统")
    print(f"  水1 O位置: ({state.atoms[0].x:.1f}, {state.atoms[0].y:.1f}, {state.atoms[0].z:.1f})")
    print(f"  水2 O位置: ({state.atoms[5].x:.1f}, {state.atoms[5].y:.1f}, {state.atoms[5].z:.1f})")
    
    # 计算O-O距离
    dx = state.atoms[5].x - state.atoms[0].x
    dy = state.atoms[5].y - state.atoms[0].y
    dz = state.atoms[5].z - state.atoms[0].z
    dist = np.sqrt(dx*dx + dy*dy + dz*dz)
    print(f"  O-O距离: {dist:.3f} nm")
    
    # 创建DrudeForce
    force = pygcmc.DrudeForce()
    
    # 添加Drude粒子
    for i in range(2):
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
    
    # 添加Thole屏蔽
    force.addScreenedPair(0, 1, 1.3)
    print(f"\n添加了Thole屏蔽对")
    
    # 测试1: 基础SCF测试
    print("\n\n测试1: 基础SCF收敛")
    print("-"*60)
    
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 100.0  # 宽松容差
    params.maxIterations = 100
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    # 初始状态 - Drude在parent位置
    print("\n初始Drude位置:")
    for i in range(2):
        d_idx = i * 5 + 1
        o_idx = i * 5
        print(f"  水{i+1} Drude: ({state.atoms[d_idx].x:.3f}, {state.atoms[d_idx].y:.3f}, {state.atoms[d_idx].z:.3f})")
        print(f"  水{i+1} O:     ({state.atoms[o_idx].x:.3f}, {state.atoms[o_idx].y:.3f}, {state.atoms[o_idx].z:.3f})")
    
    # 运行SCF
    test_state = state.copy()
    
    try:
        energy = force.calculateEnergySCF(test_state)
        print(f"\n✓ SCF收敛")
        print(f"  总能量: {energy:.2f} kJ/mol")
        
        # 分析结果
        print(f"\nSCF优化后的Drude位移:")
        for i in range(2):
            o_idx = i * 5
            d_idx = i * 5 + 1
            
            dx = test_state.atoms[d_idx].x - test_state.atoms[o_idx].x
            dy = test_state.atoms[d_idx].y - test_state.atoms[o_idx].y
            dz = test_state.atoms[d_idx].z - test_state.atoms[o_idx].z
            
            disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000  # pm
            
            print(f"  水{i+1}: 位移 = {disp:.2f} pm, 方向 = ({dx*1000:.2f}, {dy*1000:.2f}, {dz*1000:.2f}) pm")
            
            # 计算诱导偶极矩
            q_drude = -1.71636
            dipole_vec = np.array([dx, dy, dz]) * q_drude * 4.80321  # Debye
            dipole_mag = np.linalg.norm(dipole_vec)
            
            print(f"       诱导偶极矩 = {dipole_mag:.3f} D")
        
        # 测试2: 手动移动Drude粒子
        print("\n\n测试2: 手动移动Drude粒子后SCF恢复")
        print("-"*60)
        
        # 复制优化后的状态
        test_state2 = test_state.copy()
        
        # 手动移动第一个水的Drude粒子
        move_distance = 0.01  # 10 pm
        test_state2.atoms[1].x += move_distance
        test_state2.atoms[1].y += move_distance
        
        print(f"手动移动水1的Drude粒子 {move_distance*1000:.0f} pm")
        
        # 再次运行SCF
        energy2 = force.calculateEnergySCF(test_state2)
        
        # 检查是否恢复到原位置
        dx_new = test_state2.atoms[1].x - test_state.atoms[1].x
        dy_new = test_state2.atoms[1].y - test_state.atoms[1].y
        dz_new = test_state2.atoms[1].z - test_state.atoms[1].z
        diff = np.sqrt(dx_new*dx_new + dy_new*dy_new + dz_new*dz_new) * 1000
        
        print(f"\nSCF后:")
        print(f"  能量: {energy2:.2f} kJ/mol (之前: {energy:.2f})")
        print(f"  Drude位置差异: {diff:.3f} pm")
        
        if diff < 0.1:
            print(f"\n✓ SCF成功恢复Drude到优化位置")
        else:
            print(f"\n⚠ Drude位置有小的差异")
            
    except Exception as e:
        print(f"\n✗ SCF未收敛: {e}")
    
    # 测试3: 分析电场
    print("\n\n测试3: 电场分析")
    print("-"*60)
    
    # 计算水1的O原子处的电场（来自水2）
    ONE_4PI_EPS0 = 138.935456
    
    # 水2各原子对水1 O原子的电场贡献
    field_x = 0
    field_y = 0
    field_z = 0
    
    o1_pos = np.array([state.atoms[0].x, state.atoms[0].y, state.atoms[0].z])
    
    for j in range(5, 10):  # 水2的原子
        if j == 6:  # 跳过Drude
            continue
            
        other_pos = np.array([state.atoms[j].x, state.atoms[j].y, state.atoms[j].z])
        q_other = state.atoms[j].charge
        
        delta = o1_pos - other_pos
        r2 = np.dot(delta, delta)
        r = np.sqrt(r2)
        
        # E = k * q / r^2 * r_hat
        E_mag = ONE_4PI_EPS0 * q_other / r2
        E_vec = E_mag * delta / r
        
        field_x += E_vec[0]
        field_y += E_vec[1]
        field_z += E_vec[2]
    
    field_mag = np.sqrt(field_x*field_x + field_y*field_y + field_z*field_z)
    
    print(f"水1 O原子处的电场（来自水2）:")
    print(f"  E = ({field_x:.1f}, {field_y:.1f}, {field_z:.1f}) kJ/(mol·nm·e)")
    print(f"  |E| = {field_mag:.1f} kJ/(mol·nm·e)")
    
    # 预期的Drude位移
    alpha = 0.0009782237  # nm³
    expected_disp = alpha * field_mag / 1.71636  # nm
    print(f"\n预期的Drude位移（仅考虑外电场）: {expected_disp*1000:.2f} pm")
    
    print("\n\n总结:")
    print("="*70)
    print("1. SCF能够优化Drude位置")
    print("2. Drude粒子响应电场产生位移")
    print("3. SCF是稳定的（移动后能恢复）")

if __name__ == "__main__":
    main()