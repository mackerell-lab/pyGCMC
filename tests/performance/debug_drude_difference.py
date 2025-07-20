#!/usr/bin/env python3
"""
调试PyGCMC和OpenMM的Drude实现差异
"""

import numpy as np
import pygcmc

def test_simple_drude_system():
    """
    测试最简单的Drude系统
    """
    print("简单Drude系统调试")
    print("="*70)
    
    # 创建两个水分子系统
    state = pygcmc.MCState()
    
    box_size = 1.0  # nm
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.cutoff = 0.45
    
    # SWM4-NDP电荷
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    
    # 位置：两个水分子相距约0.52 nm
    positions = [
        # 水1
        [0.3, 0.3, 0.3],      # O
        [0.3, 0.3, 0.3],      # D (初始在O位置)
        [0.396, 0.3, 0.3],    # H1
        [0.252, 0.377, 0.3],  # H2
        [0.3, 0.3, 0.3],      # M
        # 水2  
        [0.7, 0.7, 0.7],      # O
        [0.7, 0.7, 0.7],      # D (初始在O位置)
        [0.796, 0.7, 0.7],    # H1
        [0.652, 0.777, 0.7],  # H2
        [0.7, 0.7, 0.7]       # M
    ]
    
    atoms = []
    residues = []
    
    for i in range(10):
        atom = pygcmc.MCAtom()
        atom.x = positions[i][0]
        atom.y = positions[i][1]
        atom.z = positions[i][2]
        atom.charge = charges[i % 5]
        atom.type = i % 5 if i % 5 < 4 else 3
        atoms.append(atom)
    
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
    
    # 力场参数（16个值的数组）
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    
    ljSigma = []
    ljEps = []
    for i in range(4):
        for j in range(4):
            if i == 0 and j == 0:
                ljSigma.append(0.318395)
                ljEps.append(0.88257)
            else:
                ljSigma.append(0.0)
                ljEps.append(0.0)
    
    state.forcefield.ljSigma = ljSigma
    state.forcefield.ljEps = ljEps
    
    # 1. 分析初始电场
    print("\n1. 手动计算电场")
    print("-"*60)
    
    ONE_4PI_EPS0 = 138.935456  # kJ·mol⁻¹·nm·e⁻²
    
    # 计算水1的O原子处的电场（来自水2）
    o1_pos = np.array(positions[0])
    field = np.zeros(3)
    
    # 只计算来自水2的带电粒子的贡献
    for j in [5, 7, 8, 9]:  # O2, H1, H2, M (跳过D2)
        other_pos = np.array(positions[j])
        q_other = charges[j % 5]
        
        delta = o1_pos - other_pos
        r2 = np.dot(delta, delta)
        r = np.sqrt(r2)
        
        E_mag = ONE_4PI_EPS0 * q_other / r2
        field += E_mag * delta / r
        
        print(f"  来自原子{j}(q={q_other:+.3f}): E_contribution = {E_mag * delta / r}")
    
    field_mag = np.linalg.norm(field)
    print(f"\n  总电场: E = ({field[0]:.2f}, {field[1]:.2f}, {field[2]:.2f}) kJ/(mol·nm·e)")
    print(f"  |E| = {field_mag:.2f} kJ/(mol·nm·e)")
    
    # 预期的Drude位移
    alpha = 0.0009782237  # nm³
    q_drude = -1.71636
    
    # F = q·E，位移 = F/k = q·E/k
    # 对于谐振子：k = q²/(2α)
    k_spring = q_drude**2 * ONE_4PI_EPS0 / (2 * alpha)
    expected_disp = abs(q_drude) * field_mag / k_spring
    
    print(f"\n  弹簧常数 k = {k_spring:.1f} kJ/(mol·nm²)")
    print(f"  预期Drude位移 = |q|·|E|/k = {expected_disp*1000:.2f} pm")
    
    # 2. 测试PyGCMC的DrudeForce
    print("\n\n2. PyGCMC DrudeForce测试")
    print("-"*60)
    
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
    
    # 设置SCF参数
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.1  # 宽松容差
    params.maxIterations = 500
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.05  # 50 pm
    force.setSCFParameters(params)
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    # 运行SCF
    state_opt = state.copy()
    
    print("  运行SCF优化...")
    try:
        energy_drude = force.calculateEnergySCF(state_opt)
        print(f"  Drude能量: {energy_drude:.6f} kJ/mol")
        
        # 检查Drude位移
        print("\n  Drude位移:")
        total_disp = 0
        for i in range(2):
            o_idx = i * 5
            d_idx = i * 5 + 1
            
            dx = state_opt.atoms[d_idx].x - state_opt.atoms[o_idx].x
            dy = state_opt.atoms[d_idx].y - state_opt.atoms[o_idx].y
            dz = state_opt.atoms[d_idx].z - state_opt.atoms[o_idx].z
            
            disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000  # pm
            total_disp += disp
            
            print(f"    水{i+1}: {disp:.4f} pm, 方向=({dx*1000:.3f}, {dy*1000:.3f}, {dz*1000:.3f}) pm")
        
        avg_disp = total_disp / 2
        print(f"\n  平均位移: {avg_disp:.4f} pm")
        
        if avg_disp < 0.01:
            print("\n  ⚠️ Drude位移几乎为0！")
            print("  可能的原因：")
            print("  1. SCF没有正确考虑外电场")
            print("  2. 电场计算有误")
            print("  3. 力常数设置不当")
            
    except Exception as e:
        print(f"  SCF失败: {e}")
    
    # 3. 调试SCF细节
    print("\n\n3. 调试SCF实现")
    print("-"*60)
    
    # 尝试不同的参数
    print("  尝试不同的SCF参数...")
    
    # 更激进的参数
    params2 = pygcmc.DrudeSCFParams()
    params2.tolerance = 1000.0  # 非常宽松
    params2.maxIterations = 10
    params2.dampingFactor = 1.0  # 无阻尼
    params2.maxDrudeDistance = 0.5  # 500 pm
    force.setSCFParameters(params2)
    
    state_opt2 = state.copy()
    
    # 手动移动Drude粒子看看会发生什么
    print("\n  手动移动Drude粒子测试...")
    state_opt2.atoms[1].x += 0.01  # 移动10 pm
    
    energy_before = force.calculateEnergySCF(state.copy())
    energy_after = force.calculateEnergySCF(state_opt2)
    
    print(f"  移动前能量: {energy_before:.6f} kJ/mol")
    print(f"  移动后能量: {energy_after:.6f} kJ/mol")
    print(f"  能量差: {energy_after - energy_before:.6f} kJ/mol")
    
    if abs(energy_after - energy_before) < 1e-6:
        print("\n  ⚠️ 移动Drude粒子后能量没有变化！")
        print("  这表明能量计算可能没有包含Drude贡献")

def main():
    """
    主函数
    """
    test_simple_drude_system()
    
    print("\n\n最终结论:")
    print("="*70)
    print("需要检查PyGCMC的DrudeForce实现：")
    print("1. 是否正确计算了外电场")
    print("2. 是否正确应用了电场力")
    print("3. SCF迭代是否正确更新位置")
    print("4. 能量计算是否包含所有项")

if __name__ == "__main__":
    main()