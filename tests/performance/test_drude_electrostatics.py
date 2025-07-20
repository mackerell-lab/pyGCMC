#!/usr/bin/env python3
"""
测试PyGCMC的Drude是否考虑静电相互作用
不使用PME，使用简单的截断方法
"""

import numpy as np
import pygcmc

def test_drude_with_point_charge():
    """
    测试Drude对点电荷的响应
    """
    print("测试Drude对点电荷的响应")
    print("="*70)
    
    # 创建系统：一个水分子 + 一个点电荷
    state = pygcmc.MCState()
    
    box_size = 2.0  # nm - 大盒子避免周期性影响
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.cutoff = 0.9  # nm
    
    # SWM4-NDP水分子
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    
    # 水分子在盒子中心
    water_pos = [
        [1.0, 1.0, 1.0],      # O
        [1.0, 1.0, 1.0],      # D (初始在O位置)
        [1.096, 1.0, 1.0],    # H1
        [0.952, 1.077, 1.0],  # H2
        [1.0, 1.0, 1.0]       # M
    ]
    
    atoms = []
    
    # 添加水分子
    for i in range(5):
        atom = pygcmc.MCAtom()
        atom.x = water_pos[i][0]
        atom.y = water_pos[i][1]
        atom.z = water_pos[i][2]
        atom.charge = charges[i]
        atom.type = i if i < 4 else 3
        atoms.append(atom)
    
    # 添加一个点电荷 (作为第6个原子，但不是Drude粒子)
    point_charge = pygcmc.MCAtom()
    point_charge.x = 1.5  # 距离氧原子0.5 nm
    point_charge.y = 1.0
    point_charge.z = 1.0
    point_charge.charge = 5.0  # 强正电荷
    point_charge.type = 0  # 随便一个类型
    atoms.append(point_charge)
    
    # 创建残基
    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 5
    res1.active = True
    res1.type = 0
    
    res2 = pygcmc.MCResidue()
    res2.atomStart = 5
    res2.atomCount = 1
    res2.active = True
    res2.type = 1
    
    state.atoms = atoms
    state.residues = [res1, res2]
    state.activeAtomCount = 6
    state.activeResidueCount = 2
    
    # 力场参数 - 需要5x5的矩阵来覆盖所有原子类型
    state.forcefield.numTotalTypes = 5  # 增加到5以包含点电荷
    state.forcefield.numMovementTypes = 5
    
    # 创建5x5的LJ参数矩阵
    ljSigma = []
    ljEps = []
    for i in range(5):
        for j in range(5):
            if i == 0 and j == 0:  # O-O
                ljSigma.append(0.318395)
                ljEps.append(0.88257)
            else:
                ljSigma.append(0.0)
                ljEps.append(0.0)
    
    state.forcefield.ljSigma = ljSigma
    state.forcefield.ljEps = ljEps
    
    print(f"系统设置：")
    print(f"  水分子O原子位置: ({water_pos[0][0]}, {water_pos[0][1]}, {water_pos[0][2]})")
    print(f"  点电荷位置: ({point_charge.x}, {point_charge.y}, {point_charge.z})")
    print(f"  点电荷大小: +{point_charge.charge}")
    print(f"  距离: 0.5 nm")
    
    # 计算理论电场和预期位移
    ONE_4PI_EPS0 = 138.935456  # kJ·mol⁻¹·nm·e⁻²
    r = 0.5  # nm
    E_field = ONE_4PI_EPS0 * point_charge.charge / (r * r)
    print(f"\n理论计算：")
    print(f"  电场强度: {E_field:.1f} kJ/(mol·nm·e)")
    
    # 预期位移：F = q*E, 位移 = F/k
    q_drude = -1.71636
    alpha = 0.0009782237  # nm³
    k_spring = q_drude**2 * ONE_4PI_EPS0 / (2 * alpha)
    expected_disp = abs(q_drude) * E_field / k_spring * 1000  # pm
    print(f"  弹簧常数: {k_spring:.1f} kJ/(mol·nm²)")
    print(f"  预期Drude位移: {expected_disp:.2f} pm (朝向正电荷)")
    
    # 创建DrudeForce - 只对水分子
    force = pygcmc.DrudeForce()
    
    force.addParticle(
        drudeIndex=1,     # D
        parentIndex=0,    # O
        aniso1Index=-1,
        aniso2Index=-1,
        aniso3Index=-1,
        aniso4Index=-1,
        charge=-1.71636,
        polarizability=0.0009782237,
        aniso12=1.0,
        aniso34=1.0
    )
    
    # 不添加Thole对（点电荷不是Drude粒子）
    
    # SCF参数
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.1
    params.maxIterations = 500
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.1  # 100 pm
    force.setSCFParameters(params)
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    # 运行SCF
    print("\n运行Drude SCF优化...")
    state_opt = state.copy()
    
    try:
        energy = force.calculateEnergySCF(state_opt)
        print(f"  Drude能量: {energy:.4f} kJ/mol")
        
        # 检查Drude位移
        dx = state_opt.atoms[1].x - state_opt.atoms[0].x
        dy = state_opt.atoms[1].y - state_opt.atoms[0].y
        dz = state_opt.atoms[1].z - state_opt.atoms[0].z
        
        disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000  # pm
        disp_vec = np.array([dx, dy, dz]) * 1000  # pm
        
        print(f"\n实际Drude位移:")
        print(f"  位移大小: {disp:.2f} pm")
        print(f"  位移向量: ({disp_vec[0]:.2f}, {disp_vec[1]:.2f}, {disp_vec[2]:.2f}) pm")
        
        # 计算诱导偶极矩
        dipole_vec = np.array([dx, dy, dz]) * q_drude * 4.80321  # Debye
        dipole_mag = np.linalg.norm(dipole_vec)
        print(f"  诱导偶极矩: {dipole_mag:.3f} D")
        
        # 判断是否响应了电场
        if disp > 5.0 and disp_vec[0] > 0:  # 应该向+x方向移动
            print(f"\n✓ Drude响应了点电荷的电场！")
            print(f"  实际/预期 = {disp/expected_disp:.2f}")
        else:
            print(f"\n✗ Drude没有正确响应电场")
            
    except Exception as e:
        print(f"  SCF失败: {e}")

def test_drude_with_cutoff_energy():
    """
    使用截断库仑能量测试
    """
    print("\n\n使用截断库仑能量测试")
    print("="*70)
    
    # 创建简单系统
    state = pygcmc.MCState()
    
    box_size = 2.0
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.cutoff = 0.9
    
    # 两个带电粒子
    atoms = []
    
    # 粒子1: 正电荷
    atom1 = pygcmc.MCAtom()
    atom1.x = 1.0
    atom1.y = 1.0
    atom1.z = 1.0
    atom1.charge = 1.0
    atom1.type = 0
    atoms.append(atom1)
    
    # 粒子2: Drude粒子（负电荷）
    atom2 = pygcmc.MCAtom()
    atom2.x = 1.0  # 初始在粒子1位置
    atom2.y = 1.0
    atom2.z = 1.0
    atom2.charge = -1.0
    atom2.type = 1
    atoms.append(atom2)
    
    # 粒子3: 外部正电荷
    atom3 = pygcmc.MCAtom()
    atom3.x = 1.5
    atom3.y = 1.0
    atom3.z = 1.0
    atom3.charge = 2.0
    atom3.type = 0
    atoms.append(atom3)
    
    res = pygcmc.MCResidue()
    res.atomStart = 0
    res.atomCount = 3
    res.active = True
    res.type = 0
    
    state.atoms = atoms
    state.residues = [res]
    state.activeAtomCount = 3
    state.activeResidueCount = 1
    
    # 力场
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 2
    state.forcefield.ljSigma = [0.0, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.0, 0.0, 0.0, 0.0]
    
    print("测试系统：")
    print("  粒子1 (parent): q=+1, 位置=(1,1,1)")
    print("  粒子2 (Drude): q=-1, 初始位置=(1,1,1)")
    print("  粒子3 (外部): q=+2, 位置=(1.5,1,1)")
    
    # 计算截断库仑能量
    print("\n使用computeSystemEnergyCutoff计算能量...")
    try:
        result = pygcmc.computeSystemEnergyCutoff(state)
        if isinstance(result, tuple) and len(result) >= 2:
            elec_energy, vdw_energy = result[:2]
            print(f"  静电能量: {elec_energy:.4f} kJ/mol")
            print(f"  LJ能量: {vdw_energy:.4f} kJ/mol")
            print(f"  总能量: {elec_energy + vdw_energy:.4f} kJ/mol")
        else:
            print(f"  能量: {result}")
    except Exception as e:
        print(f"  计算失败: {e}")
    
    # 手动移动Drude粒子测试能量变化
    print("\n手动移动Drude粒子...")
    state2 = state.copy()
    state2.atoms[1].x = 1.05  # 向外部电荷移动5 pm
    
    try:
        result2 = pygcmc.computeSystemEnergyCutoff(state2)
        if isinstance(result2, tuple) and len(result2) >= 2:
            elec_energy2 = result2[0]
            print(f"  移动后静电能量: {elec_energy2:.4f} kJ/mol")
            print(f"  能量变化: {elec_energy2 - elec_energy:.4f} kJ/mol")
            
            if elec_energy2 < elec_energy:
                print("  ✓ 向正电荷移动降低了能量（合理）")
            else:
                print("  ✗ 能量变化不合理")
    except:
        pass

def main():
    """
    主函数
    """
    test_drude_with_point_charge()
    test_drude_with_cutoff_energy()
    
    print("\n\n结论：")
    print("="*70)
    print("测试PyGCMC的Drude是否考虑了静电相互作用")
    print("如果Drude正确响应点电荷，说明SCF包含了静电力")

if __name__ == "__main__":
    main()