#!/usr/bin/env python3
"""
详细测试Thole实现，包括与解析解的比较
"""

import pygcmc
import numpy as np

def analytical_thole_energy(r, q_drude, alpha1, alpha2, thole_param=1.3):
    """
    计算两个Drude偶极之间的解析Thole能量
    
    使用OpenMM约定：
    - Drude粒子电荷 q_drude (负)
    - 母原子隐含电荷 -q_drude (正)
    """
    # 物理常数
    ONE_4PI_EPS0 = 138.935456  # kJ/mol·nm·e^-2
    
    # 屏蔽参数
    uscale = thole_param / (alpha1 * alpha2)**(1.0/6.0)
    u = r * uscale
    
    # 屏蔽函数
    exp_u = np.exp(-u)
    screening = 1.0 - (1.0 + 0.5 * u) * exp_u
    
    # 四种相互作用的能量
    # DD: q_drude * q_drude (两个负电荷)
    # DP: q_drude * (-q_drude) (负正)
    # PD: (-q_drude) * q_drude (正负)
    # PP: (-q_drude) * (-q_drude) (两个正电荷)
    
    q2 = q_drude * q_drude  # q_drude是负的，所以q2是正的
    
    # 假设偶极完全对齐（最坏情况）
    # DD距离 = r, DP距离 = r, PD距离 = r, PP距离 = r
    # 这是简化的情况，实际上偶极会有方向
    
    energy = ONE_4PI_EPS0 * q2 * screening / r * (1 - 1 - 1 + 1)
    # 注意：对于完全重合的偶极，净能量为0
    # 实际计算中偶极会分离，产生非零能量
    
    return energy

def test_thole_implementation():
    """
    详细测试Thole实现
    """
    print("详细测试Thole屏蔽实现")
    print("="*70)
    
    # 创建简单的两粒子系统
    state = pygcmc.MCState()
    atoms = []
    
    # 只创建两个相互作用的偶极（O-D对）
    # 第一个偶极
    o1 = pygcmc.MCAtom()
    o1.x, o1.y, o1.z = 0.0, 0.0, 0.0
    o1.charge = 0.0  # 暂时设为0，稍后讨论
    o1.type = 0
    atoms.append(o1)
    
    d1 = pygcmc.MCAtom()
    d1.x, d1.y, d1.z = 0.0, 0.0, 0.0  # 初始与O重合
    d1.charge = -1.71636  # Drude电荷
    d1.type = 1
    atoms.append(d1)
    
    # 第二个偶极
    o2 = pygcmc.MCAtom()
    o2.x, o2.y, o2.z = 0.5, 0.0, 0.0  # 距离0.5 nm
    o2.charge = 0.0
    o2.type = 0
    atoms.append(o2)
    
    d2 = pygcmc.MCAtom()
    d2.x, d2.y, d2.z = 0.5, 0.0, 0.0
    d2.charge = -1.71636
    d2.type = 1
    atoms.append(d2)
    
    # 创建残基
    residues = []
    for i in range(2):
        res = pygcmc.MCResidue()
        res.atomStart = 2 * i
        res.atomCount = 2
        res.active = True
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = 4
    state.activeResidueCount = 2
    
    # 大盒子
    state.info.box = np.array([10.0, 10.0, 10.0])
    state.info.cutoff = 5.0
    
    # 力场（简化）
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 2
    state.forcefield.ljSigma = [0.0, 0.0]
    state.forcefield.ljEps = [0.0, 0.0]
    
    # 创建DrudeForce
    force = pygcmc.DrudeForce()
    
    # 添加Drude粒子
    drude_charge = -1.71636
    polarizability = 0.0009782237
    
    force.addParticle(
        drudeIndex=1, parentIndex=0,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=drude_charge,
        polarizability=polarizability,
        aniso12=1.0, aniso34=1.0
    )
    
    force.addParticle(
        drudeIndex=3, parentIndex=2,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=drude_charge,
        polarizability=polarizability,
        aniso12=1.0, aniso34=1.0
    )
    
    # 添加Thole屏蔽
    force.addScreenedPair(0, 1, 1.3)
    
    # 设置SCF参数
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-8
    params.maxIterations = 1  # 只做一次，因为初始位置就是平衡位置
    params.dampingFactor = 0.0
    params.maxDrudeDistance = 1.0  # 大值，不限制
    force.setSCFParameters(params)
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    print("\n测试不同O原子电荷设置下的能量:")
    print("-"*60)
    print(f"{'O电荷':>10} {'总电荷':>10} {'Thole能量':>15} {'说明':>20}")
    print("-"*60)
    
    # 测试不同的电荷设置
    o_charges = [0.0, 1.71636, 3.43272]
    descriptions = ["错误设置", "原始电荷", "OpenMM风格"]
    
    for o_charge, desc in zip(o_charges, descriptions):
        # 设置O原子电荷
        state.atoms[0].charge = o_charge
        state.atoms[2].charge = o_charge
        
        # 计算总电荷
        total_charge = o_charge + drude_charge
        
        try:
            energy = force.calculateEnergySCF(state)
            print(f"{o_charge:10.5f} {total_charge:10.5f} {energy:15.6f} {desc:>20}")
        except Exception as e:
            print(f"{o_charge:10.5f} {total_charge:10.5f} {'错误':>15} {str(e):>20}")
    
    print("\n结论:")
    print("1. O电荷=0时，只有Drude-Drude相互作用，能量应该很小")
    print("2. O电荷=1.71636时，总电荷=0，这是物理正确的设置")
    print("3. O电荷=3.43272时，遵循OpenMM约定，但总电荷≠0")
    print("\n注意：如果我们的实现正确遵循OpenMM，那么:")
    print("- Thole能量应该只依赖于Drude电荷（通过charge参数传入）")
    print("- 不应该依赖于state.atoms中的电荷")

def main():
    test_thole_implementation()

if __name__ == "__main__":
    main()