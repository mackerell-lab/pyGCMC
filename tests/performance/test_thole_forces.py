#!/usr/bin/env python3
"""
测试Thole屏蔽力计算的正确性
通过数值导数验证解析力
"""

import pygcmc
import numpy as np

def create_two_water_system(distance=0.5):
    """
    创建两个水分子系统，用于测试Thole相互作用
    """
    state = pygcmc.MCState()
    atoms = []
    residues = []
    
    # SWM4-NDP参数
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]  # O, D, H1, H2, M
    atom_types = [0, 1, 2, 2, 3]
    
    # 第一个水分子在原点
    positions1 = [
        [0.0, 0.0, 0.0],      # O
        [0.0, 0.0, 0.0],      # D (初始与O重合)
        [0.09572, 0.0, 0.0],  # H1
        [-0.09572, 0.0, 0.0], # H2
        [0.0, 0.024034, 0.0]  # M
    ]
    
    # 第二个水分子在x轴上distance处
    positions2 = [
        [distance, 0.0, 0.0],          # O
        [distance, 0.0, 0.0],          # D (初始与O重合)
        [distance + 0.09572, 0.0, 0.0], # H1
        [distance - 0.09572, 0.0, 0.0], # H2
        [distance, 0.024034, 0.0]       # M
    ]
    
    # 创建原子
    for i, (pos, charge, atype) in enumerate(zip(positions1 + positions2, 
                                                  charges + charges, 
                                                  atom_types + atom_types)):
        atom = pygcmc.MCAtom()
        atom.x = pos[0]
        atom.y = pos[1]
        atom.z = pos[2]
        atom.charge = charge
        atom.type = atype
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
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = 2
    
    # 设置大盒子避免PBC
    state.info.box = np.array([10.0, 10.0, 10.0])
    state.info.cutoff = 5.0
    
    # 设置力场参数
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state

def numerical_force(force_obj, state, atom_idx, h=1e-5):
    """
    计算数值力（有限差分）
    """
    # 保存原始位置
    original_pos = [state.atoms[atom_idx].x, 
                    state.atoms[atom_idx].y, 
                    state.atoms[atom_idx].z]
    
    force_numerical = np.zeros(3)
    
    # 对每个坐标计算数值导数
    for i, coord in enumerate(['x', 'y', 'z']):
        # 正向扰动
        setattr(state.atoms[atom_idx], coord, original_pos[i] + h)
        energy_plus = force_obj.calculateEnergySCF(state)
        
        # 负向扰动
        setattr(state.atoms[atom_idx], coord, original_pos[i] - h)
        energy_minus = force_obj.calculateEnergySCF(state)
        
        # 数值力 = -dE/dx
        force_numerical[i] = -(energy_plus - energy_minus) / (2 * h)
        
        # 恢复原始位置
        setattr(state.atoms[atom_idx], coord, original_pos[i])
    
    return force_numerical

def test_thole_forces():
    """
    测试Thole力计算
    """
    print("测试Thole屏蔽力计算的正确性")
    print("="*60)
    
    # 创建两水分子系统
    distances = [0.4, 0.5, 0.6, 0.8]  # nm
    
    for dist in distances:
        print(f"\n测试距离: {dist} nm")
        print("-"*40)
        
        state = create_two_water_system(dist)
        
        # 创建DrudeForce
        force = pygcmc.DrudeForce()
        
        # 添加Drude粒子
        drude_charge = -1.71636
        polarizability = 0.0009782237
        
        for i in range(2):
            force.addParticle(
                drudeIndex=5*i+1,
                parentIndex=5*i,
                aniso1Index=-1,
                aniso2Index=-1,
                aniso3Index=-1,
                aniso4Index=-1,
                charge=drude_charge,
                polarizability=polarizability,
                aniso12=1.0,
                aniso34=1.0
            )
        
        # 添加Thole屏蔽对
        force.addScreenedPair(0, 1, 1.3)
        
        # 设置SCF参数
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e-6
        params.maxIterations = 200
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        force.setSCFParameters(params)
        force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        
        # 运行SCF获得平衡位置
        energy = force.calculateEnergySCF(state)
        print(f"能量: {energy:.6f} kJ/mol")
        
        # 测试关键原子的力
        test_atoms = [0, 1, 5, 6]  # O1, D1, O2, D2
        atom_names = ['O1', 'D1', 'O2', 'D2']
        
        print("\n力的验证 (kJ/mol/nm):")
        print(f"{'原子':<6} {'解析x':<12} {'数值x':<12} {'误差x(%)':<10}")
        print(f"{'':6} {'解析y':<12} {'数值y':<12} {'误差y(%)':<10}")
        print(f"{'':6} {'解析z':<12} {'数值z':<12} {'误差z(%)':<10}")
        
        max_error = 0.0
        
        for atom_idx, name in zip(test_atoms, atom_names):
            # 计算解析力
            forces_analytical = [pygcmc.Vec3() for _ in range(state.activeAtomCount)]
            _ = force.calculateScreenedCoulombEnergy(state, forces_analytical)
            _ = force.calculateHarmonicEnergy(state, forces_analytical)
            
            f_analytical = np.array([forces_analytical[atom_idx].x,
                                    forces_analytical[atom_idx].y,
                                    forces_analytical[atom_idx].z])
            
            # 计算数值力
            f_numerical = numerical_force(force, state, atom_idx)
            
            # 计算误差
            errors = np.abs((f_analytical - f_numerical) / (f_numerical + 1e-10)) * 100
            max_error = max(max_error, np.max(errors))
            
            print(f"{name:<6} {f_analytical[0]:12.6f} {f_numerical[0]:12.6f} {errors[0]:10.2f}")
            print(f"{'':6} {f_analytical[1]:12.6f} {f_numerical[1]:12.6f} {errors[1]:10.2f}")
            print(f"{'':6} {f_analytical[2]:12.6f} {f_numerical[2]:12.6f} {errors[2]:10.2f}")
        
        if max_error < 1.0:
            print(f"\n✓ 测试通过！最大误差: {max_error:.2f}%")
        else:
            print(f"\n✗ 测试失败！最大误差: {max_error:.2f}%")

def main():
    """
    主函数
    """
    test_thole_forces()
    
    print("\n\n结论:")
    print("如果所有测试通过（误差 < 1%），说明Thole力计算正确")
    print("如果测试失败，需要检查力的导数公式")

if __name__ == "__main__":
    main()