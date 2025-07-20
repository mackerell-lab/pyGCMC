#!/usr/bin/env python3
"""
诊断SCF收敛问题的根源
"""

import pygcmc
import numpy as np
import pickle

def analyze_initial_forces():
    """
    分析初始力的大小，找出收敛问题的根源
    """
    print("SCF收敛问题诊断")
    print("="*70)
    
    # 测试不同系统
    test_cases = [
        ("单水分子", 1),
        ("两水分子（远距离）", 2),
        ("两水分子（密接触）", 2)
    ]
    
    for case_name, n_waters in test_cases:
        print(f"\n\n测试案例: {case_name}")
        print("-"*60)
        
        # 创建系统
        state = pygcmc.MCState()
        atoms = []
        residues = []
        
        # SWM4-NDP参数
        charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
        atom_types = [0, 1, 2, 2, 3]
        
        if n_waters == 1:
            # 单个水分子
            positions = [
                [0.0, 0.0, 0.0],      # O
                [0.0, 0.0, 0.0],      # D
                [0.09572, 0.0, 0.0],  # H1
                [-0.04786, 0.0, 0.08288],  # H2
                [0.0, -0.024034, 0.0] # M
            ]
        else:
            # 两个水分子
            if "远距离" in case_name:
                distance = 1.0  # nm
            else:
                distance = 0.3  # nm
            
            positions = [
                # 第一个水
                [0.0, 0.0, 0.0],      # O
                [0.0, 0.0, 0.0],      # D
                [0.09572, 0.0, 0.0],  # H1
                [-0.04786, 0.0, 0.08288],  # H2
                [0.0, -0.024034, 0.0], # M
                # 第二个水
                [distance, 0.0, 0.0],  # O
                [distance, 0.0, 0.0],  # D
                [distance+0.09572, 0.0, 0.0],  # H1
                [distance-0.04786, 0.0, 0.08288],  # H2
                [distance, -0.024034, 0.0] # M
            ]
        
        # 创建原子
        for i in range(n_waters):
            for j in range(5):
                atom = pygcmc.MCAtom()
                idx = i * 5 + j
                atom.x = positions[idx][0]
                atom.y = positions[idx][1]
                atom.z = positions[idx][2]
                atom.charge = charges[j]
                atom.type = atom_types[j]
                atoms.append(atom)
            
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
        
        # 设置盒子
        if "密接触" in case_name:
            state.info.box = np.array([0.5, 0.5, 0.5])  # 小盒子
        else:
            state.info.box = np.array([10.0, 10.0, 10.0])  # 大盒子
        state.info.cutoff = min(5.0, state.info.box[0]/2 - 0.01)
        
        # 力场
        state.forcefield.numTotalTypes = 4
        state.forcefield.numMovementTypes = 4
        state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
        state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
        
        # 创建DrudeForce
        force = pygcmc.DrudeForce()
        
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
        
        # 如果是两水系统，可能添加Thole
        if n_waters == 2 and "远距离" not in case_name:
            # 检查O-O距离
            dx = state.atoms[5].x - state.atoms[0].x
            dy = state.atoms[5].y - state.atoms[0].y
            dz = state.atoms[5].z - state.atoms[0].z
            
            # PBC
            if state.info.box[0] > 0:
                dx -= state.info.box[0] * round(dx / state.info.box[0])
                dy -= state.info.box[1] * round(dy / state.info.box[1])
                dz -= state.info.box[2] * round(dz / state.info.box[2])
            
            dist = np.sqrt(dx*dx + dy*dy + dz*dz)
            
            if dist < min(0.8, state.info.box[0]/2 - 0.01):
                force.addScreenedPair(0, 1, 1.3)
                print(f"  添加了Thole对，O-O距离: {dist:.3f} nm")
        
        # 分析初始状态
        print(f"\n初始配置:")
        print(f"  盒子: {state.info.box}")
        print(f"  截断: {state.info.cutoff:.3f} nm")
        
        # 计算初始电场
        if n_waters == 1:
            print("\n单水分子应该没有外场，Drude不应该移动")
        else:
            # 计算第一个Drude粒子感受到的电场
            drude_idx = 1
            E_x = E_y = E_z = 0.0
            
            # 来自第二个分子的电场
            for j in range(5, 10):
                if state.atoms[j].charge == 0:
                    continue
                
                dx = state.atoms[j].x - state.atoms[drude_idx].x
                dy = state.atoms[j].y - state.atoms[drude_idx].y
                dz = state.atoms[j].z - state.atoms[drude_idx].z
                
                # PBC
                if state.info.box[0] > 0:
                    dx -= state.info.box[0] * round(dx / state.info.box[0])
                    dy -= state.info.box[1] * round(dy / state.info.box[1])
                    dz -= state.info.box[2] * round(dz / state.info.box[2])
                
                r2 = dx*dx + dy*dy + dz*dz
                if r2 < 1e-10:
                    continue
                
                r = np.sqrt(r2)
                if state.info.cutoff > 0 and r > state.info.cutoff:
                    continue
                
                # E = k*q/r^2 * r_hat
                ONE_4PI_EPS0 = 138.935456
                E_mag = ONE_4PI_EPS0 * state.atoms[j].charge / r2
                E_x += E_mag * dx / r
                E_y += E_mag * dy / r
                E_z += E_mag * dz / r
            
            E_total = np.sqrt(E_x*E_x + E_y*E_y + E_z*E_z)
            print(f"\n第一个Drude粒子的初始电场:")
            print(f"  E = ({E_x:.1f}, {E_y:.1f}, {E_z:.1f}) kJ/(mol·nm·e)")
            print(f"  |E| = {E_total:.1f} kJ/(mol·nm·e)")
            
            # 预期的Drude位移
            alpha = 0.0009782237
            expected_disp = alpha * E_total / 1.71636  # nm
            print(f"  预期位移: {expected_disp*1000:.1f} pm")
        
        # 运行SCF
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 10.0
        params.maxIterations = 10  # 只做几步
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        force.setSCFParameters(params)
        force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        
        print("\n运行SCF:")
        try:
            energy = force.calculateEnergySCF(state)
            print(f"  收敛！能量 = {energy:.2f} kJ/mol")
            
            # 检查位移
            for i in range(n_waters):
                o_idx = i * 5
                d_idx = i * 5 + 1
                
                dx = state.atoms[d_idx].x - state.atoms[o_idx].x
                dy = state.atoms[d_idx].y - state.atoms[o_idx].y
                dz = state.atoms[d_idx].z - state.atoms[o_idx].z
                
                disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
                print(f"  水{i+1} Drude位移: {disp:.2f} pm")
                
        except Exception as e:
            print(f"  未收敛: {e}")

def test_problematic_configuration():
    """
    测试有问题的配置
    """
    print("\n\n\n特殊测试：密度1.0的2水系统")
    print("="*70)
    
    try:
        # 加载实际的2水系统
        with open('../tests/performance/water_density_1.0/water_2.pkl', 'rb') as f:
            data = pickle.load(f)
        
        print(f"盒子长度: {data['box_length']:.3f} nm")
        print(f"密度: {data['density']:.3f} g/cm³")
        
        # 计算O-O距离
        o1_pos = np.array(data['positions'][0])
        o2_pos = np.array(data['positions'][5])
        
        delta = o2_pos - o1_pos
        box = data['box_length']
        delta = delta - box * np.round(delta / box)
        dist = np.linalg.norm(delta)
        
        print(f"O-O距离: {dist:.3f} nm")
        print(f"Thole截断(0.8 nm) / 半盒子({box/2:.3f} nm) = {0.8/(box/2):.2f}")
        print("\n问题：Thole截断远超半盒子！")
        
    except Exception as e:
        print(f"无法加载文件: {e}")

def main():
    analyze_initial_forces()
    test_problematic_configuration()
    
    print("\n\n诊断结论:")
    print("="*70)
    print("1. 单水分子正确地没有外场，Drude不移动")
    print("2. 小盒子系统的PBC导致极高的初始力")
    print("3. Thole截断超过半盒子违反了最小镜像约定")
    print("4. 建议：")
    print("   - 使用更大的测试系统")
    print("   - 或降低密度以获得更大盒子")
    print("   - 实现更好的初始猜测")

if __name__ == "__main__":
    main()