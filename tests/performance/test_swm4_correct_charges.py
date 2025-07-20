#!/usr/bin/env python3
"""
使用正确的SWM4-NDP电荷测试
"""

import pygcmc
import numpy as np
import pickle

def test_correct_swm4_charges():
    """
    测试正确的SWM4-NDP电荷设置
    """
    print("测试正确的SWM4-NDP水模型电荷")
    print("="*70)
    
    # 正确的SWM4-NDP电荷（来自PSF文件）
    charges = {
        'O': 1.71636,    # 氧核心（正电荷）
        'D': -1.71636,   # Drude粒子（负电荷）
        'H1': 0.55733,   # 氢1
        'H2': 0.55733,   # 氢2
        'M': -1.11466    # 虚拟位点M
    }
    
    print("\nSWM4-NDP电荷分配:")
    print("-"*40)
    for atom, charge in charges.items():
        print(f"  {atom:3}: {charge:+8.5f} e")
    
    total = sum(charges.values())
    print(f"  总计: {total:+8.5f} e (应该是0)")
    
    print("\n关键理解:")
    print("1. O原子电荷是+1.71636（不是0或3.43272）")
    print("2. Drude粒子电荷是-1.71636")
    print("3. O和D形成偶极：O(+) --- D(-)")
    print("4. 总电荷必须为0")
    
    # 加载之前生成的水系统测试
    print("\n\n测试实际系统:")
    print("-"*40)
    
    try:
        # 加载2水系统
        with open('../tests/performance/water_density_1.0/water_2.pkl', 'rb') as f:
            data = pickle.load(f)
        
        # 创建状态
        state = pygcmc.MCState()
        atoms = []
        residues = []
        
        # 使用正确的电荷
        correct_charges = [charges['O'], charges['D'], charges['H1'], 
                          charges['H2'], charges['M']]
        atom_types = [0, 1, 2, 2, 3]
        
        n_waters = 2
        positions = data['positions']
        
        # 创建原子
        for i in range(n_waters):
            for j in range(5):
                atom = pygcmc.MCAtom()
                idx = i * 5 + j
                atom.x = positions[idx][0]
                atom.y = positions[idx][1]
                atom.z = positions[idx][2]
                atom.charge = correct_charges[j]
                atom.type = atom_types[j]
                atoms.append(atom)
            
            # 创建残基
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
        state.info.box = np.array([data['box_length']] * 3)
        state.info.cutoff = min(1.2, data['box_length'] / 2 - 0.01)
        
        # 设置力场参数
        state.forcefield.numTotalTypes = 4
        state.forcefield.numMovementTypes = 4
        state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]  # 只有O有LJ
        state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
        
        # 验证电荷
        print(f"\n2水系统电荷验证:")
        for i in range(n_waters):
            total_charge = 0.0
            print(f"\n水分子 {i+1}:")
            atom_names = ['O', 'D', 'H1', 'H2', 'M']
            for j in range(5):
                idx = i * 5 + j
                charge = state.atoms[idx].charge
                total_charge += charge
                print(f"  {atom_names[j]:3}: {charge:+8.5f} e")
            print(f"  总计: {total_charge:+8.5f} e")
        
        # 创建DrudeForce测试
        force = pygcmc.DrudeForce()
        
        # 添加Drude粒子（使用正确的Drude电荷）
        drude_charge = charges['D']  # -1.71636
        polarizability = 0.0009782237
        
        for i in range(n_waters):
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
        
        # 不添加Thole（因为距离太近会导致错误）
        
        # SCF参数
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 10.0
        params.maxIterations = 100
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        force.setSCFParameters(params)
        force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        
        # 计算能量
        energy = force.calculateEnergySCF(state)
        print(f"\n\n能量计算:")
        print(f"  总能量: {energy:.2f} kJ/mol")
        print(f"  每水能量: {energy/n_waters:.2f} kJ/mol")
        
    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()

def main():
    test_correct_swm4_charges()
    
    print("\n\n结论:")
    print("1. 必须使用PSF文件中的正确电荷")
    print("2. 水分子总电荷必须为0")
    print("3. Drude模型中O原子有正电荷，D粒子有负电荷")

if __name__ == "__main__":
    main()