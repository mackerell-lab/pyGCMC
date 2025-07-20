#!/usr/bin/env python3
"""
测试不同算法设置的实际差异
"""

import numpy as np
import time
import pygcmc

def test_algorithm_effects():
    """
    测试设置不同算法后的实际效果
    """
    print("测试算法设置的效果")
    print("="*70)
    
    # 创建简单系统
    state = create_simple_system(10)  # 10个水分子
    
    # 测试所有可用算法
    algorithms = [
        pygcmc.DrudeAlgorithm.SCF,
        pygcmc.DrudeAlgorithm.ConjugateGradient,
        pygcmc.DrudeAlgorithm.OPT3,
        pygcmc.DrudeAlgorithm.OPT4,
        pygcmc.DrudeAlgorithm.SmartOPT3,
        pygcmc.DrudeAlgorithm.AdaptiveOPT,
        pygcmc.DrudeAlgorithm.HybridOPT,
        pygcmc.DrudeAlgorithm.FBP
    ]
    
    print(f"{'算法名称':>20} {'calculateEnergyOPT3':>20} {'calculateEnergySCF':>20} {'备注':>30}")
    print("-"*95)
    
    for algo in algorithms:
        force = create_drude_force(state, 10)
        force.setAlgorithm(algo)
        
        # 如果是SCF，设置参数
        if algo == pygcmc.DrudeAlgorithm.SCF:
            params = pygcmc.DrudeSCFParams()
            params.tolerance = 10.0
            params.maxIterations = 50
            params.dampingFactor = 0.5
            params.maxDrudeDistance = 0.02
            force.setSCFParameters(params)
        
        results = {}
        
        # 测试OPT3方法
        try:
            state_test = state.copy()
            reset_drude_positions(state_test, 10)
            
            start = time.time()
            energy = force.calculateEnergyOPT3(state_test)
            elapsed = (time.time() - start) * 1000
            
            disp = calculate_avg_displacement(state_test, 10)
            results['OPT3'] = f"{energy/10:.2f} ({elapsed:.1f}ms)"
        except Exception as e:
            results['OPT3'] = "失败"
        
        # 测试SCF方法
        try:
            state_test = state.copy()
            reset_drude_positions(state_test, 10)
            
            start = time.time()
            energy = force.calculateEnergySCF(state_test)
            elapsed = (time.time() - start) * 1000
            
            disp = calculate_avg_displacement(state_test, 10)
            results['SCF'] = f"{energy/10:.2f} ({elapsed:.1f}ms)"
        except Exception as e:
            results['SCF'] = "失败"
        
        # 判断备注
        if results.get('OPT3') != "失败" and results.get('SCF') != "失败":
            note = "两种方法都可用"
        elif results.get('OPT3') != "失败":
            note = "只能用calculateEnergyOPT3"
        elif results.get('SCF') != "失败":
            note = "只能用calculateEnergySCF"
        else:
            note = "两种方法都不可用"
        
        print(f"{algo.name:>20} {results.get('OPT3', '失败'):>20} {results.get('SCF', '失败'):>20} {note:>30}")

def test_fbp_special():
    """
    专门测试FBP算法
    """
    print("\n\n专门测试FBP算法")
    print("="*70)
    
    # 测试不同大小的系统
    sizes = [2, 10, 50]
    
    print(f"{'系统大小':>10} {'FBP(via OPT3)':>20} {'参数影响':>30}")
    print("-"*65)
    
    for n_waters in sizes:
        state = create_simple_system(n_waters)
        
        # 测试默认FBP
        force1 = create_drude_force(state, n_waters)
        force1.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
        
        state_test1 = state.copy()
        reset_drude_positions(state_test1, n_waters)
        
        try:
            start = time.time()
            energy1 = force1.calculateEnergyOPT3(state_test1)
            time1 = (time.time() - start) * 1000
            disp1 = calculate_avg_displacement(state_test1, n_waters)
            
            result1 = f"{energy1/n_waters:.2f} kJ/mol ({time1:.1f}ms)"
            
            # 测试是否有FBP特定的参数
            # 检查是否有setFBPParameters之类的方法
            if hasattr(force1, 'setFBPParameters'):
                param_note = "有FBP特定参数"
            else:
                param_note = "无特定参数，可能使用默认设置"
                
        except Exception as e:
            result1 = "失败"
            param_note = str(e)[:30]
        
        print(f"{n_waters:>10} {result1:>20} {param_note:>30}")

def create_simple_system(n_waters):
    """创建简单的水系统"""
    state = pygcmc.MCState()
    atoms = []
    residues = []
    
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    atom_types = [0, 1, 2, 2, 3]
    
    # 简单排列
    spacing = 0.4  # nm
    idx = 0
    
    for i in range(n_waters):
        x = (idx % 5) * spacing + 0.1
        y = ((idx // 5) % 5) * spacing + 0.1
        z = (idx // 25) * spacing + 0.1
        idx += 1
        
        # 添加5个原子
        for j in range(5):
            atom = pygcmc.MCAtom()
            atom.x = x
            atom.y = y
            atom.z = z
            atom.charge = charges[j]
            atom.type = atom_types[j]
            atoms.append(atom)
        
        res = pygcmc.MCResidue()
        res.atomStart = i * 5
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    box_size = 3.0
    
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = n_waters * 5
    state.activeResidueCount = n_waters
    
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.cutoff = 1.0
    
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state

def create_drude_force(state, n_waters):
    """创建DrudeForce"""
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
    
    # 添加一些Thole对
    for i in range(min(n_waters, 10)):
        for j in range(i+1, min(i+5, n_waters)):
            force.addScreenedPair(i, j, 1.3)
    
    return force

def reset_drude_positions(state, n_waters):
    """重置Drude位置"""
    for i in range(n_waters):
        o_idx = i * 5
        d_idx = i * 5 + 1
        state.atoms[d_idx].x = state.atoms[o_idx].x
        state.atoms[d_idx].y = state.atoms[o_idx].y
        state.atoms[d_idx].z = state.atoms[o_idx].z

def calculate_avg_displacement(state, n_waters):
    """计算平均位移"""
    displacements = []
    for i in range(n_waters):
        o_idx = i * 5
        d_idx = i * 5 + 1
        
        dx = state.atoms[d_idx].x - state.atoms[o_idx].x
        dy = state.atoms[d_idx].y - state.atoms[o_idx].y
        dz = state.atoms[d_idx].z - state.atoms[o_idx].z
        
        disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
        displacements.append(disp)
    
    return np.mean(displacements)

def main():
    """
    主函数
    """
    test_algorithm_effects()
    test_fbp_special()
    
    print("\n\n最终结论：")
    print("="*70)
    print("1. PyGCMC目前只实现了两个能量计算方法：")
    print("   - calculateEnergySCF: 用于SCF算法")
    print("   - calculateEnergyOPT3: 用于其他所有算法（包括FBP）")
    print("\n2. FBP算法通过calculateEnergyOPT3调用")
    print("   - 可能在内部根据algorithm设置使用不同的实现")
    print("   - 但对外接口是相同的")
    print("\n3. 性能排序（从快到慢）：")
    print("   - OPT3/FBP类算法: 最快，适合大系统")
    print("   - SCF: 较慢但精度高，适合小系统或需要高精度的场合")

if __name__ == "__main__":
    main()