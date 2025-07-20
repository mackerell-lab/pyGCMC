#!/usr/bin/env python3
"""
测试修正Thole力计算后的SCF收敛性
"""

import pygcmc
import numpy as np
import time

def test_scf_convergence_after_fix():
    """
    测试修正Thole力后的SCF收敛性
    """
    print("测试修正Thole力计算后的SCF收敛性")
    print("="*70)
    
    # 使用之前生成的密度1.0的水系统
    system_sizes = [2, 4, 8, 16, 32]
    
    print(f"\n{'系统':<10} {'Thole对':<10} {'收敛?':<8} {'迭代次数':<10} {'时间(ms)':<10} {'能量/水':<12}")
    print("-"*70)
    
    for n_waters in system_sizes:
        try:
            # 加载之前生成的系统
            import pickle
            pickle_file = f'../tests/performance/water_density_1.0/water_{n_waters}.pkl'
            with open(pickle_file, 'rb') as f:
                data = pickle.load(f)
            
            # 创建状态
            state = pygcmc.MCState()
            atoms = []
            residues = []
            
            # SWM4-NDP参数
            charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
            atom_types = [0, 1, 2, 2, 3]
            
            positions = data['positions']
            box_length = data['box_length']
            
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
            state.info.box = np.array([box_length, box_length, box_length])
            state.info.cutoff = min(1.2, box_length / 2 - 0.01)
            
            # 设置力场参数
            state.forcefield.numTotalTypes = 4
            state.forcefield.numMovementTypes = 4
            state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
            state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
            
            # 创建DrudeForce
            force = pygcmc.DrudeForce()
            
            # 添加Drude粒子
            drude_charge = -1.71636
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
            
            # 添加Thole对（动态截断）
            n_thole_pairs = 0
            thole_cutoff = min(0.8, box_length / 2.0 - 0.01)
            
            for i in range(n_waters):
                o1_idx = i * 5
                for j in range(i+1, n_waters):
                    o2_idx = j * 5
                    
                    dx = state.atoms[o1_idx].x - state.atoms[o2_idx].x
                    dy = state.atoms[o1_idx].y - state.atoms[o2_idx].y
                    dz = state.atoms[o1_idx].z - state.atoms[o2_idx].z
                    
                    # PBC
                    dx -= box_length * round(dx / box_length)
                    dy -= box_length * round(dy / box_length)
                    dz -= box_length * round(dz / box_length)
                    
                    dist = np.sqrt(dx*dx + dy*dy + dz*dz)
                    
                    if dist < thole_cutoff:
                        force.addScreenedPair(i, j, 1.3)
                        n_thole_pairs += 1
            
            # 设置SCF参数
            params = pygcmc.DrudeSCFParams()
            params.tolerance = 10.0
            params.maxIterations = 100
            params.dampingFactor = 0.5
            params.maxDrudeDistance = 0.02
            force.setSCFParameters(params)
            force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
            
            # 记录迭代次数（这需要在C++中实现）
            start_time = time.time()
            
            try:
                energy = force.calculateEnergySCF(state)
                elapsed_time = (time.time() - start_time) * 1000
                converged = True
                iterations = "收敛"  # 需要C++返回实际迭代次数
            except Exception as e:
                elapsed_time = (time.time() - start_time) * 1000
                converged = False
                energy = float('nan')
                iterations = "未收敛"
            
            print(f"{n_waters:<10} {n_thole_pairs:<10} {'是' if converged else '否':<8} "
                  f"{iterations:<10} {elapsed_time:<10.1f} {energy/n_waters if converged else 'NaN':<12.1f}")
            
        except Exception as e:
            print(f"{n_waters:<10} 错误: {str(e)}")

def main():
    """
    主函数
    """
    test_scf_convergence_after_fix()
    
    print("\n\n总结:")
    print("1. 修正Thole力计算公式后，SCF收敛性应该显著改善")
    print("2. 小系统由于PBC限制仍可能有问题")
    print("3. 如果仍有收敛问题，可能需要进一步调整阻尼因子或检查其他力的计算")

if __name__ == "__main__":
    main()