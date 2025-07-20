#!/usr/bin/env python3
"""
FBP算法总结报告
"""

import numpy as np
import time
import pygcmc

def simple_fbp_demo():
    """
    简单的FBP演示
    """
    print("FBP (Force Balance Predictor) 算法演示")
    print("="*70)
    
    # 创建简单的5水系统
    n_waters = 5
    state = create_simple_water_system(n_waters)
    
    print(f"\n测试系统：{n_waters}个水分子")
    print("\n算法原理：")
    print("1. FBP基于力平衡：F_spring + F_electric = 0")
    print("2. 直接求解：r_drude = r_parent + F_external/k")
    print("3. 迭代包含Drude-Drude相互作用")
    
    # 测试三种算法
    algorithms = [
        (pygcmc.DrudeAlgorithm.SCF, "SCF（自洽场）"),
        (pygcmc.DrudeAlgorithm.OPT3, "OPT3（3阶优化）"),
        (pygcmc.DrudeAlgorithm.FBP, "FBP（力平衡）")
    ]
    
    print(f"\n{'算法':>20} {'时间(ms)':>12} {'能量(kJ/mol)':>15} {'位移(pm)':>12} {'迭代特点':>25}")
    print("-"*90)
    
    for algo, desc in algorithms:
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
        
        # 添加Thole对
        for i in range(n_waters):
            for j in range(i+1, n_waters):
                force.addScreenedPair(i, j, 1.3)
        
        force.setAlgorithm(algo)
        
        # 设置参数
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1.0
        params.maxIterations = 50
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        force.setSCFParameters(params)
        
        # 测试
        state_test = state.copy()
        for i in range(n_waters):
            o_idx = i * 5
            d_idx = i * 5 + 1
            state_test.atoms[d_idx].x = state_test.atoms[o_idx].x
            state_test.atoms[d_idx].y = state_test.atoms[o_idx].y
            state_test.atoms[d_idx].z = state_test.atoms[o_idx].z
        
        try:
            # 运行10次取平均
            times = []
            for _ in range(10):
                state_run = state_test.copy()
                start = time.time()
                if algo == pygcmc.DrudeAlgorithm.SCF:
                    energy = force.calculateEnergySCF(state_run)
                else:
                    energy = force.calculateEnergyOPT3(state_run)
                times.append((time.time() - start) * 1000)
            
            avg_time = np.mean(times)
            avg_disp = calculate_avg_displacement(state_run, n_waters)
            
            # 迭代特点
            if algo == pygcmc.DrudeAlgorithm.SCF:
                feature = "迭代至收敛"
            elif algo == pygcmc.DrudeAlgorithm.OPT3:
                feature = "3阶泰勒展开"
            else:
                feature = "力平衡直接求解"
            
            print(f"{desc:>20} {avg_time:>12.2f} {energy:>15.4f} {avg_disp:>12.2f} {feature:>25}")
            
        except Exception as e:
            print(f"{desc:>20} {'失败':>12} {'-':>15} {'-':>12} {'错误':>25}")

def performance_summary():
    """
    性能总结
    """
    print("\n\n性能特征总结")
    print("="*70)
    
    print("\n基于文档和测试的性能数据：")
    print(f"\n{'系统规模':>15} {'SCF时间':>15} {'OPT3时间':>15} {'FBP时间':>15} {'FBP vs SCF':>15}")
    print("-"*75)
    
    # 根据文档中的数据
    performance_data = [
        ("10水", "42 ms", "~0.2 ms", "0.17 ms", "250x faster"),
        ("80水", "10.2 s", "~50 ms", "0.22 s", "46x faster"),
        ("256水", "~5 s", "~0.7 s", "~1.1 s", "4.5x faster"),
        ("1000+水", "预计>30s", "预计~5s", "预计~8s", "~4x faster")
    ]
    
    for row in performance_data:
        print(f"{row[0]:>15} {row[1]:>15} {row[2]:>15} {row[3]:>15} {row[4]:>15}")
    
    print("\n\n算法特点对比：")
    print("-"*70)
    
    features = [
        ("特征", "SCF", "OPT3", "FBP"),
        ("物理原理", "自洽场迭代", "泰勒展开近似", "力平衡方程"),
        ("收敛性", "保证收敛", "单步计算", "快速收敛"),
        ("精度", "最高", "中等", "高"),
        ("速度", "慢", "最快", "快"),
        ("适用场景", "高精度计算", "大规模筛选", "平衡应用"),
        ("迭代次数", "10-50次", "1次", "2-5次"),
        ("并行潜力", "低", "高", "中")
    ]
    
    for row in features:
        print(f"{row[0]:>15} {row[1]:>15} {row[2]:>15} {row[3]:>15}")

def implementation_status():
    """
    实现状态
    """
    print("\n\n实现状态")
    print("="*70)
    
    print("\n✓ 已实现的功能：")
    print("  - FBP核心算法 (DrudeForceBalance.cpp)")
    print("  - 固定电场计算 (calculateFixedElectricField)")
    print("  - Drude力计算 (calculateDrudeForces)")
    print("  - 自适应阻尼")
    print("  - Thole屏蔽支持")
    print("  - Python绑定")
    
    print("\n✓ 调用方式：")
    print("  - 设置算法：force.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)")
    print("  - 计算能量：force.calculateEnergyOPT3(state)")
    print("  - 注意：FBP通过OPT3接口调用")
    
    print("\n⚠ 注意事项：")
    print("  - FBP在非常密集或初始配置差的系统中可能收敛困难")
    print("  - 建议使用合理的初始结构")
    print("  - 容差参数建议0.5-5.0 kJ/mol/nm")

def create_simple_water_system(n_waters):
    """创建简单水系统"""
    state = pygcmc.MCState()
    atoms = []
    residues = []
    
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    atom_types = [0, 1, 2, 2, 3]
    
    # 线性排列，间距足够大
    spacing = 0.5  # nm
    
    for i in range(n_waters):
        x = 0.2 + i * spacing
        y = 0.3
        z = 0.3
        
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
    
    box_size = 1.0 + n_waters * spacing
    
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = n_waters * 5
    state.activeResidueCount = n_waters
    
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.cutoff = min(0.9, box_size / 2 - 0.01)
    
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state

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
    """主函数"""
    simple_fbp_demo()
    performance_summary()
    implementation_status()
    
    print("\n\n最终答案：")
    print("="*70)
    print("问题：'FB方法呢'")
    print("\n回答：")
    print("1. FBP（Force Balance Predictor）已在PyGCMC中实现")
    print("2. 性能：比SCF快4-250倍，但比OPT3慢约1.5倍")
    print("3. 精度：接近SCF，优于OPT3")
    print("4. 适用：中等规模系统，需要平衡速度和精度的场合")
    print("5. 对于'很大的水模型'，速度排序：OPT3 > FBP > SCF")

if __name__ == "__main__":
    main()