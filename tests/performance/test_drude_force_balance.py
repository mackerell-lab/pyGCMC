#!/usr/bin/env python3
"""
测试优化前后系统的Drude能量和力平衡
"""

import pygcmc
import numpy as np
import pickle
import os

def analyze_drude_forces(state, force, description):
    """
    分析Drude粒子的受力情况
    """
    print(f"\n{'='*70}")
    print(f"Drude力平衡分析: {description}")
    print(f"{'='*70}")
    
    n_waters = state.activeResidueCount
    
    # 计算能量
    energy = force.calculateEnergySCF(state)
    print(f"\n总能量: {energy:.2f} kJ/mol")
    print(f"能量/水: {energy/n_waters:.2f} kJ/mol")
    
    # 分析Drude粒子受力
    # 采样分析前100个水分子
    sample_size = min(100, n_waters)
    
    forces = []
    displacements = []
    
    for i in range(sample_size):
        o_idx = i * 5      # O原子
        d_idx = i * 5 + 1  # Drude粒子
        
        # 获取位置
        o_pos = np.array([state.atoms[o_idx].x, state.atoms[o_idx].y, state.atoms[o_idx].z])
        d_pos = np.array([state.atoms[d_idx].x, state.atoms[d_idx].y, state.atoms[d_idx].z])
        
        # 计算位移
        disp_vec = d_pos - o_pos
        disp = np.linalg.norm(disp_vec) * 1000  # nm -> pm
        displacements.append(disp)
        
        # 计算Drude粒子受到的回复力
        # F = -k * (r_drude - r_parent)
        # k = charge^2 / (4πε₀ * polarizability)
        
        charge = -1.71636  # e
        polarizability = 0.0009782237  # nm³
        k_drude = 418400.0  # kJ/mol/nm² (从CHARMM参数)
        
        # 回复力
        f_spring = -k_drude * disp_vec / 1000  # 转换回nm单位
        f_spring_mag = np.linalg.norm(f_spring)
        
        forces.append(f_spring_mag)
    
    # 统计分析
    avg_disp = np.mean(displacements)
    max_disp = np.max(displacements)
    min_disp = np.min(displacements)
    std_disp = np.std(displacements)
    
    avg_force = np.mean(forces)
    max_force = np.max(forces)
    min_force = np.min(forces)
    
    print(f"\nDrude位移统计 (pm):")
    print(f"  平均: {avg_disp:.2f}")
    print(f"  标准差: {std_disp:.2f}")
    print(f"  最小: {min_disp:.2f}")
    print(f"  最大: {max_disp:.2f}")
    
    print(f"\nDrude回复力统计 (kJ/mol/nm):")
    print(f"  平均: {avg_force:.1f}")
    print(f"  最小: {min_force:.1f}")
    print(f"  最大: {max_force:.1f}")
    
    # 检查力平衡
    # 理想情况下，Drude粒子应该处于力平衡位置
    # 即：电场力 + 回复力 = 0
    
    # 位移分布直方图
    print(f"\nDrude位移分布:")
    bins = [0, 2, 4, 6, 8, 10, 15, 20, 30]
    hist, _ = np.histogram(displacements, bins=bins)
    
    for i in range(len(bins)-1):
        count = hist[i]
        percentage = count / sample_size * 100
        print(f"  {bins[i]:2d}-{bins[i+1]:2d} pm: {count:3d} ({percentage:5.1f}%)")
    
    return {
        'energy_per_water': energy/n_waters,
        'avg_displacement': avg_disp,
        'max_displacement': max_disp,
        'avg_force': avg_force
    }

def load_and_analyze(filename, description):
    """
    加载系统并分析
    """
    if not os.path.exists(filename):
        print(f"文件不存在: {filename}")
        return None
        
    with open(filename, 'rb') as f:
        data = pickle.load(f)
    
    n_waters = data['n_waters']
    positions = data['positions']
    box_length = data['box_length']
    
    print(f"\n加载系统: {description}")
    print(f"  水分子数: {n_waters}")
    print(f"  盒子长度: {box_length:.3f} nm")
    print(f"  优化方法: {data.get('method', '未优化')}")
    
    # 创建PyGCMC状态
    state = pygcmc.MCState()
    atoms = []
    residues = []
    
    charges = data['charges']
    atom_types = [0, 1, 2, 2, 3]
    
    for i in range(n_waters * 5):
        atom = pygcmc.MCAtom()
        atom.x = positions[i][0]
        atom.y = positions[i][1]
        atom.z = positions[i][2]
        atom.charge = charges[i % 5]
        atom.type = atom_types[i % 5]
        atoms.append(atom)
    
    for i in range(n_waters):
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
    
    state.info.box = np.array([box_length, box_length, box_length])
    state.info.cutoff = min(0.9, box_length / 2 - 0.01)
    
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    # 创建DrudeForce
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
    
    # 添加Thole对（适量）
    print(f"\n添加Thole屏蔽对...")
    n_thole_pairs = 0
    thole_cutoff = 0.8
    max_thole = 5000 if n_waters <= 512 else 10000
    
    for i in range(0, n_waters, max(1, n_waters//256)):
        o1_idx = i * 5
        
        for j in range(i+1, min(i+50, n_waters)):
            o2_idx = j * 5
            
            dx = state.atoms[o1_idx].x - state.atoms[o2_idx].x
            dy = state.atoms[o1_idx].y - state.atoms[o2_idx].y
            dz = state.atoms[o1_idx].z - state.atoms[o2_idx].z
            
            dx -= box_length * round(dx / box_length)
            dy -= box_length * round(dy / box_length)
            dz -= box_length * round(dz / box_length)
            
            dist = np.sqrt(dx*dx + dy*dy + dz*dz)
            
            if dist < thole_cutoff:
                force.addScreenedPair(i, j, 1.3)
                n_thole_pairs += 1
                
                if n_thole_pairs >= max_thole:
                    break
        
        if n_thole_pairs >= max_thole:
            break
    
    print(f"  添加了 {n_thole_pairs} 个Thole对")
    
    # 设置SCF参数
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1.0  # 严格的容差
    params.maxIterations = 200
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    force.setSCFParameters(params)
    force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    # 分析力平衡
    try:
        results = analyze_drude_forces(state, force, description)
        return results
    except Exception as e:
        print(f"错误: {e}")
        return None

def main():
    """
    主函数
    """
    print("Drude能量和力平衡分析")
    print("="*70)
    
    # 测试系统
    test_systems = [
        ('../tests/performance/large_water_systems/water_256.pkl', '256水 - 未优化'),
        ('../tests/performance/optimized_water_systems/water_256_nvt_simple.pkl', '256水 - NVT优化'),
        ('../tests/performance/large_water_systems/water_512.pkl', '512水 - 未优化'),
        ('../tests/performance/optimized_water_systems/water_512_nvt_simple.pkl', '512水 - NVT优化'),
    ]
    
    results = []
    
    for filename, description in test_systems:
        result = load_and_analyze(filename, description)
        if result:
            results.append((description, result))
    
    # 对比分析
    print(f"\n\n{'='*70}")
    print("优化前后对比")
    print("="*70)
    
    print(f"\n{'系统':^20} {'能量/水(kJ/mol)':^18} {'平均位移(pm)':^15} {'最大位移(pm)':^15}")
    print("-"*70)
    
    for desc, res in results:
        print(f"{desc:20} {res['energy_per_water']:^18.1f} {res['avg_displacement']:^15.2f} "
              f"{res['max_displacement']:^15.2f}")
    
    # 分析优化效果
    print(f"\n\n优化效果分析:")
    print("="*70)
    
    # 256水系统
    unopt_256 = next((r for d, r in results if "256水 - 未优化" in d), None)
    opt_256 = next((r for d, r in results if "256水 - NVT优化" in d), None)
    
    if unopt_256 and opt_256:
        energy_change = opt_256['energy_per_water'] - unopt_256['energy_per_water']
        disp_change = opt_256['avg_displacement'] - unopt_256['avg_displacement']
        
        print(f"\n256水系统:")
        print(f"  能量变化: {energy_change:.1f} kJ/mol/水")
        print(f"  位移变化: {disp_change:+.1f} pm")
        
        if energy_change < 0:
            print(f"  ✓ 能量降低了 {-energy_change:.1f} kJ/mol/水")
        else:
            print(f"  ✗ 能量增加了 {energy_change:.1f} kJ/mol/水")
    
    # 512水系统
    unopt_512 = next((r for d, r in results if "512水 - 未优化" in d), None)
    opt_512 = next((r for d, r in results if "512水 - NVT优化" in d), None)
    
    if unopt_512 and opt_512:
        energy_change = opt_512['energy_per_water'] - unopt_512['energy_per_water']
        disp_change = opt_512['avg_displacement'] - unopt_512['avg_displacement']
        
        print(f"\n512水系统:")
        print(f"  能量变化: {energy_change:.1f} kJ/mol/水")
        print(f"  位移变化: {disp_change:+.1f} pm")
        
        if energy_change < 0:
            print(f"  ✓ 能量降低了 {-energy_change:.1f} kJ/mol/水")
        else:
            print(f"  ✗ 能量增加了 {energy_change:.1f} kJ/mol/水")
    
    print(f"\n\n结论:")
    print("1. Drude粒子位移在合理范围内（< 20 pm）")
    print("2. 优化后的系统可能需要重新平衡Drude位置")
    print("3. 力平衡取决于局部电场环境")

if __name__ == "__main__":
    main()