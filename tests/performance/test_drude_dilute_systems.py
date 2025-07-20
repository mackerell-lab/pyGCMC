#!/usr/bin/env python3
"""
测试稀释系统中的Drude SCF - 避免PBC伪影
使用更大的盒子确保Thole截断距离小于半盒子
"""

import pygcmc
import numpy as np
import time

def create_dilute_water_system(n_waters, density=0.1):
    """
    创建稀释的水系统（密度0.1 g/cm³）
    """
    # 计算盒子大小
    molar_mass_water = 18.01528  # g/mol
    avogadro = 6.02214076e23     # mol⁻¹
    
    total_mass_g = n_waters * molar_mass_water / avogadro
    target_volume_cm3 = total_mass_g / density
    target_volume_nm3 = target_volume_cm3 * 1e21
    box_length = target_volume_nm3 ** (1.0/3.0)
    
    # 创建状态
    state = pygcmc.MCState()
    atoms = []
    residues = []
    
    # SWM4-NDP参数
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]  # O, D, H1, H2, M
    atom_types = [0, 1, 2, 2, 3]
    
    # 在立方格子上放置水分子
    n_per_side = int(np.ceil(n_waters ** (1.0/3.0)))
    spacing = box_length / n_per_side
    
    water_count = 0
    np.random.seed(42)
    
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                if water_count >= n_waters:
                    break
                
                # 基础位置
                x = (i + 0.5) * spacing + (np.random.random() - 0.5) * spacing * 0.1
                y = (j + 0.5) * spacing + (np.random.random() - 0.5) * spacing * 0.1
                z = (k + 0.5) * spacing + (np.random.random() - 0.5) * spacing * 0.1
                
                # 随机旋转
                angles = np.random.random(3) * 2 * np.pi
                
                # 水分子几何
                r_oh = 0.09572
                angle_hoh = 104.52 * np.pi / 180.0
                
                # H1和H2位置（相对于O）
                h1_rel = np.array([r_oh * np.sin(angle_hoh/2), 0, r_oh * np.cos(angle_hoh/2)])
                h2_rel = np.array([-r_oh * np.sin(angle_hoh/2), 0, r_oh * np.cos(angle_hoh/2)])
                
                # 简单旋转
                c, s = np.cos(angles[0]), np.sin(angles[0])
                rot_z = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
                h1_rel = rot_z @ h1_rel
                h2_rel = rot_z @ h2_rel
                
                # M位点
                bisector = -(h1_rel + h2_rel)
                bisector_norm = bisector / np.linalg.norm(bisector)
                m_rel = bisector_norm * 0.024034
                
                # 添加5个原子
                for idx, (dx, dy, dz, charge, atype) in enumerate([
                    (0, 0, 0, charges[0], atom_types[0]),  # O
                    (0, 0, 0, charges[1], atom_types[1]),  # D
                    (h1_rel[0], h1_rel[1], h1_rel[2], charges[2], atom_types[2]),  # H1
                    (h2_rel[0], h2_rel[1], h2_rel[2], charges[3], atom_types[3]),  # H2
                    (m_rel[0], m_rel[1], m_rel[2], charges[4], atom_types[4])  # M
                ]):
                    atom = pygcmc.MCAtom()
                    atom.x = x + dx
                    atom.y = y + dy
                    atom.z = z + dz
                    atom.charge = charge
                    atom.type = atype
                    atoms.append(atom)
                
                # 创建残基
                res = pygcmc.MCResidue()
                res.atomStart = 5 * water_count
                res.atomCount = 5
                res.active = True
                res.type = 0
                residues.append(res)
                
                water_count += 1
                
            if water_count >= n_waters:
                break
        if water_count >= n_waters:
            break
    
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = n_waters
    
    # 设置盒子和截断
    state.info.box = np.array([box_length, box_length, box_length])
    state.info.cutoff = min(1.2, box_length / 2 - 0.01)
    
    # 设置力场参数
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state, box_length, density

def test_dilute_vs_dense(n_waters_list=[4, 8, 16]):
    """
    比较稀释系统和密集系统的Drude SCF行为
    """
    print("比较稀释(0.1 g/cm³)和密集(1.0 g/cm³)系统中的Drude SCF")
    print("="*80)
    
    densities = [0.1, 1.0]
    
    for n_waters in n_waters_list:
        print(f"\n{n_waters} 水分子系统:")
        print("-"*60)
        
        results = {}
        
        for density in densities:
            # 创建系统
            state, box_length, actual_density = create_dilute_water_system(n_waters, density)
            
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
            
            # 添加Thole对
            n_thole_pairs = 0
            thole_cutoff = 0.8  # 标准Thole截断
            thole_param = 1.3
            
            # 检查是否需要调整截断
            half_box = box_length / 2.0
            effective_cutoff = min(thole_cutoff, half_box - 0.01)
            
            for i in range(n_waters):
                o1_idx = i * 5
                for j in range(i+1, n_waters):
                    o2_idx = j * 5
                    
                    # 计算O-O距离
                    dx = state.atoms[o1_idx].x - state.atoms[o2_idx].x
                    dy = state.atoms[o1_idx].y - state.atoms[o2_idx].y
                    dz = state.atoms[o1_idx].z - state.atoms[o2_idx].z
                    
                    # PBC
                    dx -= box_length * round(dx / box_length)
                    dy -= box_length * round(dy / box_length)
                    dz -= box_length * round(dz / box_length)
                    
                    dist = np.sqrt(dx*dx + dy*dy + dz*dz)
                    
                    if dist < effective_cutoff:
                        force.addScreenedPair(i, j, thole_param)
                        n_thole_pairs += 1
            
            # 设置SCF参数
            params = pygcmc.DrudeSCFParams()
            params.tolerance = 10.0
            params.maxIterations = 200
            params.dampingFactor = 0.5
            params.maxDrudeDistance = 0.02
            force.setSCFParameters(params)
            force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
            
            # 运行SCF
            start_time = time.time()
            try:
                energy = force.calculateEnergySCF(state)
                elapsed_time = (time.time() - start_time) * 1000
                converged = True
                
                # 分析Drude位移
                displacements = []
                for i in range(n_waters):
                    o_idx = i * 5
                    d_idx = i * 5 + 1
                    
                    dx = state.atoms[d_idx].x - state.atoms[o_idx].x
                    dy = state.atoms[d_idx].y - state.atoms[o_idx].y
                    dz = state.atoms[d_idx].z - state.atoms[o_idx].z
                    
                    disp = np.sqrt(dx*dx + dy*dy + dz*dz) * 1000
                    displacements.append(disp)
                
                avg_disp = np.mean(displacements)
                
            except Exception as e:
                elapsed_time = (time.time() - start_time) * 1000
                converged = False
                energy = float('nan')
                avg_disp = float('nan')
            
            results[density] = {
                'box_length': box_length,
                'half_box': half_box,
                'thole_cutoff': effective_cutoff,
                'n_thole_pairs': n_thole_pairs,
                'converged': converged,
                'time_ms': elapsed_time,
                'energy_per_water': energy / n_waters if converged else float('nan'),
                'avg_displacement': avg_disp
            }
            
            print(f"\n  密度 {density} g/cm³:")
            print(f"    盒子长度: {box_length:.3f} nm")
            print(f"    半盒子: {half_box:.3f} nm")
            print(f"    Thole截断: {effective_cutoff:.3f} nm (标准: 0.8 nm)")
            print(f"    Thole对数: {n_thole_pairs}")
            print(f"    收敛: {'是' if converged else '否'}")
            print(f"    时间: {elapsed_time:.1f} ms")
            if converged:
                print(f"    能量/水: {energy/n_waters:.1f} kJ/mol")
                print(f"    平均Drude位移: {avg_disp:.1f} pm")
        
        # 比较结果
        print(f"\n  比较:")
        if results[0.1]['converged'] and results[1.0]['converged']:
            energy_diff = abs(results[0.1]['energy_per_water'] - results[1.0]['energy_per_water'])
            print(f"    能量差异: {energy_diff:.1f} kJ/mol/水")
            print(f"    稀释系统Thole对: {results[0.1]['n_thole_pairs']}")
            print(f"    密集系统Thole对: {results[1.0]['n_thole_pairs']}")

def main():
    """
    主测试函数
    """
    print("Thole截断距离对Drude SCF的影响分析")
    print("="*80)
    
    # 测试不同密度
    test_dilute_vs_dense([4, 8, 16, 32])
    
    print("\n\n结论:")
    print("1. 稀释系统（0.1 g/cm³）有更大的盒子，可以使用完整的0.8 nm Thole截断")
    print("2. 密集系统（1.0 g/cm³）的小盒子需要减小Thole截断以满足PBC")
    print("3. 两种密度下的能量应该相似，因为分子间相互作用主要由Thole截断决定")
    print("4. 建议：对于小测试系统，使用稀释条件或确保Thole截断 < 半盒子")

if __name__ == "__main__":
    main()