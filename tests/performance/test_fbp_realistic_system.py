#!/usr/bin/env python3
"""
使用真实的水系统测试FBP性能
从预平衡的配置开始，避免收敛问题
"""

import pygcmc
import numpy as np
import time

def create_equilibrated_water_box(n_waters=128):
    """
    创建一个预平衡的水分子系统
    使用更大的间距和合理的初始配置
    """
    print(f"\n创建{n_waters}个水分子系统（预平衡配置）...")
    
    # 使用实际的水密度计算盒子大小
    # 水的密度约1 g/cm³，分子量18 g/mol
    # 1 mol水占体积 = 18 cm³ = 18e21 nm³
    # 1个水分子占体积 = 18e21 / 6.022e23 = 29.9 nm³
    volume_per_water = 30.0  # nm³
    total_volume = n_waters * volume_per_water
    box_length = np.cbrt(total_volume)
    
    # 但为了避免初始重叠，使用更大的盒子
    box_length *= 1.5  # 给更多空间
    
    print(f"盒子尺寸: {box_length:.2f} nm")
    
    atoms = []
    residues = []
    
    # 在立方体中均匀分布水分子
    n_per_side = int(np.ceil(n_waters**(1/3)))
    spacing = box_length / n_per_side
    
    water_count = 0
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                if water_count >= n_waters:
                    break
                
                # 中心位置加小扰动
                x_base = (i + 0.5 + 0.1*np.random.randn()) * spacing
                y_base = (j + 0.5 + 0.1*np.random.randn()) * spacing
                z_base = (k + 0.5 + 0.1*np.random.randn()) * spacing
                
                # 随机旋转水分子
                angle = np.random.rand() * 2 * np.pi
                
                # SWM4-NDP水模型
                positions = [
                    (0.0, 0.0, 0.0, 1.71636, 0),   # O
                    (0.0, 0.0, 0.0, -1.71636, 1),  # D (初始与O重合)
                    (0.09572, 0.0, 0.0, 0.55733, 2),  # H1
                    (-0.04786, 0.08288, 0.0, 0.55733, 2),  # H2
                    (0.0, -0.024034, 0.0, -1.11466, 3)  # M
                ]
                
                for dx, dy, dz, charge, typ in positions:
                    # 旋转
                    x_rot = dx * np.cos(angle) - dy * np.sin(angle)
                    y_rot = dx * np.sin(angle) + dy * np.cos(angle)
                    
                    atom = pygcmc.MCAtom()
                    atom.x = x_base + x_rot
                    atom.y = y_base + y_rot
                    atom.z = z_base + dz
                    atom.charge = charge
                    atom.type = typ
                    atoms.append(atom)
                
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
    
    state = pygcmc.MCState()
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = n_waters
    
    # 周期性边界
    state.info.box = np.array([box_length, box_length, box_length])
    state.info.cutoff = min(1.2, box_length/2 - 0.1)  # 标准截断1.2 nm
    
    # SWM4-NDP力场
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    print(f"系统: {n_waters}水, {len(atoms)}原子, 截断{state.info.cutoff:.2f} nm")
    
    return state

def benchmark_algorithms(n_waters=128):
    """基准测试各算法"""
    print(f"\n{'='*80}")
    print(f"基准测试: {n_waters}水分子系统")
    print(f"{'='*80}")
    
    # 创建系统
    state = create_equilibrated_water_box(n_waters)
    
    # 测试配置
    tests = [
        ("SCF-0.001", "SCF", 0.001, 2000),    # 超高精度
        ("SCF-0.01", "SCF", 0.01, 1000),      # 高精度
        ("SCF-0.1", "SCF", 0.1, 500),         # 中精度
        ("SCF-1.0", "SCF", 1.0, 500),         # 低精度
        ("FBP-0.1", "FBP", 0.1, 200),         # FBP高精度
        ("FBP-1.0", "FBP", 1.0, 200),         # FBP标准
        ("FBP-10.0", "FBP", 10.0, 200),       # FBP快速
    ]
    
    results = []
    ref_energy = None
    
    for test_name, algo, tol, max_iter in tests:
        print(f"\n测试 {test_name}...")
        
        force = pygcmc.DrudeForce()
        
        # SWM4-NDP参数
        charge = -1.71636
        k_spring = 418400.0
        polarizability = 1.71636**2 * 138.935456 / k_spring
        
        # 添加Drude
        for i in range(n_waters):
            force.addParticle(
                drudeIndex=5*i+1, parentIndex=5*i,
                aniso1Index=-1, aniso2Index=-1,
                aniso3Index=-1, aniso4Index=-1,
                charge=charge, polarizability=polarizability,
                aniso12=1.0, aniso34=1.0
            )
        
        # Thole屏蔽
        for i in range(n_waters):
            for j in range(i+1, n_waters):
                force.addScreenedPair(i, j, 1.3)
        
        # 参数设置
        params = pygcmc.DrudeSCFParams()
        params.tolerance = tol
        params.maxIterations = max_iter
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        force.setSCFParameters(params)
        
        if algo == "SCF":
            force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        else:
            force.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
        
        # 重置Drude
        state_copy = state.copy()
        for i in range(n_waters):
            state_copy.atoms[5*i+1].x = state_copy.atoms[5*i].x
            state_copy.atoms[5*i+1].y = state_copy.atoms[5*i].y
            state_copy.atoms[5*i+1].z = state_copy.atoms[5*i].z
        
        # 计时
        start = time.time()
        try:
            energy = force.calculateEnergySCF(state_copy)
            elapsed = (time.time() - start) * 1000
            
            # 第一个成功的能量作为参考
            if ref_energy is None:
                ref_energy = energy
            
            # 计算位移
            disps = []
            for i in range(n_waters):
                dx = state_copy.atoms[5*i+1].x - state_copy.atoms[5*i].x
                dy = state_copy.atoms[5*i+1].y - state_copy.atoms[5*i].y
                dz = state_copy.atoms[5*i+1].z - state_copy.atoms[5*i].z
                disps.append(np.sqrt(dx*dx + dy*dy + dz*dz))
            
            avg_disp = np.mean(disps) * 1000  # pm
            max_disp = np.max(disps) * 1000
            
            results.append({
                'name': test_name,
                'energy': energy,
                'time': elapsed,
                'avg_disp': avg_disp,
                'max_disp': max_disp,
                'converged': True
            })
            
            print(f"  能量: {energy:.2f} kJ/mol")
            print(f"  时间: {elapsed:.1f} ms")
            print(f"  平均位移: {avg_disp:.2f} pm")
            
        except Exception as e:
            print(f"  失败: {e}")
            results.append({
                'name': test_name,
                'energy': float('nan'),
                'time': float('nan'),
                'converged': False
            })
    
    # 汇总结果
    print(f"\n\n{'='*100}")
    print(f"结果汇总 ({n_waters}水分子)")
    print(f"{'='*100}")
    print(f"{'算法':<12} {'能量(kJ/mol)':<15} {'与参考差':<12} {'误差(%)':<10} {'时间(ms)':<10} {'速度比':<8} {'平均位移(pm)':<12}")
    print(f"{'-'*100}")
    
    # 找到第一个收敛的SCF作为时间基准
    ref_time = None
    for r in results:
        if r['converged'] and r['name'].startswith('SCF'):
            ref_time = r['time']
            break
    
    for r in results:
        if r['converged']:
            error = abs(r['energy'] - ref_energy)
            error_pct = error / abs(ref_energy) * 100 if abs(ref_energy) > 0.1 else 0
            speedup = ref_time / r['time'] if ref_time else 1.0
            
            print(f"{r['name']:<12} {r['energy']:<15.2f} {error:<12.2f} {error_pct:<10.3f} {r['time']:<10.1f} {speedup:<8.1f}x {r['avg_disp']:<12.2f}")
        else:
            print(f"{r['name']:<12} {'未收敛':<15} {'-':<12} {'-':<10} {'-':<10} {'-':<8} {'-':<12}")
    
    # 专门比较FBP和SCF
    print(f"\n\nFBP vs SCF 详细对比:")
    print(f"{'-'*60}")
    
    # 找到对应的结果
    scf_01 = next((r for r in results if r['name'] == 'SCF-0.1'), None)
    fbp_10 = next((r for r in results if r['name'] == 'FBP-1.0'), None)
    
    if scf_01 and fbp_10 and scf_01['converged'] and fbp_10['converged']:
        speedup = scf_01['time'] / fbp_10['time']
        error = abs(fbp_10['energy'] - scf_01['energy']) / abs(scf_01['energy']) * 100
        
        print(f"SCF (tol=0.1) vs FBP (tol=1.0):")
        print(f"  速度提升: {speedup:.1f}x")
        print(f"  能量误差: {error:.2f}%")
        print(f"  FBP时间: {fbp_10['time']:.1f} ms")
        print(f"  SCF时间: {scf_01['time']:.1f} ms")

def test_different_sizes():
    """测试不同大小系统"""
    print("\n\n不同系统大小的性能测试")
    print("="*80)
    
    sizes = [32, 64, 128]
    
    for n in sizes:
        benchmark_algorithms(n)

if __name__ == "__main__":
    # 先测试128水系统
    benchmark_algorithms(128)
    
    # 再测试扩展性
    test_different_sizes()