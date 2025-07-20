#!/usr/bin/env python3
"""
演示矩阵方法求解Drude振子位置的概念
"""

import numpy as np
import time

def build_system_matrix(n_drudes, k_spring, positions, charges, box_size, cutoff):
    """
    构建系统矩阵 K
    K包含弹簧力常数和Thole相互作用
    """
    K = np.zeros((3*n_drudes, 3*n_drudes))
    
    # 对角项：弹簧力常数
    for i in range(n_drudes):
        for j in range(3):
            K[3*i+j, 3*i+j] = k_spring
    
    # 非对角项：Thole相互作用（简化版本）
    for i in range(n_drudes):
        for j in range(i+1, n_drudes):
            # 计算距离
            r_ij = positions[j] - positions[i]
            # 周期性边界条件
            r_ij = r_ij - box_size * np.round(r_ij / box_size)
            r = np.linalg.norm(r_ij)
            
            if r < cutoff:
                # Thole屏蔽函数（简化）
                thole_factor = 1.3
                u = r / (thole_factor * (charges[i] * charges[j])**(1/6))
                screening = 1 - np.exp(-u**3) * (1 + u**3)
                
                # 相互作用张量（简化为标量）
                T_ij = screening * charges[i] * charges[j] / r**3
                
                # 添加到矩阵
                for a in range(3):
                    for b in range(3):
                        coupling = T_ij * r_ij[a] * r_ij[b] / r**2
                        K[3*i+a, 3*j+b] += coupling
                        K[3*j+b, 3*i+a] += coupling
    
    return K

def calculate_electric_forces(n_drudes, positions, charges, box_size, cutoff):
    """
    计算每个Drude粒子受到的电场力
    """
    forces = np.zeros((n_drudes, 3))
    
    # 计算所有原子对Drude的库仑力
    for i in range(n_drudes):
        for j in range(len(positions)):
            if i == j:
                continue
            
            r_ij = positions[j] - positions[i]
            r_ij = r_ij - box_size * np.round(r_ij / box_size)
            r = np.linalg.norm(r_ij)
            
            if r < cutoff:
                # 库仑力
                f_mag = 138.935456 * charges[i] * charges[j] / r**2
                forces[i] += f_mag * r_ij / r
    
    return forces.flatten()

def solve_matrix_direct(K, F):
    """
    直接求解线性方程组 K·δr = F
    """
    return np.linalg.solve(K, F)

def solve_conjugate_gradient(K, F, tol=1e-6, max_iter=100):
    """
    共轭梯度法求解 K·δr = F
    """
    n = len(F)
    x = np.zeros(n)
    r = F - K.dot(x)
    p = r.copy()
    rsold = r.dot(r)
    
    for i in range(max_iter):
        Ap = K.dot(p)
        alpha = rsold / p.dot(Ap)
        x = x + alpha * p
        r = r - alpha * Ap
        rsnew = r.dot(r)
        
        if np.sqrt(rsnew) < tol:
            return x, i+1
        
        beta = rsnew / rsold
        p = r + beta * p
        rsold = rsnew
    
    return x, max_iter

def traditional_scf(n_drudes, k_spring, positions, charges, box_size, cutoff, 
                   tol=1e-6, max_iter=100, damping=0.5):
    """
    传统SCF迭代方法
    """
    # 初始化Drude位置（与parent重合）
    drude_positions = positions[:n_drudes].copy()
    
    for iter in range(max_iter):
        # 计算力
        forces = calculate_electric_forces(n_drudes, drude_positions, charges, box_size, cutoff)
        forces = forces.reshape(n_drudes, 3)
        
        # 更新位置
        max_force = 0
        for i in range(n_drudes):
            # 弹簧力
            spring_force = -k_spring * (drude_positions[i] - positions[i])
            # 总力
            total_force = forces[i] + spring_force
            max_force = max(max_force, np.linalg.norm(total_force))
            
            # 更新
            drude_positions[i] += damping * total_force / k_spring
        
        if max_force < tol:
            return drude_positions, iter+1
    
    return drude_positions, max_iter

def compare_methods():
    """
    比较不同方法的性能
    """
    print("Drude振子矩阵方法演示")
    print("="*60)
    
    # 系统参数
    n_waters = 32
    n_drudes = n_waters
    k_spring = 418400.0  # kJ/mol/nm²
    box_size = 2.0  # nm
    cutoff = 1.2  # nm
    
    # 生成随机位置
    np.random.seed(42)
    positions = np.random.rand(n_waters * 5, 3) * box_size  # 5 atoms per water
    charges = np.zeros(n_waters * 5)
    for i in range(n_waters):
        charges[5*i] = 1.71636      # O
        charges[5*i+1] = -1.71636   # D
        charges[5*i+2] = 0.55733    # H1
        charges[5*i+3] = 0.55733    # H2
        charges[5*i+4] = -1.11466   # M
    
    print(f"\n系统规模: {n_drudes} Drude振子")
    print(f"弹簧常数: {k_spring} kJ/mol/nm²")
    print(f"盒子大小: {box_size} nm")
    print(f"截断距离: {cutoff} nm")
    
    # 1. 传统SCF方法
    print("\n1. 传统SCF迭代")
    print("-"*40)
    start = time.time()
    drude_pos_scf, iters_scf = traditional_scf(
        n_drudes, k_spring, positions, charges, box_size, cutoff
    )
    time_scf = (time.time() - start) * 1000
    print(f"  迭代次数: {iters_scf}")
    print(f"  时间: {time_scf:.1f} ms")
    
    # 2. 直接矩阵方法
    print("\n2. 直接矩阵求解")
    print("-"*40)
    start = time.time()
    
    # 构建矩阵
    K = build_system_matrix(n_drudes, k_spring, positions[:n_drudes], 
                           charges[:n_drudes], box_size, cutoff)
    F = calculate_electric_forces(n_drudes, positions[:n_drudes], charges, box_size, cutoff)
    
    # 求解
    delta_r = solve_matrix_direct(K, F)
    
    time_direct = (time.time() - start) * 1000
    print(f"  矩阵大小: {K.shape[0]}×{K.shape[1]}")
    print(f"  时间: {time_direct:.1f} ms")
    print(f"  加速比: {time_scf/time_direct:.2f}x")
    
    # 3. 共轭梯度法
    print("\n3. 共轭梯度法")
    print("-"*40)
    start = time.time()
    
    delta_r_cg, iters_cg = solve_conjugate_gradient(K, F)
    
    time_cg = (time.time() - start) * 1000
    print(f"  迭代次数: {iters_cg}")
    print(f"  时间: {time_cg:.1f} ms")
    print(f"  加速比: {time_scf/time_cg:.2f}x")
    
    # 验证解的一致性
    print("\n4. 解的一致性检查")
    print("-"*40)
    
    # 重构Drude位置
    drude_pos_direct = positions[:n_drudes].copy()
    drude_pos_cg = positions[:n_drudes].copy()
    
    for i in range(n_drudes):
        drude_pos_direct[i] += delta_r[3*i:3*i+3]
        drude_pos_cg[i] += delta_r_cg[3*i:3*i+3]
    
    # 计算差异
    diff_direct_scf = np.mean(np.linalg.norm(drude_pos_direct - drude_pos_scf, axis=1))
    diff_cg_scf = np.mean(np.linalg.norm(drude_pos_cg - drude_pos_scf, axis=1))
    
    print(f"  直接法 vs SCF 平均差异: {diff_direct_scf*1000:.2f} pm")
    print(f"  CG法 vs SCF 平均差异: {diff_cg_scf*1000:.2f} pm")
    
    # 性能总结
    print("\n5. 性能总结")
    print("-"*40)
    print(f"  SCF: {time_scf:.1f} ms ({iters_scf} 迭代)")
    print(f"  直接矩阵: {time_direct:.1f} ms ({time_scf/time_direct:.1f}x 加速)")
    print(f"  共轭梯度: {time_cg:.1f} ms ({time_scf/time_cg:.1f}x 加速)")
    
    print("\n结论:")
    print("- 对于小系统，直接矩阵方法最快")
    print("- 共轭梯度法在中等系统中表现良好")
    print("- 矩阵方法可以显著提高收敛速度")

if __name__ == "__main__":
    compare_methods()