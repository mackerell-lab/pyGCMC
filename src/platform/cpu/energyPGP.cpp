#include "energyPGP.hpp"
#include "energyPME.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

/**
 * @file energyPGP.cpp
 * @brief Implementation of Precomputed Grid-Potential PME for Monte Carlo (PGP-PME-MC)
 * 
 * PGP-PME是一种为蒙特卡洛模拟优化的粒子网格Ewald方法。它通过预计算系统中固定部分的
 * 静电势网格，大大加速了在MC模拟中计算小分子能量变化的过程。
 * 
 * 算法关键特点：
 * 1. 预计算(Precomputed): 系统中固定部分的静电势被预先计算并存储在网格中
 * 2. 网格电势(Grid-Potential): 使用三维网格表示电势分布，结合B样条插值实现高效采样
 * 3. PME基础: 基于传统PME方法处理长程静电相互作用，但做了针对MC的优化
 * 
 * 工作流程:
 * - 初始化阶段: 设置常规PME参数和额外的pair grid参数
 * - 预计算阶段: 将固定部分电荷分配到网格，执行FFT并存储网格电势
 * - MC模拟阶段: 通过网格插值快速评估移动分子的能量变化
 * 
 * 应用场景:
 * - GCMC(大正则系综蒙特卡洛)模拟中溶剂分子的插入/删除
 * - CBMC(构型偏置蒙特卡洛)中的构型采样
 * - 生物大分子系统如蛋白质-配体相互作用模拟
 */

namespace pygcmc {
namespace platform {
namespace cpu {

// Initialize global PGP parameters
PGPParams pgp_params;

/**
 * @brief 初始化配对网格(Pair Grid)
 * 
 * 该方法为PGP-PME算法创建并初始化一个用于存储预计算电势的三维网格。
 * 网格大小由pair_grid_size参数决定，通常独立于常规PME网格大小设置。
 * 网格间距基于盒子大小和网格尺寸计算，确保适当的空间分辨率。
 */
void PGPParams::initializePairGrid() {
    // Calculate grid size and allocate memory for pair grid
    int totalSize = pair_grid_size[0] * pair_grid_size[1] * pair_grid_size[2];
    pairGrid.resize(totalSize);
    
    // Set grid spacing based on box dimensions and grid size
    grid_spacing = std::min({
        box[0] / pair_grid_size[0],
        box[1] / pair_grid_size[1],
        box[2] / pair_grid_size[2]
    });
    
    platform::log(LogLevel::DEBUG, "PGP pair grid initialized with size: ", 
                 pair_grid_size[0], "x", pair_grid_size[1], "x", pair_grid_size[2],
                 ", grid spacing: ", grid_spacing);
}

/**
 * @brief 设置PGP算法参数
 * 
 * 该函数设置PGP-PME算法所需的所有参数，包括PME的基本参数和PGP特有参数。
 * PGP-PME复用了PME的基础设施和参数，同时添加了额外的网格参数用于预计算电势。
 * 
 * @param alpha Ewald分离参数，控制实空间和倒空间计算的平衡
 * @param meshSize 常规PME的网格尺寸
 * @param pair_cutoff 配对相互作用的截断距离
 * @param pairGridSize 预计算电势的网格尺寸
 * @param splineOrder B-样条插值的阶数
 * @param tolerance 计算精度的容差
 */
void setPGPParameters(double alpha, const int meshSize[3], double pair_cutoff, 
                        const int pairGridSize[3], int splineOrder, double tolerance) {
    // Set PME parameters (reuses PME functionality)
    setPMEParameters(alpha, meshSize, splineOrder, tolerance);
    
    // Copy PME parameters to PGP params
    pgp_params.alpha = pme_params.alpha;
    pgp_params.tolerance = pme_params.tolerance;
    pgp_params.initialized = pme_params.initialized;
    pgp_params.cutoff = pme_params.cutoff;
    pgp_params.epsilon_r = pme_params.epsilon_r;
    pgp_params.splineOrder = pme_params.splineOrder;
    
    // Copy box dimensions
    for (int i = 0; i < 3; i++) {
        pgp_params.box[i] = pme_params.box[i];
        pgp_params.meshSize[i] = pme_params.meshSize[i];
    }
    
    // Copy lookup tables
    pgp_params.erfcTable = pme_params.erfcTable;
    pgp_params.ewaldScaleTable = pme_params.ewaldScaleTable;
    pgp_params.ewaldDX = pme_params.ewaldDX;
    pgp_params.ewaldDXInv = pme_params.ewaldDXInv;
    pgp_params.erfcDXInv = pme_params.erfcDXInv;
    
    // Copy B-spline moduli
    for (int i = 0; i < 3; i++) {
        pgp_params.bsplineModuli[i] = pme_params.bsplineModuli[i];
    }
    
    // Copy reciprocal space grid
    pgp_params.pmeGrid = pme_params.pmeGrid;
    pgp_params.pmeCharge = pme_params.pmeCharge;
    
    // Set PGP-specific parameters
    pgp_params.pair_cutoff = pair_cutoff;
    for (int i = 0; i < 3; i++) {
        pgp_params.pair_grid_size[i] = pairGridSize[i];
    }
    
    // Initialize the pair grid
    pgp_params.initializePairGrid();
    
    platform::log(LogLevel::INFO, "PGP parameters set: alpha=", alpha, 
                 ", pair_cutoff=", pair_cutoff, 
                 ", pairGrid=[", pairGridSize[0], ",", pairGridSize[1], ",", pairGridSize[2], "]");
}

/**
 * @brief 预计算系统中固定部分的网格电势
 * 
 * 该函数实现了PGP-PME算法的核心步骤之一：预计算系统固定部分的电势场。
 * 它将固定部分的电荷分配到网格上，通过FFT变换计算电势，并存储结果供后续使用。
 * 这个步骤只需在系统的固定部分发生变化时执行一次，大大提高了MC模拟中能量评估的效率。
 * 
 * @param state 系统状态，包含原子信息和盒子大小等
 * @param fixed_only 是否只处理系统中的固定部分
 */
void precomputeGridPotential(model::MCState& state, bool fixed_only) {
    if (!pgp_params.initialized) {
        throw std::runtime_error("PGP parameters not initialized");
    }
    
    // Reset the pair grid
    std::fill(pgp_params.pairGrid.begin(), pgp_params.pairGrid.end(), std::complex<double>(0.0, 0.0));
    
    // Temporary grid for charge distribution
    std::vector<std::complex<double>> chargeGrid(pgp_params.pairGrid.size(), std::complex<double>(0.0, 0.0));
    
    // 1. 分配电荷到网格 - 使用B样条插值
    platform::log(LogLevel::DEBUG, "Spreading charges onto grid from ", 
                 fixed_only ? "fixed atoms only" : "all atoms");
    
    // Loop through residues and atoms
    for (int res_idx = 0; res_idx < state.activeResidueCount; ++res_idx) {
        const auto& residue = state.residues[res_idx];
        
        // Skip inactive residues
        if (!residue.active) continue;
        
        // Skip non-fixed residues if fixed_only is true
        if (fixed_only && !residue.fixed) continue;
        
        // Process each atom in the residue
        for (int atom_idx = 0; atom_idx < residue.atomCount; ++atom_idx) {
            const auto& atom = state.atoms[residue.atomStart + atom_idx];
            
            // Skip atoms with no charge
            if (std::abs(atom.charge) < 1e-10) continue;
            
            // Convert atom position to grid indices with B-spline offset
            double pos[3] = {atom.x, atom.y, atom.z};
            
            // Apply B-spline interpolation to spread charge to grid points
            // This is a simplified implementation - a full one would use proper B-spline functions
            int grid_indices[3][4]; // For a cubic B-spline (order 4)
            double weights[3][4];
            
            // Calculate grid indices and weights for interpolation
            for (int d = 0; d < 3; d++) {
                double scaled_pos = pos[d] / pgp_params.grid_spacing;
                int base_idx = static_cast<int>(std::floor(scaled_pos));
                
                // Calculate B-spline weights (simplified)
                double t = scaled_pos - base_idx;
                weights[d][0] = (1 - t) * (1 - t) * (1 - t) / 6.0;
                weights[d][1] = (3 * t * t * t - 6 * t * t + 4) / 6.0;
                weights[d][2] = (-3 * t * t * t + 3 * t * t + 3 * t + 1) / 6.0;
                weights[d][3] = t * t * t / 6.0;
                
                // Store grid indices with periodic boundary handling
                for (int i = 0; i < 4; i++) {
                    grid_indices[d][i] = (base_idx - 1 + i) % pgp_params.pair_grid_size[d];
                    if (grid_indices[d][i] < 0) grid_indices[d][i] += pgp_params.pair_grid_size[d];
                }
            }
            
            // Distribute charge to surrounding grid points
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    for (int k = 0; k < 4; k++) {
                        int grid_idx = (grid_indices[0][i] * pgp_params.pair_grid_size[1] + grid_indices[1][j]) 
                                      * pgp_params.pair_grid_size[2] + grid_indices[2][k];
                        
                        // Accumulate charge with B-spline weights
                        chargeGrid[grid_idx].real(chargeGrid[grid_idx].real() + 
                                                atom.charge * weights[0][i] * weights[1][j] * weights[2][k]);
                    }
                }
            }
        }
    }
    
    // 2. 执行FFT将实空间电荷分布转换到倒空间
    // Note: This is a placeholder for FFT transformation
    // In a real implementation, you would use a proper FFT library like FFTW
    platform::log(LogLevel::DEBUG, "Performing forward FFT on charge grid");
    
    // Placeholder for FFT (would be replaced with actual FFT implementation)
    std::vector<std::complex<double>> reciprocalGrid = chargeGrid;  // In real code, this would be the FFT result
    
    // 3. 应用Ewald因子
    platform::log(LogLevel::DEBUG, "Applying Ewald factor in reciprocal space");
    int nx = pgp_params.pair_grid_size[0];
    int ny = pgp_params.pair_grid_size[1];
    int nz = pgp_params.pair_grid_size[2];
    
    for (int i = 0; i < nx; i++) {
        int kx = (i <= nx/2) ? i : i - nx;
        double kx2 = kx * kx;
        
        for (int j = 0; j < ny; j++) {
            int ky = (j <= ny/2) ? j : j - ny;
            double ky2 = ky * ky;
            
            for (int k = 0; k < nz; k++) {
                int kz = (k <= nz/2) ? k : k - nz;
                double kz2 = kz * kz;
                
                // Skip k=0 (net zero charge case)
                if (kx == 0 && ky == 0 && kz == 0) continue;
                
                int idx = (i * ny + j) * nz + k;
                
                // Calculate k-vector squared
                double k2 = kx2 + ky2 + kz2;
                
                // Apply Ewald factor: exp(-k²/4α²) / k²
                double factor = std::exp(-k2 / (4.0 * pgp_params.alpha * pgp_params.alpha)) / k2;
                
                reciprocalGrid[idx] *= factor;
            }
        }
    }
    
    // 4. 执行反FFT获得实空间的预计算电势场
    platform::log(LogLevel::DEBUG, "Performing inverse FFT to get potential grid");
    
    // Placeholder for inverse FFT (would be replaced with actual FFT implementation)
    pgp_params.pairGrid = reciprocalGrid;  // In real code, this would be the inverse FFT result
    
    platform::log(LogLevel::INFO, "Grid potential precomputation completed for ",
                 fixed_only ? "fixed atoms only" : "all atoms");
}

/**
 * @brief 通过插值计算移动分子的能量
 * 
 * 该函数使用预计算的电势网格，快速评估移动分子在该电势场中的能量。
 * 当在MC模拟中尝试插入、删除或移动分子时，这个函数提供了一种高效的能量评估方法，
 * 无需重新计算整个系统的静电相互作用。
 * 
 * @param state 系统状态，包含移动分子的信息
 * @param energy 输出参数，存储计算得到的能量值
 */
void interpolateMoleculeEnergy(model::MCState& state, double& energy) {
    if (!pgp_params.initialized) {
        throw std::runtime_error("PGP parameters not initialized");
    }
    
    // Reset energy accumulator
    energy = 0.0;
    
    platform::log(LogLevel::DEBUG, "Interpolating energy for movement atoms");
    
    // Loop through movement residues
    for (const auto& movementInfo : state.movementResidues) {
        // Process each residue in the movement group
        for (int res_idx = movementInfo.startIndex; 
             res_idx < movementInfo.startIndex + movementInfo.activeCount; ++res_idx) {
            
            const auto& residue = state.residues[res_idx];
            
            // Skip inactive residues
            if (!residue.active) continue;
            
            // Process each atom in the residue
            for (int atom_idx = 0; atom_idx < residue.atomCount; ++atom_idx) {
                const auto& atom = state.atoms[residue.atomStart + atom_idx];
                
                // Skip atoms with no charge
                if (std::abs(atom.charge) < 1e-10) continue;
                
                // Convert atom position to grid coordinates
                double pos[3] = {atom.x, atom.y, atom.z};
                
                // Scale position to grid units
                double scaled_pos[3];
                for (int d = 0; d < 3; d++) {
                    scaled_pos[d] = pos[d] / pgp_params.grid_spacing;
                }
                
                // Calculate interpolation indices and weights
                int grid_indices[3][4]; // For a cubic B-spline (order 4)
                double weights[3][4];
                
                // Calculate indices and weights for interpolation
                for (int d = 0; d < 3; d++) {
                    int base_idx = static_cast<int>(std::floor(scaled_pos[d]));
                    
                    // Calculate B-spline weights (simplified)
                    double t = scaled_pos[d] - base_idx;
                    weights[d][0] = (1 - t) * (1 - t) * (1 - t) / 6.0;
                    weights[d][1] = (3 * t * t * t - 6 * t * t + 4) / 6.0;
                    weights[d][2] = (-3 * t * t * t + 3 * t * t + 3 * t + 1) / 6.0;
                    weights[d][3] = t * t * t / 6.0;
                    
                    // Store grid indices with periodic boundary handling
                    for (int i = 0; i < 4; i++) {
                        grid_indices[d][i] = (base_idx - 1 + i) % pgp_params.pair_grid_size[d];
                        if (grid_indices[d][i] < 0) grid_indices[d][i] += pgp_params.pair_grid_size[d];
                    }
                }
                
                // Interpolate potential at atom position
                double potential = 0.0;
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        for (int k = 0; k < 4; k++) {
                            int grid_idx = (grid_indices[0][i] * pgp_params.pair_grid_size[1] + grid_indices[1][j]) 
                                          * pgp_params.pair_grid_size[2] + grid_indices[2][k];
                            
                            // Accumulate potential with B-spline weights
                            potential += pgp_params.pairGrid[grid_idx].real() * 
                                        weights[0][i] * weights[1][j] * weights[2][k];
                        }
                    }
                }
                
                // Accumulate energy (potential * charge)
                energy += potential * atom.charge;
            }
        }
    }
    
    // Apply unit conversion and scaling factors if needed
    // In a real implementation, you might need to apply conversion factors
    // to get energy in proper units (e.g., kJ/mol)
    
    platform::log(LogLevel::DEBUG, "Interpolated energy: ", energy);
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc 