#include "energyPGP.hpp"
#include "energyPME.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

/**
 * @file energyPGP.cpp
 * @brief Implementation of Precomputed Grid-Potential Particle Mesh Ewald for Monte Carlo (PGP-PME-MC)
 * 
 * PGP-PME (Precomputed Grid-Potential Particle Mesh Ewald)是一种为蒙特卡洛模拟优化的粒子网格Ewald方法。
 * 它通过预计算系统中固定部分的静电势网格，大大加速了在MC模拟中计算小分子能量变化的过程。
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
 * @brief 初始化预计算电势的三维网格
 * 
 * 这是PGP-PME算法的基础步骤，负责创建和初始化用于存储预计算电势的三维网格。
 * 该函数根据pair_grid_size参数分配网格内存，并计算合适的网格间距。
 */
void PGPParams::initializePairGrid() {
    // 计算网格总大小并分配内存
    // 这一步决定了存储预计算电势所需的内存大小
    // 网格大小影响计算精度和内存消耗，需要权衡
    int totalSize = pair_grid_size[0] * pair_grid_size[1] * pair_grid_size[2];
    pairGrid.resize(totalSize);
    
    // 设置网格间距，取三个维度中的最小值
    // 这确保了在各个方向上的分辨率至少达到指定精度
    // 网格间距对插值精度至关重要，间距越小精度越高但内存消耗越大
    grid_spacing = std::min({
        box[0] / pair_grid_size[0],
        box[1] / pair_grid_size[1],
        box[2] / pair_grid_size[2]
    });
    
    // 输出调试信息，帮助用户确认网格设置是否合理
    // 在大型模拟中，网格大小需要谨慎选择以平衡精度和性能
    platform::log(LogLevel::DEBUG, "PGP pair grid initialized with size: ", 
                 pair_grid_size[0], "x", pair_grid_size[1], "x", pair_grid_size[2],
                 ", grid spacing: ", grid_spacing);
}

/**
 * @brief 设置PGP-PME算法的所有参数
 * 
 * 该函数是PGP-PME算法的入口点，用于配置算法运行所需的所有参数。
 * 它首先设置标准PME参数，然后添加PGP特有的参数，最后初始化配对网格。
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
    // 首先设置标准PME参数
    // 这里复用了PME的参数设置函数，避免代码重复
    // 标准PME参数包括alpha、网格大小、样条阶数和误差容限
    setPMEParameters(alpha, meshSize, splineOrder, tolerance);
    
    // 将标准PME参数复制到PGP参数结构中
    // 这确保了PGP方法与PME方法使用相同的基础参数
    // 参数复制是深拷贝，避免意外修改原始PME参数
    pgp_params.alpha = pme_params.alpha;
    pgp_params.tolerance = pme_params.tolerance;
    pgp_params.initialized = pme_params.initialized;
    pgp_params.cutoff = pme_params.cutoff;
    pgp_params.epsilon_r = pme_params.epsilon_r;
    pgp_params.splineOrder = pme_params.splineOrder;
    
    // 复制盒子尺寸和网格尺寸
    // 盒子尺寸决定了周期性边界条件的应用方式
    // 网格尺寸影响倒空间计算的精度
    for (int i = 0; i < 3; i++) {
        pgp_params.box[i] = pme_params.box[i];
        pgp_params.meshSize[i] = pme_params.meshSize[i];
    }
    
    // 复制PME查找表
    // 这些表格用于加速erfc函数和Ewald缩放因子的计算
    // 预计算表格可以大幅减少运行时的计算量
    pgp_params.erfcTable = pme_params.erfcTable;
    pgp_params.ewaldScaleTable = pme_params.ewaldScaleTable;
    pgp_params.ewaldDX = pme_params.ewaldDX;
    pgp_params.ewaldDXInv = pme_params.ewaldDXInv;
    pgp_params.erfcDXInv = pme_params.erfcDXInv;
    
    // 复制B样条模块
    // B样条模块用于在网格上平滑分配电荷和插值电势
    // 模块复制确保了PGP方法使用与PME相同的插值精度
    for (int i = 0; i < 3; i++) {
        pgp_params.bsplineModuli[i] = pme_params.bsplineModuli[i];
    }
    
    // 复制PME网格
    // 这些网格用于标准PME计算，PGP会在此基础上添加额外的预计算网格
    pgp_params.pmeGrid = pme_params.pmeGrid;
    pgp_params.pmeCharge = pme_params.pmeCharge;
    
    // 设置PGP特有参数
    // pair_cutoff定义了预计算电势的截断距离，通常小于PME的实空间截断
    // pair_grid_size定义了预计算电势网格的尺寸，影响插值精度
    pgp_params.pair_cutoff = pair_cutoff;
    for (int i = 0; i < 3; i++) {
        pgp_params.pair_grid_size[i] = pairGridSize[i];
    }
    
    // 初始化预计算电势的网格
    // 这一步分配网格内存并计算网格间距
    pgp_params.initializePairGrid();
    
    // 输出参数设置信息
    // 这有助于调试和确认参数设置是否符合预期
    platform::log(LogLevel::INFO, "PGP parameters set: alpha=", alpha, 
                 ", pair_cutoff=", pair_cutoff, 
                 ", pairGrid=[", pairGridSize[0], ",", pairGridSize[1], ",", pairGridSize[2], "]");
}

/**
 * @brief 预计算系统中固定部分的网格电势
 * 
 * 这是PGP-PME算法的核心函数之一，负责预计算系统中固定部分的静电势场。
 * 该函数将固定部分的电荷分配到网格上，通过FFT变换计算电势，并存储结果供后续使用。
 * 预计算步骤只需在系统固定部分发生变化时执行，显著提高了MC模拟的效率。
 * 
 * @param state 系统状态，包含原子坐标、电荷和盒子信息
 * @param fixed_only 是否只处理系统中的固定部分(true)，还是处理所有部分(false)
 */
void precomputeGridPotential(model::MCState& state, bool fixed_only) {
    // 检查参数是否已初始化
    // 这是安全检查，确保在使用前已正确设置了PGP参数
    if (!pgp_params.initialized) {
        throw std::runtime_error("PGP parameters not initialized");
    }
    
    // 重置配对网格
    // 在重新计算前清空网格，避免旧数据的影响
    std::fill(pgp_params.pairGrid.begin(), pgp_params.pairGrid.end(), std::complex<double>(0.0, 0.0));
    
    // 创建临时电荷网格
    // 这个网格用于存储电荷分布，之后会通过FFT转换为电势
    std::vector<std::complex<double>> chargeGrid(pgp_params.pairGrid.size(), std::complex<double>(0.0, 0.0));
    
    // 1. 分配电荷到网格 - 使用B样条插值
    // 这一步将原子点电荷平滑地分配到网格点上
    // 日志输出帮助跟踪处理的是固定部分还是全部原子
    platform::log(LogLevel::DEBUG, "Spreading charges onto grid from ", 
                 fixed_only ? "fixed atoms only" : "all atoms");
    
    // 遍历所有残基和原子
    // 这个循环是最耗时的部分之一，处理每个原子的电荷分配
    for (int res_idx = 0; res_idx < state.activeResidueCount; ++res_idx) {
        const auto& residue = state.residues[res_idx];
        
        // 跳过非活跃残基
        // 非活跃残基不参与能量计算
        if (!residue.active) continue;
        
        // 如果只处理固定部分，则跳过非固定残基
        // 这是PGP方法的关键优化点，只预计算固定部分的电势
        if (fixed_only && !residue.fixed) continue;
        
        // 处理残基中的每个原子
        // 循环遍历残基中的所有原子，将其电荷分配到网格上
        for (int atom_idx = 0; atom_idx < residue.atomCount; ++atom_idx) {
            const auto& atom = state.atoms[residue.atomStart + atom_idx];
            
            // 跳过无电荷原子
            // 无电荷原子不贡献静电势，可以跳过以提高效率
            if (std::abs(atom.charge) < 1e-10) continue;
            
            // 将原子位置转换为网格索引
            // 这个转换考虑了B样条插值需要的偏移
            double pos[3] = {atom.x, atom.y, atom.z};
            
            // 应用B样条插值将电荷分布到网格点上
            // 这是一个简化实现，完整实现应使用完整的B样条函数
            // B样条插值确保了电荷分布的平滑性和连续性
            int grid_indices[3][4]; // 用于三次B样条(阶数4)
            double weights[3][4];
            
            // 计算插值的网格索引和权重
            // 这个循环为三个维度分别计算
            for (int d = 0; d < 3; d++) {
                // 将原子位置缩放到网格单位
                // 这个转换考虑了网格间距
                double scaled_pos = pos[d] / pgp_params.grid_spacing;
                int base_idx = static_cast<int>(std::floor(scaled_pos));
                
                // 计算B样条权重(简化版)
                // 这些权重决定了原子电荷如何分配到相邻网格点
                // 三次B样条使用四个点进行插值
                double t = scaled_pos - base_idx;
                weights[d][0] = (1 - t) * (1 - t) * (1 - t) / 6.0;
                weights[d][1] = (3 * t * t * t - 6 * t * t + 4) / 6.0;
                weights[d][2] = (-3 * t * t * t + 3 * t * t + 3 * t + 1) / 6.0;
                weights[d][3] = t * t * t / 6.0;
                
                // 存储网格索引并处理周期性边界
                // 周期性边界处理确保了即使原子在盒子边缘也能正确计算
                for (int i = 0; i < 4; i++) {
                    grid_indices[d][i] = (base_idx - 1 + i) % pgp_params.pair_grid_size[d];
                    if (grid_indices[d][i] < 0) grid_indices[d][i] += pgp_params.pair_grid_size[d];
                }
            }
            
            // 将电荷分配到周围网格点
            // 这是三维B样条插值的核心，将电荷按权重分配到64个相邻网格点
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    for (int k = 0; k < 4; k++) {
                        // 计算三维网格索引
                        // 将三维索引转换为一维数组索引
                        int grid_idx = (grid_indices[0][i] * pgp_params.pair_grid_size[1] + grid_indices[1][j]) 
                                      * pgp_params.pair_grid_size[2] + grid_indices[2][k];
                        
                        // 按B样条权重累加电荷
                        // 原子电荷乘以三个维度的权重积
                        chargeGrid[grid_idx].real(chargeGrid[grid_idx].real() + 
                                                atom.charge * weights[0][i] * weights[1][j] * weights[2][k]);
                    }
                }
            }
        }
    }
    
    // 2. 执行FFT将实空间电荷分布转换到倒空间
    // 使用CustomFFT进行前向FFT变换
    platform::log(LogLevel::DEBUG, "Performing forward FFT on charge grid");
    
    // 使用CustomFFT执行前向FFT变换
    // 从energyPME.cpp中使用同样的FFT实现
    int nx = pgp_params.pair_grid_size[0];
    int ny = pgp_params.pair_grid_size[1];
    int nz = pgp_params.pair_grid_size[2];
    
    // 备份chargeGrid以便调试
    std::vector<std::complex<double>> reciprocalGrid = chargeGrid;
    
    // 执行3D前向FFT变换
    CustomFFT::fft3D_forward(reciprocalGrid.data(), nx, ny, nz);
    
    // 3. 应用Ewald因子
    // 在倒空间中对每个k向量应用Ewald因子
    // 这一步计算了长程静电相互作用
    platform::log(LogLevel::DEBUG, "Applying Ewald factor in reciprocal space");
    
    // 遍历所有倒空间网格点
    // 这个三重循环应用了Ewald因子到每个k向量
    for (int i = 0; i < nx; i++) {
        // 计算k向量的x分量
        // 考虑了周期性边界条件
        int kx = (i <= nx/2) ? i : i - nx;
        double kx2 = kx * kx;
        
        for (int j = 0; j < ny; j++) {
            // 计算k向量的y分量
            int ky = (j <= ny/2) ? j : j - ny;
            double ky2 = ky * ky;
            
            for (int k = 0; k < nz; k++) {
                // 计算k向量的z分量
                int kz = (k <= nz/2) ? k : k - nz;
                double kz2 = kz * kz;
                
                // 跳过k=0(净零电荷情况)
                // k=0对应于系统的总电荷，通常需要特殊处理
                if (kx == 0 && ky == 0 && kz == 0) continue;
                
                // 计算一维数组索引
                int idx = (i * ny + j) * nz + k;
                
                // 计算k向量的平方
                // 这决定了Ewald因子的大小
                double k2 = kx2 + ky2 + kz2;
                
                // 应用Ewald因子: exp(-k²/4α²) / k²
                // 这是PME方法中的标准Ewald因子
                // 它控制了倒空间中不同波长的贡献
                double factor = std::exp(-k2 / (4.0 * pgp_params.alpha * pgp_params.alpha)) / k2;
                
                // 将因子应用到倒空间网格点
                reciprocalGrid[idx] *= factor;
            }
        }
    }
    
    // 4. 执行反FFT获得实空间的预计算电势场
    // 反FFT将倒空间的电势转换回实空间
    platform::log(LogLevel::DEBUG, "Performing inverse FFT to get potential grid");
    
    // 使用CustomFFT执行反向FFT变换
    CustomFFT::fft3D_backward(reciprocalGrid.data(), nx, ny, nz);
    
    // 保存结果到预计算电势网格
    pgp_params.pairGrid = reciprocalGrid;
    
    // 记录预计算完成信息
    // 这有助于确认预计算过程已成功完成
    platform::log(LogLevel::INFO, "Grid potential precomputation completed for ",
                 fixed_only ? "fixed atoms only" : "all atoms");
}

/**
 * @brief 通过插值计算移动分子的能量
 * 
 * 该函数是PGP-PME算法的另一个核心函数，用于在预计算的电势场中快速评估移动分子的能量。
 * 通过B样条插值从预计算的网格电势中获取能量，避免了直接计算分子间相互作用，
 * 大大提高了MC模拟中能量评估的效率。
 * 
 * @param state 系统状态，包含移动分子的信息和预计算的电势网格
 * @param energy 输出参数，存储计算得到的能量值
 */
void interpolateMoleculeEnergy(model::MCState& state, double& energy) {
    // 检查参数是否已初始化
    // 这是安全检查，确保在使用前已正确设置了PGP参数
    if (!pgp_params.initialized) {
        throw std::runtime_error("PGP parameters not initialized");
    }
    
    // 重置能量累加器
    // 初始化能量为0，准备累加每个原子的贡献
    energy = 0.0;
    
    // 记录开始计算的调试信息
    platform::log(LogLevel::DEBUG, "Interpolating energy for movement atoms");
    
    // 遍历移动残基
    // 这个循环只处理标记为移动的残基，大大减少了计算量
    for (const auto& movementInfo : state.movementResidues) {
        // 处理移动组中的每个残基
        // 一个移动组可能包含多个残基
        for (int res_idx = movementInfo.startIndex; 
             res_idx < movementInfo.startIndex + movementInfo.activeCount; ++res_idx) {
            
            const auto& residue = state.residues[res_idx];
            
            // 跳过非活跃残基
            // 只计算活跃残基的能量
            if (!residue.active) continue;
            
            // 处理残基中的每个原子
            // 循环计算残基中每个原子的能量贡献
            for (int atom_idx = 0; atom_idx < residue.atomCount; ++atom_idx) {
                const auto& atom = state.atoms[residue.atomStart + atom_idx];
                
                // 跳过无电荷原子
                // 无电荷原子不贡献静电能量
                if (std::abs(atom.charge) < 1e-10) continue;
                
                // 将原子位置转换为网格坐标
                // 提取原子的三维坐标
                double pos[3] = {atom.x, atom.y, atom.z};
                
                // 将位置缩放到网格单位
                // 这个转换考虑了网格间距
                double scaled_pos[3];
                for (int d = 0; d < 3; d++) {
                    scaled_pos[d] = pos[d] / pgp_params.grid_spacing;
                }
                
                // 计算插值索引和权重
                // 使用与电荷分配相同的B样条插值方法
                int grid_indices[3][4]; // 用于三次B样条(阶数4)
                double weights[3][4];
                
                // 计算插值的索引和权重
                // 这个循环为三个维度分别计算
                for (int d = 0; d < 3; d++) {
                    int base_idx = static_cast<int>(std::floor(scaled_pos[d]));
                    
                    // 计算B样条权重(简化版)
                    // 这些权重用于电势插值
                    double t = scaled_pos[d] - base_idx;
                    weights[d][0] = (1 - t) * (1 - t) * (1 - t) / 6.0;
                    weights[d][1] = (3 * t * t * t - 6 * t * t + 4) / 6.0;
                    weights[d][2] = (-3 * t * t * t + 3 * t * t + 3 * t + 1) / 6.0;
                    weights[d][3] = t * t * t / 6.0;
                    
                    // 存储网格索引并处理周期性边界
                    // 确保索引在有效范围内
                    for (int i = 0; i < 4; i++) {
                        grid_indices[d][i] = (base_idx - 1 + i) % pgp_params.pair_grid_size[d];
                        if (grid_indices[d][i] < 0) grid_indices[d][i] += pgp_params.pair_grid_size[d];
                    }
                }
                
                // 在原子位置插值电势
                // 这是PGP方法的核心优势，通过插值快速获取电势
                double potential = 0.0;
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        for (int k = 0; k < 4; k++) {
                            // 计算三维网格索引
                            // 将三维索引转换为一维数组索引
                            int grid_idx = (grid_indices[0][i] * pgp_params.pair_grid_size[1] + grid_indices[1][j]) 
                                          * pgp_params.pair_grid_size[2] + grid_indices[2][k];
                            
                            // 使用B样条权重累加电势
                            // 三个维度的权重乘积乘以网格点的电势值
                            potential += pgp_params.pairGrid[grid_idx].real() * 
                                        weights[0][i] * weights[1][j] * weights[2][k];
                        }
                    }
                }
                
                // 累加能量(电势*电荷)
                // 能量等于电势乘以电荷
                energy += potential * atom.charge;
            }
        }
    }
    
    // 应用单位转换和缩放因子(如需要)
    // 在实际实现中，可能需要应用转换因子以获得正确单位的能量(例如，kJ/mol)
    // 例如: energy *= ONE_4PI_EPS0 / pgp_params.epsilon_r;
    
    // 记录计算结果
    platform::log(LogLevel::DEBUG, "Interpolated energy: ", energy);
}

/**
 * @brief 通过插值计算移动分子的能量，并返回计算结果
 * 
 * 这是interpolateMoleculeEnergy的包装函数，直接返回计算得到的能量值
 * 方便Python调用和测试。
 * 
 * @param state 系统状态，包含移动分子的信息及预计算的电势网格
 * @return 计算得到的能量值
 */
double calculateMoleculeEnergy(model::MCState& state) {
    double energy = 0.0;
    interpolateMoleculeEnergy(state, energy);
    return energy;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc 