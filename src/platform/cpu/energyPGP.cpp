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
    
    platform::log(LogLevel::INFO, "PGP parameters set: alpha=", alpha, 
                 ", pair_cutoff=", pair_cutoff, 
                 ", pairGrid=[", pairGridSize[0], ",", pairGridSize[1], ",", pairGridSize[2], "]");
}

/**
 * @brief 自动调整PGP-PME参数
 * 
 * 根据系统大小和指定精度，自动设置最佳PGP-PME参数。
 * 该函数首先调用PME的自动参数调整，然后添加PGP所需的额外参数。
 * 网格大小会根据盒子尺寸和截断距离进行优化，以平衡计算效率和精度。
 * 
 * @param error_tolerance 计算精度的容差
 * @param cutoff_distance 实空间计算的截断距离
 * @param pair_cutoff 配对相互作用的截断距离
 * @param box 模拟盒子的尺寸
 */
void autoAdjustPGPParameters(double error_tolerance, double cutoff_distance, 
                               double pair_cutoff, const double box[3]) {
    // First, auto-adjust the PME parameters
    autoAdjustPMEParameters(error_tolerance, cutoff_distance, box);
    
    // Copy PME parameters to PGP params (same as in setPGPParameters)
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
    
    // Calculate appropriate pair grid size based on box and pair cutoff
    for (int i = 0; i < 3; i++) {
        // Simple heuristic: ratio of box size to cutoff, with minimum size
        pgp_params.pair_grid_size[i] = std::max(32, static_cast<int>(box[i] / pair_cutoff * 2.0));
        // Ensure it's a power of 2 for FFT efficiency
        pgp_params.pair_grid_size[i] = 1 << static_cast<int>(std::ceil(std::log2(pgp_params.pair_grid_size[i])));
    }
    
    // Initialize the pair grid
    pgp_params.initializePairGrid();
    
    platform::log(LogLevel::INFO, "PGP parameters auto-adjusted: ", 
                 "alpha=", pgp_params.alpha, 
                 ", pair_cutoff=", pgp_params.pair_cutoff, 
                 ", pairGrid=[", pgp_params.pair_grid_size[0], ",", 
                 pgp_params.pair_grid_size[1], ",", pgp_params.pair_grid_size[2], "]");
}

/**
 * @brief 将电荷分布扩散到配对网格上
 * 
 * PGP-PME算法的核心步骤之一。该函数将固定部分的电荷分布扩散到网格上，生成预计算的电势场。
 * 在标准使用中，这一步骤只需在系统的固定部分发生变化时执行，而不需要在每次MC尝试中重新计算。
 * 
 * @param state 模拟系统状态
 * @param movement_only 是否只处理移动部分的原子
 */
void spreadPairsOntoGrid([[maybe_unused]] model::MCState& state, [[maybe_unused]] bool movement_only) {
    // Reset the pair grid
    std::fill(pgp_params.pairGrid.begin(), pgp_params.pairGrid.end(), std::complex<double>(0.0, 0.0));
    
    // Implementation will depend on the specific PGP algorithm
    platform::log(LogLevel::DEBUG, "Spreading pairs onto grid. Movement only: ", movement_only);
    
    // TODO: Implement the pair grid spreading algorithm
    // 1. 分离系统中的固定部分和移动部分
    // 2. 使用B样条插值将固定部分的电荷扩散到网格上
    // 3. 执行FFT将实空间电荷分布转换到倒空间
    // 4. 应用Ewald因子
    // 5. 执行反FFT获得实空间的预计算电势场
}

/**
 * @brief 从预计算的配对网格中计算能量
 * 
 * 使用预计算的电势网格快速评估移动分子的能量。这是PGP-PME算法的核心优势部分。
 * 在MC模拟中，当尝试插入或移动分子时，不需要重新计算整个系统的静电相互作用，
 * 而是通过查询预计算的电势网格快速获得能量变化。
 * 
 * @param energy 输出参数，存储计算得到的能量
 */
void computeEnergyFromPairGrid([[maybe_unused]] double& energy) {
    // Implementation will depend on the specific PGP algorithm
    platform::log(LogLevel::DEBUG, "Computing energy from pair grid");
    
    // TODO: Implement energy calculation from pair grid
    // 1. 将移动分子的电荷通过B样条插值映射到网格点上
    // 2. 计算移动分子电荷与预计算电势的乘积
    // 3. 对所有网格点求和获得总能量
}

/**
 * @brief 使用PGP方法计算倒空间能量
 * 
 * 该函数结合了传统PME的倒空间计算和PGP特有的配对网格贡献。
 * 对于长程相互作用，复用PME的高效计算；对于中程相互作用，
 * 使用预计算的配对网格加速计算。
 * 
 * @param state 模拟系统状态
 * @param movement_only 是否只计算移动部分的能量
 * @return 倒空间总能量
 */
double computeReciprocalPGP(model::MCState& state, bool movement_only) {
    // Reuse the PME reciprocal calculation for the long-range part
    double reciprocal_energy = computeReciprocalPME(state, movement_only);
    
    // Add the pair-grid contribution
    spreadPairsOntoGrid(state, movement_only);
    double pair_grid_energy = 0.0;
    computeEnergyFromPairGrid(pair_grid_energy);
    
    platform::log(LogLevel::DEBUG, "PGP reciprocal energy: PME=", reciprocal_energy, 
                 ", pair-grid=", pair_grid_energy, 
                 ", total=", reciprocal_energy + pair_grid_energy);
    
    return reciprocal_energy + pair_grid_energy;
}

// Compute self energy (mostly the same as PME)
double computeSelfEnergyPGP(model::MCState& state, bool movement_only) {
    // The self energy calculation is the same as in PME
    return computeSelfEnergyPME(state, movement_only);
}

// Compute real space energy using the pair grid approach
void computeRealSpacePGP(model::MCState& state, bool movement_only, bool store_in_residues) {
    // This is a modified version of the real space calculation
    // TODO: Implement the PGP real space calculation
    
    // For now, just use the PME real space calculation
    computeRealSpacePME(state, movement_only, store_in_residues);
    
    platform::log(LogLevel::DEBUG, "PGP real space energy calculated");
}

// Compute pair grid interactions
void computePairGridPGP([[maybe_unused]] model::MCState& state, [[maybe_unused]] bool movement_only) {
    // This is a new function specific to PGP
    // TODO: Implement the pair grid calculations
    
    platform::log(LogLevel::DEBUG, "PGP pair grid calculation");
}

/**
 * @brief 计算整个系统的PGP能量
 * 
 * 这是PGP-PME方法的主要入口点，用于计算整个系统的能量。
 * 该函数协调实空间计算、倒空间计算和自能项计算，并将结果
 * 存储在系统状态中。
 * 
 * @param state 模拟系统状态
 */
void computeSystemEnergyPGP(model::MCState& state) {
    if (!pgp_params.initialized) {
        throw std::runtime_error("PGP parameters not initialized");
    }
    
    // Check system neutrality
    checkSystemNeutrality(state);
    
    // Reset the ewald_energy structure
    state.ewald_energy = {0.0, 0.0, 0.0};
    
    // Compute real space energy (includes vdw)
    computeRealSpacePGP(state, false, true);
    
    // Compute reciprocal space energy
    state.ewald_energy.reciprocal = computeReciprocalPGP(state, false);
    
    // Compute self energy
    state.ewald_energy.self = computeSelfEnergyPGP(state, false);
    
    platform::log(LogLevel::DEBUG, "PGP system energy: real=", state.ewald_energy.real_space,
                 ", reciprocal=", state.ewald_energy.reciprocal,
                 ", self=", state.ewald_energy.self,
                 ", total=", state.ewald_energy.real_space + state.ewald_energy.reciprocal + state.ewald_energy.self);
}

/**
 * @brief 计算移动残基的PGP能量
 * 
 * 这是PGP-PME方法用于MC模拟中能量评估的关键函数。它只计算
 * 标记为移动的残基的能量变化，大大提高了MC模拟的效率。
 * 在GCMC和CBMC场景下特别有用，例如当尝试插入新分子或重构分子构型时。
 * 
 * @param state 模拟系统状态，包含移动残基的信息
 */
void computeMovementEnergyPGP(model::MCState& state) {
    if (!pgp_params.initialized) {
        throw std::runtime_error("PGP parameters not initialized");
    }
    
    // Check system neutrality
    checkSystemNeutrality(state);
    
    // Reset the ewald_energy structure
    state.ewald_energy = {0.0, 0.0, 0.0};
    
    // Compute real space energy for movement residues
    computeRealSpacePGP(state, true, true);
    
    // Compute reciprocal space energy for movement
    state.ewald_energy.reciprocal = computeReciprocalPGP(state, true);
    
    // Compute self energy for movement
    state.ewald_energy.self = computeSelfEnergyPGP(state, true);
    
    platform::log(LogLevel::DEBUG, "PGP movement energy: real=", state.ewald_energy.real_space,
                 ", reciprocal=", state.ewald_energy.reciprocal,
                 ", self=", state.ewald_energy.self,
                 ", total=", state.ewald_energy.real_space + state.ewald_energy.reciprocal + state.ewald_energy.self);
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc 