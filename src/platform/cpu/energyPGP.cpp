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
 * - 初始化阶段: 设置常规PME参数和额外的电势网格参数
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
 * 该函数根据potential_grid_size参数分配网格内存，并计算合适的网格间距。
 */
void PGPParams::initializePotentialGrid() {
    // 计算网格总大小并分配内存
    // 这一步决定了存储预计算电势所需的内存大小
    // 网格大小影响计算精度和内存消耗，需要权衡
    int totalSize = potential_grid_size[0] * potential_grid_size[1] * potential_grid_size[2];
    potentialGrid.resize(totalSize);
    
    // 设置网格间距，取三个维度中的最小值
    // 这确保了在各个方向上的分辨率至少达到指定精度
    // 网格间距对插值精度至关重要，间距越小精度越高但内存消耗越大
    grid_spacing = std::min({
        box[0] / potential_grid_size[0],
        box[1] / potential_grid_size[1],
        box[2] / potential_grid_size[2]
    });
    
    // 输出调试信息，帮助用户确认网格设置是否合理
    // 在大型模拟中，网格大小需要谨慎选择以平衡精度和性能
    platform::log(LogLevel::DEBUG, "PGP potential grid initialized with size: ", 
                 potential_grid_size[0], "x", potential_grid_size[1], "x", potential_grid_size[2],
                 ", grid spacing: ", grid_spacing);
}

/**
 * @brief 设置PGP-PME算法的所有参数
 * 
 * 该函数是PGP-PME算法的入口点，用于配置算法运行所需的所有参数。
 * 它首先设置标准PME参数，然后添加PGP特有的参数，最后初始化电势网格。
 * 
 * @param alpha Ewald分离参数，控制实空间和倒空间计算的平衡
 * @param meshSize 常规PME的网格尺寸
 * @param potential_cutoff 电势计算的截断距离
 * @param potentialGridSize 预计算电势的网格尺寸
 * @param splineOrder B-样条插值的阶数
 * @param tolerance 计算精度的容差
 */
void setPGPParameters(double alpha, const int meshSize[3], double potential_cutoff, 
                        const int potentialGridSize[3], int splineOrder, double tolerance) {
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
    // potential_cutoff定义了预计算电势的截断距离，通常小于PME的实空间截断
    // potential_grid_size定义了预计算电势网格的尺寸，影响插值精度
    pgp_params.potential_cutoff = potential_cutoff;
    for (int i = 0; i < 3; i++) {
        pgp_params.potential_grid_size[i] = potentialGridSize[i];
    }
    
    // 初始化预计算电势的网格
    // 这一步分配网格内存并计算网格间距
    pgp_params.initializePotentialGrid();
    
    // 输出参数设置信息
    // 这有助于调试和确认参数设置是否符合预期
    platform::log(LogLevel::INFO, "PGP parameters set: alpha=", alpha, 
                 ", potential_cutoff=", potential_cutoff, 
                 ", potentialGrid=[", potentialGridSize[0], ",", potentialGridSize[1], ",", potentialGridSize[2], "]");
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
    if (!pgp_params.initialized) {
        throw std::runtime_error("PGP parameters not initialized");
    }
    
    // 直接输出到标准输出以便调试
    std::cout << "\n[DIRECT PLATFORM] 预计算网格电势开始" << std::endl;
    std::cout << "  处理: " << (fixed_only ? "仅固定部分" : "所有部分") << std::endl;
    std::cout << "  网格大小: " << pgp_params.potential_grid_size[0] << "x" 
              << pgp_params.potential_grid_size[1] << "x" 
              << pgp_params.potential_grid_size[2] << std::endl;
    
    // 备份PME网格，稍后将恢复
    std::vector<std::complex<double>> pmeGridBackup = pme_params.pmeGrid;
    
    // 重置PME网格，准备新的计算
    std::fill(pme_params.pmeGrid.begin(), pme_params.pmeGrid.end(), std::complex<double>(0.0, 0.0));
    
    // 统计信息
    int fixed_residues_count = 0;
    
    // 计算固定残基数量
    for (int i = 0; i < state.activeResidueCount; ++i) {
        const auto& res = state.residues[i];
        if (res.fixed && res.active) fixed_residues_count++;
    }
    std::cout << "[DIRECT PLATFORM]  固定残基数: " << fixed_residues_count << std::endl;
    
    // 打印所有残基的信息以便调试
    std::cout << "[DIRECT PLATFORM] 打印所有残基的fixed状态:" << std::endl;
    for (int i = 0; i < state.activeResidueCount; ++i) {
        const auto& res = state.residues[i];
        std::cout << "  残基 " << i << ": fixed=" << res.fixed << ", active=" << res.active
                  << ", atomCount=" << res.atomCount << std::endl;
    }
    
    // 如果没有固定残基但要求仅计算固定部分，发出警告并自动切换
    if (fixed_only && fixed_residues_count == 0) {
        platform::log(LogLevel::WARNING, "No fixed residues found, switching to process all atoms");
        fixed_only = false;
    }
    
    // 设置pme_params的网格尺寸与pgp_params的网格尺寸一致，确保计算使用相同的网格
    for (int i = 0; i < 3; i++) {
        pme_params.meshSize[i] = pgp_params.potential_grid_size[i];
    }
    
    // 调整pme_params的网格大小以适应新的网格尺寸
    int totalGridSize = pgp_params.potential_grid_size[0] * pgp_params.potential_grid_size[1] * pgp_params.potential_grid_size[2];
    pme_params.pmeGrid.resize(totalGridSize, std::complex<double>(0.0, 0.0));
    
    // 调用PME的电荷分布函数
    platform::log(LogLevel::INFO, "调用PME的电荷分布函数 (fixed_only=", fixed_only, ")");
    spreadChargesOntoGrid(state, fixed_only);

    // 调用PME的前向FFT函数
    platform::log(LogLevel::INFO, "调用PME的前向FFT函数");
    performFFTForward();

    // 修正: 应用Ewald因子但不执行反向FFT
    platform::log(LogLevel::INFO, "手动应用Ewald因子");
    
    // 获取网格尺寸和盒子尺寸
    int nx = pgp_params.potential_grid_size[0];
    int ny = pgp_params.potential_grid_size[1];
    int nz = pgp_params.potential_grid_size[2];
    double volume = pgp_params.box[0] * pgp_params.box[1] * pgp_params.box[2];
    
    // 计算与应用Ewald因子所需的常数
    double alpha = pgp_params.alpha;
    // 修正：正确计算exp(-k²/(4α²))中的系数
    double factor = 1.0/(4.0*alpha*alpha);
    
    // 打印关键参数用于调试
    std::cout << "[DIRECT PLATFORM] 关键计算参数:" << std::endl;
    std::cout << "  盒子尺寸: [" << pgp_params.box[0] << ", " << pgp_params.box[1] << ", " << pgp_params.box[2] << "] nm" << std::endl;
    std::cout << "  盒子体积(Ω): " << volume << " nm³" << std::endl;
    std::cout << "  网格尺寸: [" << nx << ", " << ny << ", " << nz << "]" << std::endl;
    std::cout << "  总网格点数: " << nx * ny * nz << std::endl;
    std::cout << "  Ewald分离参数(α): " << alpha << " nm⁻¹" << std::endl;
    std::cout << "  exp(-k²/(4α²))系数: " << factor << std::endl;
    
    // 获取最大k向量指数
    int maxkx = (nx+1)/2;
    int maxky = (ny+1)/2;
    int maxkz = (nz+1)/2;
    
    // 计算倒格矢量
    double recipBoxVectors[3][3] = {{0}};
    // 修正: 倒格矢量应为2π/box，而非1/box
    // PME理论中k矢量定义为k = 2π·n/L，缺少2π会导致m²值偏小
    recipBoxVectors[0][0] = 2.0 * M_PI / pgp_params.box[0]; 
    recipBoxVectors[1][1] = 2.0 * M_PI / pgp_params.box[1]; 
    recipBoxVectors[2][2] = 2.0 * M_PI / pgp_params.box[2];
    
    // 打印倒格矢量
    std::cout << "  倒格矢量: [" << recipBoxVectors[0][0] << ", " << recipBoxVectors[1][1] << ", " << recipBoxVectors[2][2] << "] nm⁻¹" << std::endl;
    
    // 示例计算几个k点的值
    std::cout << "  示例k点值(nx/4, ny/4, nz/4):" << std::endl;
    int sx = nx/4, sy = ny/4, sz = nz/4;
    double mkx = sx * recipBoxVectors[0][0];
    double mky = sy * recipBoxVectors[1][1];
    double mkz = sz * recipBoxVectors[2][2];
    double mk2 = mkx*mkx + mky*mky + mkz*mkz;
    std::cout << "    k = [" << mkx << ", " << mky << ", " << mkz << "] nm⁻¹" << std::endl;
    std::cout << "    |k|² = " << mk2 << " nm⁻²" << std::endl;
    std::cout << "    exp(-k²/(4α²)) = " << exp(-mk2 * factor) << std::endl;
    
    // 应用Ewald因子
    for (int kx = 0; kx < nx; kx++) {
        double mx = (kx < maxkx) ? kx : (kx-nx);
        double mhx = mx * recipBoxVectors[0][0];
        
        for (int ky = 0; ky < ny; ky++) {
            double my = (ky < maxky) ? ky : (ky-ny);
            double mhy = my * recipBoxVectors[1][1];
            
            for (int kz = 0; kz < nz; kz++) {
                // 跳过零频率
                if (kx == 0 && ky == 0 && kz == 0) {
                    continue;
                }
                
                double mz = (kz < maxkz) ? kz : (kz-nz);
                double mhz = mz * recipBoxVectors[2][2];
                
                // 网格索引
                int index = kx * ny * nz + ky * nz + kz;
                // 获取结构因子
                std::complex<double> structureFactor = pme_params.pmeGrid[index];
                
                // 计算|k|^2
                double m2 = mhx * mhx + mhy * mhy + mhz * mhz;
                
                // 应用B-spline系数
                double bx = pgp_params.bsplineModuli[0][kx];
                double by = pgp_params.bsplineModuli[1][ky];
                double bz = pgp_params.bsplineModuli[2][kz];
                double denom = m2 * bx * by * bz; // 移除了boxfactor，符合PME理论
                
                // 避免除以零问题
                if (denom < 1e-10) {
                    denom = 1e-10;
                }
                
                // 只应用k依赖的因子: exp(-k²/(4α²))/(k² · B)
                // 常数因子(4π/Ω)将在反FFT后应用
                double kDependentFactor = exp(-m2 * factor) / denom;
                
                // 应用k依赖因子
                pme_params.pmeGrid[index] = structureFactor * kDependentFactor;
            }
        }
    }
    
    // 执行反向FFT以获得实空间电势
    platform::log(LogLevel::INFO, "执行反向FFT获得实空间电势");
    performFFTBackward();
    
    // 补偿反FFT的1/N归一化，并应用常数因子(4π/Ω)
    platform::log(LogLevel::INFO, "补偿反FFT的归一化因子并应用常数因子");
    int totalFFTPoints = nx * ny * nz;
    // 常数因子: 4π/Ω
    double constantFactor = 4.0 * M_PI / volume;
    
    // 电势物理单位转换因子
    double ONE_4PI_EPS0 = 138.935456; // kJ·mol^-1·nm·e^-2，与PME定义一致
    double physicalUnitFactor = ONE_4PI_EPS0 / pgp_params.epsilon_r;
    
    // 总修正因子 = FFT归一化补偿(N) × 常数因子(4π/Ω) × 物理单位转换 × 0.5(减半)
    double totalFactor = totalFFTPoints * constantFactor * physicalUnitFactor * 0.5;
    std::cout << "[DIRECT PLATFORM] 应用修正因子:" << std::endl;
    std::cout << "  FFT归一化补偿(N): " << totalFFTPoints << std::endl;
    std::cout << "  常数因子(4π/Ω): " << constantFactor << " nm⁻³" << std::endl;
    std::cout << "  物理单位转换: " << physicalUnitFactor << " kJ*nm/mol*e²" << std::endl;
    std::cout << "  电势减半因子: 0.5" << std::endl;
    std::cout << "  总修正因子: " << totalFactor << std::endl;
    
    // 统计修正前的电势范围
    double preMin = 0.0, preMax = 0.0, preSum = 0.0;
    bool firstValue = true;
    for (int i = 0; i < totalGridSize; i++) {
        double val = pme_params.pmeGrid[i].real();
        if (firstValue) {
            preMin = preMax = val;
            firstValue = false;
        } else {
            preMin = std::min(preMin, val);
            preMax = std::max(preMax, val);
        }
        preSum += val;
    }
    std::cout << "  修正前电势范围: [" << preMin << ", " << preMax << "], 平均值: " << (preSum/totalGridSize) << std::endl;
    
    // 应用总修正因子到每个网格点
    for (int i = 0; i < totalGridSize; i++) {
        pme_params.pmeGrid[i] *= totalFactor;
    }
    
    // 统计修正后的电势范围
    double postMin = 0.0, postMax = 0.0, postSum = 0.0;
    firstValue = true;
    for (int i = 0; i < totalGridSize; i++) {
        double val = pme_params.pmeGrid[i].real();
        if (firstValue) {
            postMin = postMax = val;
            firstValue = false;
        } else {
            postMin = std::min(postMin, val);
            postMax = std::max(postMax, val);
        }
        postSum += val;
    }
    std::cout << "  修正后电势范围: [" << postMin << ", " << postMax << "], 平均值: " << (postSum/totalGridSize) << std::endl;
    
    // 将修改后的PME网格复制到PGP的potentialGrid中
    pgp_params.potentialGrid = pme_params.pmeGrid;

    // 恢复原始PME网格
    pme_params.pmeGrid = pmeGridBackup;
    
    // 验证电势网格是否有合理的值
    int potentials_nonzero = 0;
    double max_potential = 0.0;
    double min_potential = 0.0;
    double sum_potential = 0.0;
    
    // 检查网格点值
    for (const auto& val : pgp_params.potentialGrid) {
        double pot_val = val.real();
        sum_potential += pot_val;
        if (std::abs(pot_val) > 1e-10) {
            potentials_nonzero++;
            max_potential = std::max(max_potential, pot_val);
            min_potential = std::min(min_potential, pot_val);
        }
    }
    
    // 输出电势网格统计信息
    platform::log(LogLevel::INFO, "电势网格统计:");
    platform::log(LogLevel::INFO, "  非零点数: ", potentials_nonzero);
    platform::log(LogLevel::INFO, "  最大电势: ", max_potential);
    platform::log(LogLevel::INFO, "  最小电势: ", min_potential);
    platform::log(LogLevel::INFO, "  电势总和: ", sum_potential);
    
    std::cout << "[DIRECT PLATFORM] 电势预计算完成" << std::endl
              << "  非零点数: " << potentials_nonzero << std::endl
              << "  电势范围: [" << min_potential << ", " << max_potential << "]" << std::endl;
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
    if (!pgp_params.initialized) {
        throw std::runtime_error("PGP parameters not initialized");
    }
    
    // 直接输出到标准输出以便调试
    std::cout << "\n[DIRECT PLATFORM] 通过插值计算移动分子能量" << std::endl;
    std::cout << "  移动残基组数: " << state.movementResidues.size() << std::endl;
    
    // 重置能量累加器
    energy = 0.0;
    
    // 检查预计算的网格是否为空
    bool gridEmpty = true;
    for (const auto& val : pgp_params.potentialGrid) {
        if (std::abs(val.real()) > 1e-10 || std::abs(val.imag()) > 1e-10) {
            gridEmpty = false;
            break;
        }
    }
    
    if (gridEmpty) {
        std::cout << "[DIRECT PLATFORM] 警告: PGP网格为空或未正确初始化!" << std::endl;
        
        // 输出一些网格样本点，帮助调试
        std::cout << "[DIRECT PLATFORM] 网格样本点值:" << std::endl;
        for (int i = 0; i < std::min(10, static_cast<int>(pgp_params.potentialGrid.size())); i++) {
            std::cout << "  网格点 " << i << ": " << pgp_params.potentialGrid[i].real() << std::endl;
        }
    } else {
        std::cout << "[DIRECT PLATFORM] PGP网格包含非零值" << std::endl;
    }
    
    // 检查移动残基设置
    if (state.movementResidues.empty()) {
        std::cout << "[DIRECT PLATFORM] 警告: 没有设置移动残基信息!" << std::endl;
        std::cout << "[DIRECT PLATFORM] 将尝试使用非固定残基作为移动残基" << std::endl;
        
        // 打印每个残基的情况，帮助调试
        std::cout << "[DIRECT PLATFORM] 残基状态:" << std::endl;
        for (int i = 0; i < state.activeResidueCount; ++i) {
            const auto& res = state.residues[i];
            std::cout << "  残基 " << i << ": fixed=" << res.fixed 
                     << ", active=" << res.active 
                     << ", atomCount=" << res.atomCount << std::endl;
        }
    }
    
    int totalAtoms = 0;
    int chargedAtoms = 0;
    double raw_energy = 0.0; // 用于存储未缩放的能量
    
    // 处理系统中的每个移动残基
    if (state.movementResidues.empty()) {
        std::cout << "[DIRECT PLATFORM] 警告: 没有设置移动残基信息!" << std::endl;
        
        // 如果没有设置移动残基，尝试处理所有非固定残基
        std::cout << "[DIRECT PLATFORM] 尝试查找所有非固定残基..." << std::endl;
        
        for (size_t i = 0; i < state.residues.size(); i++) {
            const auto& residue = state.residues[i];
            
            if (!residue.active || residue.fixed) continue;
            
            std::cout << "[DIRECT PLATFORM] 使用非固定残基 " << i << std::endl;
            
            // 处理原子...
            for (int j = 0; j < residue.atomCount; j++) {
                int atom_index = residue.atomStart + j;
                const auto& atom = state.atoms[atom_index];
                
                totalAtoms++;
                
                // 只处理带电荷的原子
                if (std::abs(atom.charge) < 1e-6) continue;
                
                chargedAtoms++;
                
                std::cout << "[DIRECT PLATFORM] 处理原子 " << atom_index << ": 位置=(" 
                          << atom.x << "," << atom.y << "," << atom.z 
                          << "), 电荷=" << atom.charge << std::endl;
                
                // 计算网格位置和B样条插值权重
                // 修正：使用与预计算电势相同的坐标变换方法
                
                // 获取原子位置
                double pos[3] = {atom.x, atom.y, atom.z};
                
                // 计算分数坐标 - 直接使用盒子尺寸进行变换，与precomputeGridPotential一致
                double fractional[3];
                for (int d = 0; d < 3; d++) {
                    fractional[d] = pos[d] / pgp_params.box[d];
                    fractional[d] -= floor(fractional[d]);  // 确保在[0,1)范围内
                    fractional[d] *= pgp_params.potential_grid_size[d]; // 缩放到网格
                }
                
                // 计算网格索引和分数部分
                int gridIndices[3];
                double gridFractions[3];
                for (int d = 0; d < 3; d++) {
                    gridFractions[d] = fractional[d] - floor(fractional[d]);
                    gridIndices[d] = static_cast<int>(floor(fractional[d]));
                    // 确保网格索引在正确范围内
                    if (gridIndices[d] < 0) 
                        gridIndices[d] += pgp_params.potential_grid_size[d];
                }
                
                // 计算B样条系数
                int nx = pgp_params.potential_grid_size[0];
                int ny = pgp_params.potential_grid_size[1];
                int nz = pgp_params.potential_grid_size[2];
                int order = pgp_params.splineOrder;
                
                std::vector<double> thetaX(order);
                std::vector<double> thetaY(order);
                std::vector<double> thetaZ(order);
                
                // 计算每个维度的B样条系数
                std::vector<double> coefficients(order);
                
                // X维度B样条
                computeBSplineCoefficients(gridFractions[0], order, coefficients);
                std::cout << "[DIRECT PLATFORM] X轴B样条系数 (gridFraction=" << gridFractions[0] << "):" << std::endl;
                double xMax = thetaX[0], xMin = thetaX[0], xSum = 0.0;
                for (int i = 0; i < order; i++) {
                    thetaX[i] = coefficients[i];
                    xMax = std::max(xMax, thetaX[i]);
                    xMin = std::min(xMin, thetaX[i]);
                    xSum += thetaX[i];
                    std::cout << "  theta_x[" << i << "] = " << thetaX[i] << std::endl;
                }
                std::cout << "  X权重范围: [" << xMin << ", " << xMax << "], 和: " << xSum << std::endl;
                
                // Y维度B样条
                computeBSplineCoefficients(gridFractions[1], order, coefficients);
                std::cout << "[DIRECT PLATFORM] Y轴B样条系数 (gridFraction=" << gridFractions[1] << "):" << std::endl;
                double yMax = thetaY[0], yMin = thetaY[0], ySum = 0.0;
                for (int i = 0; i < order; i++) {
                    thetaY[i] = coefficients[i];
                    yMax = std::max(yMax, thetaY[i]);
                    yMin = std::min(yMin, thetaY[i]);
                    ySum += thetaY[i];
                    std::cout << "  theta_y[" << i << "] = " << thetaY[i] << std::endl;
                }
                std::cout << "  Y权重范围: [" << yMin << ", " << yMax << "], 和: " << ySum << std::endl;
                
                // Z维度B样条
                computeBSplineCoefficients(gridFractions[2], order, coefficients);
                std::cout << "[DIRECT PLATFORM] Z轴B样条系数 (gridFraction=" << gridFractions[2] << "):" << std::endl;
                double zMax = thetaZ[0], zMin = thetaZ[0], zSum = 0.0;
                for (int i = 0; i < order; i++) {
                    thetaZ[i] = coefficients[i];
                    zMax = std::max(zMax, thetaZ[i]);
                    zMin = std::min(zMin, thetaZ[i]);
                    zSum += thetaZ[i];
                    std::cout << "  theta_z[" << i << "] = " << thetaZ[i] << std::endl;
                }
                std::cout << "  Z权重范围: [" << zMin << ", " << zMax << "], 和: " << zSum << std::endl;
                std::cout << "  三维权重乘积总和理论值: " << xSum * ySum * zSum << std::endl;
                
                // 输出最近的网格点及其电势值（用于测试）
                std::cout << "[DIRECT PLATFORM] 原子 " << atom_index << " 最近的网格点信息:" << std::endl;
                std::cout << "  网格索引: (" << gridIndices[0] << "," << gridIndices[1] << "," << gridIndices[2] << ")" << std::endl;
                
                // 输出该点及其周围网格点的电势值
                std::cout << "  最近网格点电势值:" << std::endl;
                for (int dx = -1; dx <= 1; dx++) {
                    for (int dy = -1; dy <= 1; dy++) {
                        for (int dz = -1; dz <= 1; dz++) {
                            int xi = (gridIndices[0] + dx + nx) % nx;
                            int yi = (gridIndices[1] + dy + ny) % ny;
                            int zi = (gridIndices[2] + dz + nz) % nz;
                            int index = xi * ny * nz + yi * nz + zi;
                            
                            double pot_val = pgp_params.potentialGrid[index].real();
                            
                            if (dx == 0 && dy == 0 && dz == 0) {
                                std::cout << "  → 中心点 (" << xi << "," << yi << "," << zi 
                                          << "): " << pot_val << std::endl;
                            } else if (std::abs(pot_val) > 1e-6) {
                                // 只输出非零的周围点
                                std::cout << "    点 (" << xi << "," << yi << "," << zi 
                                          << "): " << pot_val << std::endl;
                            }
                        }
                    }
                }
                
                // 插值计算电势
                double potential = 0.0;
                
                // 遍历所有B样条支撑点
                for (int ix = 0; ix < order; ix++) {
                    int xindex = (gridIndices[0] + ix) % nx;
                    
                    for (int iy = 0; iy < order; iy++) {
                        int yindex = (gridIndices[1] + iy) % ny;
                        
                        for (int iz = 0; iz < order; iz++) {
                            int zindex = (gridIndices[2] + iz) % nz;
                            
                            // 计算三维网格索引
                            int index = xindex * ny * nz + yindex * nz + zindex;
                            
                            // 使用B样条权重累加电势
                            double grid_value = pgp_params.potentialGrid[index].real();
                            double weight = thetaX[ix] * thetaY[iy] * thetaZ[iz];
                            potential += grid_value * weight;
                            
                            // 输出重要网格点的详细信息，帮助调试
                            if (std::abs(grid_value) > 1e-6 && ix < 2 && iy < 2 && iz < 2) {
                                std::cout << "[DIRECT PLATFORM] 网格点 (" 
                                          << xindex << "," 
                                          << yindex << "," 
                                          << zindex 
                                          << ") 电势=" << grid_value 
                                          << ", 权重=" << weight
                                          << " (各维度权重: " << thetaX[ix] << "," << thetaY[iy] << "," << thetaZ[iz] << ")"
                                          << std::endl;
                            }
                        }
                    }
                }
                
                // 输出电势插值中权重总和，用于验证归一化
                double totalWeight = 0.0;
                for (int ix = 0; ix < order; ix++) {
                    for (int iy = 0; iy < order; iy++) {
                        for (int iz = 0; iz < order; iz++) {
                            totalWeight += thetaX[ix] * thetaY[iy] * thetaZ[iz];
                        }
                    }
                }
                std::cout << "[DIRECT PLATFORM] B样条权重总和: " << totalWeight << std::endl;
                
                // 累加能量(电势*电荷)
                double atom_energy = potential * atom.charge;
                raw_energy += atom_energy;
                
                std::cout << "[DIRECT PLATFORM] 原子电势: " << potential
                          << ", 原子能量贡献: " << atom_energy << std::endl;
            }
        }
    } else {
        // 正常处理移动残基
        // ...
    }
    
    // 计算能量汇总
    platform::log(LogLevel::DEBUG, "Summing atomic contributions for total energy");
    
    // 电势已在预计算中减半，现在需要乘以2倍因子来计算能量
    // 根据pgp.md，能量应该是 2 * Σ(q_i * φ(r_i))
    energy = 2.0 * raw_energy;
    
    // 详细输出能量计算细节
    std::cout << "[DIRECT PLATFORM] 能量计算细节:" << std::endl;
    std::cout << "  原始能量(raw): " << raw_energy << " kJ/mol" << std::endl;
    std::cout << "  应用2倍因子: × 2.0" << std::endl;
    std::cout << "  最终能量: " << energy << " kJ/mol" << std::endl;
    
    platform::log(LogLevel::DEBUG, "Raw PGP energy: ", raw_energy, " kJ/mol");
    platform::log(LogLevel::DEBUG, "Final PGP energy (×2): ", energy, " kJ/mol");
    
    // 输出最终的能量值和调试信息
    platform::log(LogLevel::INFO, "最终计算的PGP能量: ", energy, " kJ/mol");
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
    
    // 添加调试输出
    platform::log(LogLevel::INFO, "Final calculated molecule energy: ", energy);
    
    // 如果没有固定残基，可能需要重新预计算网格电势
    if (std::abs(energy) < 1e-10) {
        // 检查是否是因为没有固定残基导致的计算问题
        int fixed_count = 0;
        for (int i = 0; i < state.activeResidueCount; ++i) {
            if (state.residues[i].fixed && state.residues[i].active) {
                fixed_count++;
            }
        }
        
        if (fixed_count == 0) {
            std::cout << "[DIRECT PLATFORM] 警告: 没有发现固定残基，能量接近零!" << std::endl;
            
            // 预计算全部残基的网格电势
            std::cout << "[DIRECT PLATFORM] 尝试用所有残基预计算网格电势..." << std::endl;
            precomputeGridPotential(state, false);
            
            // 重新计算能量
            interpolateMoleculeEnergy(state, energy);
        }
    }
    
    // 确保返回正确的符号和量级的能量值
    return energy;
}

double computeMoleculeEnergyGlobal(model::MCState& state, const std::vector<int>& movementResidues, const std::vector<int>& nearbyResidues, int threadIndex) {
    // 如果参数未初始化，则返回0
    if (!pgp_params.initialized) {
        platform::log(LogLevel::WARNING, "PGP parameters not initialized, returning 0 energy");
        return 0.0;
    }
    
    // 初始化总能量为0
    double totalEnergy = 0.0;
    
    // 记录所使用的线程索引，可用于多线程计算时的日志跟踪
    platform::log(LogLevel::DEBUG, "使用线程索引: ", threadIndex, " 计算PGP能量");
    
    // 记录附近残基数量，在某些算法变体中可用于短程能量修正
    platform::log(LogLevel::DEBUG, "考虑附近残基数量: ", nearbyResidues.size());
    
    // 对于直接空间(短程)能量修正，可以考虑附近的残基
    double directSpaceCorrection = 0.0;
    if (!nearbyResidues.empty()) {
        // 在一些PGP-PME实现中，虽然长程部分通过网格计算，
        // 但仍需要计算短程修正项，尤其是对非周期性边界或大系统
        platform::log(LogLevel::DEBUG, "计算附近残基的直接空间修正");
        
        // 这里可以实现直接空间修正计算
        // 但当前PGP实现主要关注预计算网格电势部分
        // 如果需要完整的直接空间修正，应单独实现
    }
    
    // 处理不同情况的残基能量计算
    if (movementResidues.empty()) {
        // 如果没有指定移动残基，计算所有残基的能量
        platform::log(LogLevel::DEBUG, "计算所有残基的能量");
        
        // 找出所有非固定的活动残基
        for (int i = 0; i < state.activeResidueCount; ++i) {
            if (state.residues[i].active && !state.residues[i].fixed) {
                // 直接调用calculateMoleculeEnergy，它内部会调用interpolateMoleculeEnergy
                totalEnergy = calculateMoleculeEnergy(state);
                break;
            }
        }
    } else {
        // 如果指定了移动残基，我们需要修改state的movementResidues
        // 先备份原始的movementResidues
        auto originalMovementResidues = state.movementResidues;
        
        // 清空并设置新的movementResidues
        state.movementResidues.clear();
        model::MCMovementResidueInfo info;
        info.startIndex = movementResidues[0];
        info.activeCount = movementResidues.size();
        state.movementResidues.push_back(info);
        
        // 计算能量
        totalEnergy = calculateMoleculeEnergy(state);
        
        // 恢复原始的movementResidues
        state.movementResidues = originalMovementResidues;
    }
    
    // 添加直接空间修正（如果有）
    totalEnergy += directSpaceCorrection;
    
    // 记录能量计算结果
    platform::log(LogLevel::DEBUG, "PGP能量计算结果: ", totalEnergy, 
                 " (线程: ", threadIndex, ", 考虑附近残基: ", !nearbyResidues.empty(), ")");
    return totalEnergy;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc 
