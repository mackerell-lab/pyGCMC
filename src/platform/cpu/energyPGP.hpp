#pragma once

#include "model/montecarlo.hpp"
#include "platform/platform.hpp"
#include "energyCommon.hpp"
#include <array>
#include <vector>
#include <complex>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief PGP-PME算法参数结构体 (Precomputed Grid-Potential Particle Mesh Ewald)
 * 
 * 该结构体包含了Precomputed Grid-Potential Particle Mesh Ewald算法所需的所有参数和数据结构。
 * PGP-PME是一种优化的PME方法，通过预计算电势网格加速蒙特卡洛模拟中的能量评估。
 * 
 * 主要包含以下几类参数：
 * 1. 常规PME参数：alpha、网格大小、插值阶数等
 * 2. 预计算网格参数：配对截断距离、网格大小和间距等
 * 3. 数据存储：预计算电势网格、B样条模数等
 */
struct PGPParams {
    bool initialized = false;      // 参数是否已初始化
    double alpha;                  // Ewald分离参数，平衡实空间和倒空间计算
    double tolerance;              // 误差容限
    double cutoff;                 // 实空间截断距离
    double epsilon_r;              // 相对介电常数
    int splineOrder;               // B样条插值阶数(通常为4，即三次B样条)
    std::array<double, 3> box;     // 模拟盒子的三维尺寸
    std::array<int, 3> meshSize;   // PME网格的三维尺寸
    
    // 预计算电势网格参数
    double pair_cutoff;                    // 配对相互作用截断距离
    std::array<int, 3> pair_grid_size;     // 预计算电势网格尺寸
    double grid_spacing;                   // 网格间距
    std::vector<std::complex<double>> pairGrid;  // 预计算电势网格数据
    
    // PME算法参数(主要由setPMEParameters函数设置)
    std::vector<double> erfcTable;         // erfc函数查找表
    std::vector<double> ewaldScaleTable;   // Ewald缩放因子查找表
    double ewaldDX;                        // Ewald表步长
    double ewaldDXInv;                     // Ewald表步长倒数
    double erfcDXInv;                      // erfc表步长倒数
    std::vector<double> bsplineModuli[3];  // B样条模数
    std::vector<std::complex<double>> pmeGrid;   // PME网格
    std::vector<double> pmeCharge;        // PME电荷网格，类型需与PME结构体匹配

    /**
     * @brief 初始化预计算电势的三维网格
     * 
     * 该方法创建并初始化用于存储预计算电势的三维网格结构。网格大小由pair_grid_size参数决定，
     * 通常根据所需精度和计算资源进行设置。网格间距自动计算，以确保在所有维度上获得足够的分辨率。
     * 
     * 算法原理：
     * 1. 根据指定的网格维度计算总网格点数
     * 2. 分配内存空间用于存储预计算的电势值
     * 3. 计算网格间距，确保在所有维度上至少达到所需分辨率
     * 
     * 计算复杂度：
     * - 空间复杂度: O(nx*ny*nz)，其中nx/ny/nz为网格在各维度的大小
     * - 时间复杂度: O(1)
     * 
     * 使用场景：
     * - 在设置PGP参数后自动调用
     * - 当系统尺寸变化时需要重新初始化
     */
    void initializePairGrid();
};

// Global PGP parameters
extern PGPParams pgp_params;

/**
 * @brief 设置PGP-PME (Precomputed Grid-Potential Particle Mesh Ewald)算法的所有参数
 * 
 * 该函数是PGP-PME算法配置的入口点，设置所有运行PGP-PME所需的参数。
 * 它首先配置标准PME参数，然后添加PGP特有的网格和截断参数，最后初始化预计算网格。
 * 
 * 算法原理：
 * 1. 调用setPMEParameters设置基础PME参数
 * 2. 将PME参数复制到PGP参数结构中
 * 3. 添加PGP特有参数如配对截断和网格大小
 * 4. 初始化预计算电势网格结构
 * 
 * 参数含义：
 * @param alpha Ewald分离参数，控制实空间和倒空间计算的平衡，典型值为0.2-0.3 Å^-1
 * @param meshSize PME计算的网格尺寸，数组形式[nx,ny,nz]，通常与盒子尺寸成比例
 * @param pair_cutoff 配对相互作用的截断距离，通常小于或等于PME的实空间截断
 * @param pairGridSize 预计算电势的网格尺寸，数组形式[nx,ny,nz]，决定插值精度
 * @param splineOrder B样条插值的阶数，通常为4(三次B样条)，影响精度和计算速度
 * @param tolerance 计算精度的容差，用于优化参数选择
 * 
 * 使用场景：
 * - 在开始模拟前调用此函数进行初始设置
 * - 当模拟条件(如盒子大小、精度要求)变化时需重新配置
 * - 在每次新的蒙特卡洛模拟开始前设置
 */
void setPGPParameters(double alpha, const int meshSize[3], double pair_cutoff, 
                        const int pairGridSize[3], int splineOrder, double tolerance);

/**
 * @brief 预计算系统中固定部分的网格电势
 * 
 * 这是Precomputed Grid-Potential Particle Mesh Ewald算法的核心函数之一，
 * 负责预计算系统中固定部分的静电势场。它将固定部分的电荷分布到网格上，
 * 通过FFT变换计算电势，并存储结果供后续能量计算使用。
 * 预计算步骤只需在系统固定部分发生变化时执行一次，大大提高了蒙特卡洛模拟的效率。
 * 
 * 算法原理：
 * 1. 将固定部分的点电荷通过B样条插值分布到网格上
 * 2. 对电荷网格执行正向FFT变换到倒空间
 * 3. 在倒空间应用Ewald因子进行长程修正
 * 4. 执行反向FFT获得实空间中的电势分布
 * 5. 将结果存储在预计算的电势网格中
 * 
 * 计算复杂度：
 * - 电荷分配: O(N*p^3)，其中N为原子数，p为B样条阶数
 * - FFT: O(M*log(M))，其中M为网格点总数(nx*ny*nz)
 * - 应用Ewald因子: O(M)
 * 
 * @param state 系统状态，包含原子坐标、电荷和盒子信息
 * @param fixed_only 是否只处理系统中的固定部分(true)，还是处理所有部分(false)
 * 
 * 使用场景：
 * - 系统初始化时预计算固定部分电势
 * - 在GCMC模拟中，固定部分(如蛋白质)的电势场可预先计算
 * - 当固定部分构型变化时，需要重新调用此函数更新电势场
 */
void precomputeGridPotential(model::MCState& state, bool fixed_only = true);

/**
 * @brief 通过插值计算移动分子的能量
 * 
 * 该函数是Precomputed Grid-Potential Particle Mesh Ewald算法的另一个核心函数，
 * 用于在预计算的电势场中快速评估移动分子的能量。利用预计算的电势网格，
 * 通过B样条插值方法高效计算移动分子在该电势场中的能量，避免了直接计算分子间相互作用，
 * 大大加速了蒙特卡洛模拟中的能量评估。
 * 
 * 算法原理：
 * 1. 遍历所有标记为移动的残基和原子
 * 2. 对每个带电原子，通过B样条插值从预计算的电势网格获取其位置的电势值
 * 3. 将电势值乘以原子电荷并累加得到总能量
 * 
 * 计算复杂度：
 * - O(M*p^3)，其中M为移动原子数，p为B样条阶数
 * - 相比传统方法O(M*N)大幅降低，N为固定原子数(通常N >> M)
 * 
 * @param state 系统状态，包含移动分子的信息及预计算的电势网格
 * @param energy 输出参数，存储计算得到的能量值
 * 
 * 使用场景：
 * - GCMC模拟中分子插入/删除的能量评估
 * - CBMC模拟中不同构型的能量比较
 * - MC移动试探中评估新构型的能量变化
 */
void interpolateMoleculeEnergy(model::MCState& state, double& energy);

/**
 * @brief 通过插值计算移动分子的能量，并返回计算结果
 * 
 * 这是interpolateMoleculeEnergy的包装函数，直接返回计算得到的能量值，
 * 方便Python调用和测试。该函数内部创建能量变量并调用原始的interpolateMoleculeEnergy函数。
 * 
 * @param state 系统状态，包含移动分子的信息及预计算的电势网格
 * @return 计算得到的能量值
 */
double calculateMoleculeEnergy(model::MCState& state);

} // namespace cpu
} // namespace platform
} // namespace pygcmc 