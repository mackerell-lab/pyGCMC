// src/platform/cpu/energyEwald.hpp

#pragma once

#include "model/montecarlo.hpp"
#include "platform/platform.hpp"
#include "energyCommon.hpp"
#include <vector>
#include <complex>

namespace pygcmc {
namespace platform {
namespace cpu {

// Ewald计算的常量
static const int NUM_TABLE_POINTS = 20000;  // 从2048增加到20000以提高精度
static const double TWO_OVER_SQRT_PI = 2.0/std::sqrt(M_PI);

// Ewald参数结构体
struct EwaldParams {
    double alpha{1.0};     // 改为double以提高精度
    int kmax[3]{15,15,15}; // 从6,6,6增加到15,15,15以提高收敛性
    double tolerance{1e-5f};
    bool initialized{false};
    double cutoff{0.0};    // 改为double
    
    // 查找表，用于优化计算
    std::vector<double> erfcTable;      // 改为double
    std::vector<double> ewaldScaleTable;
    double ewaldDX;                     // 改为double
    double ewaldDXInv;
    double erfcDXInv;
    
    // 倒空间优化的exp(ikr)表
    std::vector<std::complex<double>> expIkrTable;  // 改为double
    std::vector<std::complex<double>> expIkrXY;
    int maxK;
    
    // 表格管理方法
    void initializeTables(double cutoff);
    void initializeExpIkrTable(int numAtoms);
    double erfcApprox(double r) const;
    double ewaldScaleApprox(double r) const;
    
    // 误差估计方法
    double estimateRealSpaceError() const {
        return std::erfc(alpha * cutoff);
    }
    
    double estimateReciprocalSpaceError(const double box[3]) const {
        double minBoxSize = std::min(box[0], std::min(box[1], box[2]));
        double minKmax = std::min(kmax[0], std::min(kmax[1], kmax[2]));
        double error = minKmax * std::sqrt(alpha * minBoxSize) / 20.0;
        error *= std::exp(-(M_PI * minKmax / (alpha * minBoxSize)) * 
                         (M_PI * minKmax / (alpha * minBoxSize)));
        return error;
    }
    
    double estimateTotalError(const double box[3]) const {
        return estimateRealSpaceError() + estimateReciprocalSpaceError(box);
    }
};

// 全局Ewald参数
extern EwaldParams ewald_params;

// 函数声明
void setEwaldParameters(double alpha, const int kmax[3], double tolerance = 1e-5);
void autoAdjustParameters(double error_tolerance, double cutoff_distance, const double box[3]);

/**
 * @brief 初始化Ewald参数
 * 
 * @param cutoff 截断距离
 * @param box 盒子尺寸
 * @param alpha Ewald分离参数（如果<=0则自动计算）
 * @param tolerance 精度控制参数
 */
inline void initializeEwaldParameters(double cutoff, const double box[3], 
                                     double alpha = 0.0, double tolerance = 1e-5) {
    // 如果未指定alpha，自动计算最优值
    if (alpha <= 0.0) {
        autoAdjustParameters(tolerance, cutoff, box);
    } else {
        // 使用用户指定的alpha值
        int kmax[3] = {15, 15, 15}; // 默认值
        setEwaldParameters(alpha, kmax, tolerance);
        ewald_params.initializeTables(cutoff);
    }
}

/**
 * @brief 使用Ewald方法计算系统能量
 * 
 * @param state MC状态
 */
void computeSystemEnergyEwald(model::MCState& state);

/**
 * @brief 使用Ewald方法计算运动残基的能量
 * 
 * @param state MC状态
 */
void computeMovementEnergyEwald(model::MCState& state);

// 内部计算函数声明
std::pair<double, double> calcPairEnergyEwald(double r2, double sigma, double eps, double q1, double q2);
double computeReciprocalEnergy(model::MCState& state, bool movement_only);
double computeSelfEnergy(model::MCState& state, bool movement_only);
void computeRealSpaceEwald(model::MCState& state, bool movement_only, bool store_in_residues = true);

} // namespace cpu
} // namespace platform
} // namespace pygcmc 