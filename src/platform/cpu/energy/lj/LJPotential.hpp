#pragma once

#include "../common/EnergyConstants.hpp"
#include "model/montecarlo.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace lj {

/**
 * @brief 检查距离平方，确保不小于最小安全距离的平方
 *
 * @param r2 距离平方 (nm²)
 * @param min_safe_distance 最小安全距离 (nm)
 * @return 安全的距离平方 (nm²)
 * 
 * @note 此函数是模板函数，必须在头文件中定义，以便编译器能够在各种调用点生成对应的特化版本
 */
template <typename T>
T checkLJDistance(T r2, T min_safe_distance = LJ_MIN_SAFE_DISTANCE) {
    const T min_r2 = min_safe_distance * min_safe_distance;
    return (r2 < min_r2) ? min_r2 : r2;
}

/**
 * @brief 限制LJ能量在安全范围内
 * 
 * @param energy LJ能量 (kJ/mol)
 * @param max_safe_energy 最大安全能量 (kJ/mol)
 * @return 限制后的LJ能量 (kJ/mol)
 * 
 * @note 此函数是模板函数，必须在头文件中定义，以便编译器能够在各种调用点生成对应的特化版本
 */
template <typename T>
T capLJEnergy(T energy, T max_safe_energy = LJ_MAX_SAFE_ENERGY) {
    return std::min(std::max(energy, -max_safe_energy), max_safe_energy);
}

/**
 * @brief Calculate basic Lennard-Jones energy
 * 
 * @param r2 Squared distance in nm²
 * @param sigma LJ sigma parameter in nm
 * @param eps LJ epsilon parameter in kJ/mol
 * @return LJ energy in kJ/mol
 * 
 * @note 此函数是模板函数，必须在头文件中定义，以便编译器能够在各种调用点生成对应的特化版本
 */
template <typename T>
T calculateBasicLJEnergy(T r2, T sigma, T eps) {
    // 使用与energyCommon.hpp中相同的实现方式
    T sigma_r2 = (sigma * sigma) / r2;  // (σ/r)²
    T sigma_r6 = sigma_r2 * sigma_r2 * sigma_r2;  // (σ/r)⁶
    T sigma_r12 = sigma_r6 * sigma_r6;  // (σ/r)¹²
    T vdw_energy = 4.0 * eps * (sigma_r12 - sigma_r6);  // kJ/mol
    
    return vdw_energy;
}

/**
 * @brief Calculate LJ energy without switching function, but with safety checks
 * 
 * @param r2 Squared distance in nm²
 * @param sigma LJ sigma parameter in nm
 * @param eps LJ epsilon parameter in kJ/mol
 * @param min_safe_distance Minimum safe distance to prevent numerical instability
 * @param max_safe_energy Maximum allowed energy value
 * @return LJ energy in kJ/mol
 * 
 * @note 此函数是模板函数，必须在头文件中定义，以便编译器能够在各种调用点生成对应的特化版本
 */
template <typename T>
T calculateLJEnergyNoSwitch(
    T r2, 
    T sigma, 
    T eps, 
    T min_safe_distance = LJ_MIN_SAFE_DISTANCE,
    T max_safe_energy = LJ_MAX_SAFE_ENERGY
) {
    // Apply minimum safe distance for numerical stability
    r2 = checkLJDistance(r2, min_safe_distance);
    
    // Calculate basic LJ energy
    T vdw_energy = calculateBasicLJEnergy(r2, sigma, eps);
    
    // Apply energy capping for numerical stability
    vdw_energy = capLJEnergy(vdw_energy, max_safe_energy);
    
    return vdw_energy;
}

/**
 * @brief 简便的double版本调用，使用常量作为安全参数
 * 
 * @param r2 距离平方 (nm²)
 * @param sigma LJ sigma参数 (nm)
 * @param eps LJ epsilon参数 (kJ/mol)
 * @return LJ能量 (kJ/mol)
 */
double calcLJEnergyBasic(double r2, double sigma, double eps);

} // namespace lj
} // namespace cpu
} // namespace platform
} // namespace pygcmc 