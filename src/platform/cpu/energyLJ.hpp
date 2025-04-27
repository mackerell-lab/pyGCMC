// src/platform/cpu/energyLJ.hpp

#pragma once

#include "model/montecarlo.hpp"
#include "platform/platform.hpp"
#include "energyCommon.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Calculate CHARMM switching function value
 * 
 * @param r Distance in nm
 * @param info MC information containing switching parameters
 * @return Switching function value (between 0 and 1)
 */
float calculateSwitchingFunction(float r, const model::MCInfo& info);

/**
 * @brief Calculate basic Lennard-Jones energy
 * 
 * @param r2 Squared distance in nm²
 * @param sigma LJ sigma parameter in nm
 * @param eps LJ epsilon parameter in kJ/mol
 * @return LJ energy in kJ/mol
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
 */
template <typename T>
T calculateLJEnergyNoSwitch(
    T r2, 
    T sigma, 
    T eps, 
    T min_safe_distance = 0.01,
    T max_safe_energy = 1e6
) {
    // Apply minimum safe distance for numerical stability
    if (r2 < min_safe_distance * min_safe_distance) {
        r2 = min_safe_distance * min_safe_distance;
    }
    
    // Calculate basic LJ energy
    T vdw_energy = calculateBasicLJEnergy(r2, sigma, eps);
    
    // Apply energy capping for numerical stability
    vdw_energy = std::min(vdw_energy, max_safe_energy);
    vdw_energy = std::max(vdw_energy, -max_safe_energy);
    
    return vdw_energy;
}

/**
 * @brief Calculate LJ energy with optional switching and safety checks
 * 
 * @param r2 Squared distance in nm²
 * @param sigma LJ sigma parameter in nm
 * @param eps LJ epsilon parameter in kJ/mol
 * @param info MC information containing switching parameters
 * @param min_safe_distance Minimum safe distance to prevent numerical instability
 * @param max_safe_energy Maximum allowed energy value
 * @return LJ energy in kJ/mol
 */
template <typename T>
T calculateLJEnergy(
    T r2, 
    T sigma, 
    T eps, 
    const model::MCInfo& info,
    T min_safe_distance,
    T max_safe_energy
) {
    // 快速路径：如果不使用switching function，直接调用简化版本
    if (!info.use_switching) {
        return calculateLJEnergyNoSwitch(r2, sigma, eps, min_safe_distance, max_safe_energy);
    }
    
    // Apply minimum safe distance for numerical stability
    if (r2 < min_safe_distance * min_safe_distance) {
        r2 = min_safe_distance * min_safe_distance;
    }
    
    T r = std::sqrt(r2);
    
    // Calculate basic LJ energy
    T vdw_energy = calculateBasicLJEnergy(r2, sigma, eps);
    
    // Apply CHARMM switching function to LJ energy
    // Apply switching if distance is between r_on and r_off
    if (r > info.r_on && r < info.r_off) {
        T switch_val = static_cast<T>(calculateSwitchingFunction(static_cast<float>(r), info));
        vdw_energy *= switch_val;
    } else if (r >= info.r_off) {
        // Beyond outer cutoff radius, set to zero
        vdw_energy = 0;
    }
    
    // Apply energy capping for numerical stability
    vdw_energy = std::min(vdw_energy, max_safe_energy);
    vdw_energy = std::max(vdw_energy, -max_safe_energy);
    
    return vdw_energy;
}

// 特化的函数重载，用于处理混合类型参数的情况
inline double calculateLJEnergy(
    double r2, 
    double sigma, 
    double eps, 
    const model::MCInfo& info,
    float min_safe_distance,
    float max_safe_energy
) {
    return calculateLJEnergy(r2, sigma, eps, info, 
                           static_cast<double>(min_safe_distance), 
                           static_cast<double>(max_safe_energy));
}

// 特化的函数重载，用于处理混合类型参数的情况
inline float calculateLJEnergy(
    float r2, 
    float sigma, 
    float eps, 
    const model::MCInfo& info,
    double min_safe_distance,
    double max_safe_energy
) {
    return calculateLJEnergy(r2, sigma, eps, info, 
                           static_cast<float>(min_safe_distance), 
                           static_cast<float>(max_safe_energy));
}

// 带默认参数的模板函数
template <typename T>
T calculateLJEnergy(
    T r2, 
    T sigma, 
    T eps, 
    const model::MCInfo& info
) {
    return calculateLJEnergy(r2, sigma, eps, info, T(0.01), T(1e6));
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc 