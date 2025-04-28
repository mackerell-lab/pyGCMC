// src/platform/cpu/energyLJ.hpp

#pragma once

#include "model/montecarlo.hpp"
#include "platform/platform.hpp"
#include "energyCommon.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

// LJ能量计算相关常量
static const float LJ_MIN_SAFE_DISTANCE = 0.01f;  // 最小安全距离 (nm)
static const float LJ_MAX_SAFE_ENERGY = 1e6f;     // 最大能量值 (kJ/mol)

/**
 * @brief Calculate CHARMM switching function value
 * 
 * @param r Distance in nm
 * @param info MC information containing switching parameters
 * @return Switching function value (between 0 and 1)
 */
float calculateSwitchingFunction(float r, const model::MCInfo& info);

/**
 * @brief 检查距离平方，确保不小于最小安全距离的平方
 *
 * @param r2 距离平方 (nm²)
 * @param min_safe_distance 最小安全距离 (nm)
 * @return 安全的距离平方 (nm²)
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
    T min_safe_distance = T(LJ_MIN_SAFE_DISTANCE),
    T max_safe_energy = T(LJ_MAX_SAFE_ENERGY)
) {
    // 快速路径：如果不使用switching function，直接调用简化版本
    if (!info.use_switching) {
        return calculateLJEnergyNoSwitch(r2, sigma, eps, min_safe_distance, max_safe_energy);
    }
    
    // Apply minimum safe distance for numerical stability
    r2 = checkLJDistance(r2, min_safe_distance);
    
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
    vdw_energy = capLJEnergy(vdw_energy, max_safe_energy);
    
    return vdw_energy;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc 