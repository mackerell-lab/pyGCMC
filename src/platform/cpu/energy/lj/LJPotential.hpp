#pragma once

#include "../common/EnergyConstants.hpp"
#include "model/ModelModule.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace lj {

/**
 * @brief Check squared distance, ensure it's not smaller than the square of minimum safe distance
 *
 * @param r2 Squared distance (nm²)
 * @param min_safe_distance Minimum safe distance (nm)
 * @return Safe squared distance (nm²)
 * 
 * @note This function is a template function that must be defined in the header file so the compiler can generate specialized versions at various call points
 */
template <typename T>
T checkLJDistance(T r2, T min_safe_distance = LJ_MIN_SAFE_DISTANCE) {
    const T min_r2 = min_safe_distance * min_safe_distance;
    return (r2 < min_r2) ? min_r2 : r2;
}

/**
 * @brief Limit LJ energy within safe range
 * 
 * @param energy LJ energy (kJ/mol)
 * @param max_safe_energy Maximum safe energy (kJ/mol)
 * @return Limited LJ energy (kJ/mol)
 * 
 * @note This function is a template function that must be defined in the header file so the compiler can generate specialized versions at various call points
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
 * @note This function is a template function that must be defined in the header file so the compiler can generate specialized versions at various call points
 */
template <typename T>
T calculateBasicLJEnergy(T r2, T sigma, T eps) {
    // Use the same implementation as in energyCommon.hpp
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
 * @note This function is a template function that must be defined in the header file so the compiler can generate specialized versions at various call points
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
 * @brief Convenient double version call, using constants as safety parameters
 * 
 * @param r2 Squared distance (nm²)
 * @param sigma LJ sigma parameter (nm)
 * @param eps LJ epsilon parameter (kJ/mol)
 * @return LJ energy (kJ/mol)
 */
double calcLJEnergyBasic(double r2, double sigma, double eps);

} // namespace lj
} // namespace cpu
} // namespace platform
} // namespace pygcmc 