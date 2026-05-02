#pragma once

#include "LJPotential.hpp"
#include "model/ModelModule.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace lj {

/**
 * @brief Calculate CHARMM switching function value
 *
 * @param r Distance in nm
 * @param info MC information containing switching parameters
 * @return Switching function value (between 0 and 1)
 */
float calculateSwitchingFunction(float r, const model::MCInfo& info);

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
 *
 * @note This function is a template function that must be defined in the header file so the compiler can generate specialized versions at various call points
 */
template <typename T>
T calculateLJEnergyWithSwitching(
    T r2,
    T sigma,
    T eps,
    const model::MCInfo& info,
    T min_safe_distance = T(LJ_MIN_SAFE_DISTANCE),
    T max_safe_energy = T(LJ_MAX_SAFE_ENERGY)
) {
    // Fast path: if not using switching function, directly call simplified version
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

/**
 * @brief Convenient double version call, using constants as safety parameters
 *
 * @param r2 Squared distance (nm²)
 * @param sigma LJ sigma parameter (nm)
 * @param eps LJ epsilon parameter (kJ/mol)
 * @param info MC information containing switching function parameters
 * @return LJ energy (kJ/mol)
 */
double calcLJEnergyWithSwitching(double r2, double sigma, double eps, const model::MCInfo& info);

/**
 * @brief Configure CHARMM-style switching function
 *
 * @param state MC state to configure
 * @param use_switching Whether to enable switching
 * @param r_on Inner switching radius
 * @param r_off Outer switching radius
 */
void setSwitchingFunction(model::MCState& state, bool use_switching, float r_on, float r_off);

} // namespace lj
} // namespace cpu
} // namespace platform
} // namespace pygcmc
