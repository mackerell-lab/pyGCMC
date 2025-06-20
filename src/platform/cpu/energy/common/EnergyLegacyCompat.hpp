#pragma once

// Legacy compatibility layer for the energy refactoring
// This file provides backward compatibility for code that still uses the old interface

#include "../lj/LJSwitching.hpp"
#include "../coulomb/PairEnergyCalculation.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

// Legacy aliases for LJ calculations
using calculateSwitchingFunction = lj::calculateSwitchingFunction;
using calcLJEnergy = lj::calcLJEnergyWithSwitching;
using setSwitchingFunction = lj::setSwitchingFunction;

// Legacy alias for pair energy calculation
using calcPairEnergy = coulomb::calcPairEnergy;

// Legacy templates (forward to new implementations)
template <typename T>
T checkLJDistance(T r2, T min_safe_distance = LJ_MIN_SAFE_DISTANCE) {
    return lj::checkLJDistance(r2, min_safe_distance);
}

template <typename T>
T capLJEnergy(T energy, T max_safe_energy = LJ_MAX_SAFE_ENERGY) {
    return lj::capLJEnergy(energy, max_safe_energy);
}

template <typename T>
T calculateBasicLJEnergy(T r2, T sigma, T eps) {
    return lj::calculateBasicLJEnergy(r2, sigma, eps);
}

template <typename T>
T calculateLJEnergyNoSwitch(
    T r2, 
    T sigma, 
    T eps, 
    T min_safe_distance = LJ_MIN_SAFE_DISTANCE,
    T max_safe_energy = LJ_MAX_SAFE_ENERGY
) {
    return lj::calculateLJEnergyNoSwitch(r2, sigma, eps, min_safe_distance, max_safe_energy);
}

template <typename T>
T calculateLJEnergy(
    T r2, 
    T sigma, 
    T eps, 
    const model::MCInfo& info,
    T min_safe_distance = T(LJ_MIN_SAFE_DISTANCE),
    T max_safe_energy = T(LJ_MAX_SAFE_ENERGY)
) {
    return lj::calculateLJEnergyWithSwitching(r2, sigma, eps, info, min_safe_distance, max_safe_energy);
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc 