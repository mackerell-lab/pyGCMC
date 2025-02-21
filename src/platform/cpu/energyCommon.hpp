// src/platform/cpu/energyCommon.hpp

#pragma once

#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

// Constants for energy calculations
extern const float COULOMB;      // Coulomb constant in GROMACS MD units [kJ·nm/mol/e²]
extern const float MIN_SAFE_DISTANCE;  // Minimum allowed distance (1% of sigma)
extern const float MAX_SAFE_ENERGY;    // Maximum allowed energy per interaction

// Common utility functions
inline float capEnergy(float energy) {
    return std::min(std::max(energy, -MAX_SAFE_ENERGY), MAX_SAFE_ENERGY);
}

inline float checkDistance(float r2) {
    const float min_r2 = MIN_SAFE_DISTANCE * MIN_SAFE_DISTANCE;
    return (r2 < min_r2) ? min_r2 : r2;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc 