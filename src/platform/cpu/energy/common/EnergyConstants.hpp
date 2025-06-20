#pragma once

#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

// Energy calculation constants
extern const float COULOMB;           // Coulomb constant in GROMACS MD units [kJ·nm/mol/e²]
extern const float MIN_SAFE_DISTANCE; // Minimum allowed distance (1% of sigma)
extern const float MAX_SAFE_ENERGY;   // Maximum allowed energy per interaction

// Constants for Ewald and PME calculations
static const int NUM_TABLE_POINTS = 20000;  // High precision table size for approximations
static const double TWO_OVER_SQRT_PI = 2.0/std::sqrt(M_PI);  // Constant for Ewald calculations

// LJ energy calculation related constants
static const float LJ_MIN_SAFE_DISTANCE = 0.01f;  // Minimum safe distance (nm)
static const float LJ_MAX_SAFE_ENERGY = 1e6f;     // Maximum energy value (kJ/mol)

} // namespace cpu
} // namespace platform
} // namespace pygcmc 