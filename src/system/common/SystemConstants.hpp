#pragma once
#ifndef PYGCMC_SYSTEM_COMMON_SYSTEMCONSTANTS_HPP
#define PYGCMC_SYSTEM_COMMON_SYSTEMCONSTANTS_HPP

namespace pygcmc {
namespace system {
namespace common {

/**
 * @brief Physical and unit conversion constants for the system module
 */
class SystemConstants {
public:
    // Unit conversions
    static constexpr float ANGSTROM_TO_NM = 0.1f;
    static constexpr float NM_TO_ANGSTROM = 10.0f;
    static constexpr float KCAL_TO_KJ = 4.184f;
    static constexpr float KJ_TO_KCAL = 1.0f / 4.184f;

    // Physical constants
    static constexpr float BOLTZMANN_CONSTANT = 8.314e-3f; // kJ/(mol·K)
    static constexpr float AVOGADRO_NUMBER = 6.02214076e23f;

    // Numerical constants
    static constexpr float PI = 3.14159265358979323846f;
    static constexpr float TWO_PI = 2.0f * PI;
    static constexpr float SQRT_PI = 1.77245385090551602730f;

    // Cutoff and tolerance defaults
    static constexpr float DEFAULT_CUTOFF = 12.0f;
    static constexpr float DEFAULT_TOLERANCE = 1e-6f;
    static constexpr float DEFAULT_PAIRLIST_BUFFER = 3.0f;

    // Monte Carlo defaults
    static constexpr float DEFAULT_TEMPERATURE = 298.15f;
    static constexpr int DEFAULT_MC_TRIALS = 1000;
    static constexpr float DEFAULT_TRANSLATION_DISTANCE = 1.0f;
};

} // namespace common
} // namespace system
} // namespace pygcmc

#endif // PYGCMC_SYSTEM_COMMON_SYSTEMCONSTANTS_HPP
