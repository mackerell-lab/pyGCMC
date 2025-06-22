#pragma once

#include "../common/EnergyConstants.hpp"
#include "model/ModelModule.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace coulomb {

/**
 * @brief Calculate basic Coulomb energy
 * 
 * @param r Distance in nm
 * @param q1 Charge of first atom in e
 * @param q2 Charge of second atom in e
 * @return Coulomb energy in kJ/mol
 * 
 * @note This function is a template function that must be defined in the header file so the compiler can generate specialized versions at various call points
 */
template <typename T>
T calculateBasicCoulombEnergy(T r, T q1, T q2) {
    return COULOMB * q1 * q2 / r;  // kJ/mol
}

/**
 * @brief Calculate Coulomb energy with safety checks
 * 
 * @param r Distance in nm
 * @param q1 Charge of first atom in e
 * @param q2 Charge of second atom in e
 * @param max_safe_energy Maximum allowed energy value
 * @return Coulomb energy in kJ/mol
 * 
 * @note This function is a template function that must be defined in the header file so the compiler can generate specialized versions at various call points
 */
template <typename T>
T calculateCoulombEnergy(T r, T q1, T q2, T max_safe_energy = T(MAX_SAFE_ENERGY)) {
    T elec_energy = calculateBasicCoulombEnergy(r, q1, q2);
    
    // Apply energy capping for numerical stability
    elec_energy = std::min(elec_energy, max_safe_energy);
    elec_energy = std::max(elec_energy, -max_safe_energy);
    
    return elec_energy;
}

/**
 * @brief Convenient double version call, using constants as safety parameters
 * 
 * @param r Distance (nm)
 * @param q1 Charge of first atom (e)
 * @param q2 Charge of second atom (e)
 * @return Coulomb energy (kJ/mol)
 */
double calcCoulombEnergy(double r, double q1, double q2);

} // namespace coulomb
} // namespace cpu
} // namespace platform
} // namespace pygcmc 