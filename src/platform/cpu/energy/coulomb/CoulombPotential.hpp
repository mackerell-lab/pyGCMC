#pragma once

#include "../common/EnergyConstants.hpp"
#include "model/montecarlo.hpp"
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
 * @note 此函数是模板函数，必须在头文件中定义，以便编译器能够在各种调用点生成对应的特化版本
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
 * @note 此函数是模板函数，必须在头文件中定义，以便编译器能够在各种调用点生成对应的特化版本
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
 * @brief 简便的double版本调用，使用常量作为安全参数
 * 
 * @param r 距离 (nm)
 * @param q1 第一个原子的电荷 (e)
 * @param q2 第二个原子的电荷 (e)
 * @return 库伦能量 (kJ/mol)
 */
double calcCoulombEnergy(double r, double q1, double q2);

} // namespace coulomb
} // namespace cpu
} // namespace platform
} // namespace pygcmc 