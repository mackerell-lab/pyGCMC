#pragma once

#include "EnergyConstants.hpp"
#include "EnergyUtils.hpp"
#include "EnergyLJCalculation.hpp"
#include "model/montecarlo.hpp"
#include <utility>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Calculate LJ and Coulomb energy with safety checks
 * 
 * @param r2 Squared distance in nm²
 * @param sigma LJ sigma parameter in nm
 * @param eps LJ epsilon parameter in kJ/mol
 * @param q1 Charge of first atom in e
 * @param q2 Charge of second atom in e
 * @param info MC information containing switching parameters
 * @param calc_coulomb Whether to calculate Coulomb energy
 * @return Pair of energies: {vdw_energy, elec_energy} in kJ/mol
 */
std::pair<double, double> calcPairEnergy(
    double r2, double sigma, double eps, double q1, double q2, 
    const model::MCInfo& info,
    bool calc_coulomb = true);

} // namespace cpu
} // namespace platform
} // namespace pygcmc 