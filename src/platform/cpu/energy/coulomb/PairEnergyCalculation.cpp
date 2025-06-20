#include "PairEnergyCalculation.hpp"
#include "../common/EnergyUtils.hpp"
#include "platform/platform.hpp"
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace coulomb {

std::pair<double, double> calcPairEnergy(
    double r2, double sigma, double eps, double q1, double q2, 
    const model::MCInfo& info,
    bool calc_coulomb) {
    
    if (getEnergyDebugOutput()) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(6);
        ss << "\n=== calcPairEnergy called ===";
        ss << "\nInput parameters:"
           << "\n  Distance² = " << r2 << " nm²"
           << "\n  Sigma = " << sigma << " nm"
           << "\n  Epsilon = " << eps << " kJ/mol"
           << "\n  q1 = " << q1 << " e"
           << "\n  q2 = " << q2 << " e";
        platform::log(LogLevel::DEBUG, ss.str());
    }

    // Calculate LJ energy using the new LJ module
    double vdw_energy = lj::calcLJEnergyWithSwitching(r2, sigma, eps, info);
    
    double r = std::sqrt(r2);
    
    if (getEnergyDebugOutput()) {
        platform::log(LogLevel::DEBUG, "Distance r = ", r, " nm");
    }
    
    // Calculate Coulomb energy only if requested
    double elec_energy = 0.0;
    if (calc_coulomb) {
        elec_energy = calcCoulombEnergy(r, q1, q2);
    }

    if (getEnergyDebugOutput()) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(6);
        ss << "\nEnergy calculation details:";
        ss << "\n  COULOMB constant = " << COULOMB;
        ss << "\n  q1*q2 = " << (q1 * q2);
        ss << "\nInitial energies:";
        ss << "\n  VDW energy = " << vdw_energy << " kJ/mol";
        ss << "\n  Electrostatic energy = " << elec_energy << " kJ/mol";
        platform::log(LogLevel::DEBUG, ss.str());
    }
    
    // Apply energy capping for numerical stability
    if (getEnergyDebugOutput() && (std::abs(vdw_energy) > MAX_SAFE_ENERGY || std::abs(elec_energy) > MAX_SAFE_ENERGY)) {
        std::stringstream ss;
        ss << "\nEnergy capping applied:";
        ss << "\n  Original VDW energy = " << vdw_energy << " kJ/mol";
        ss << "\n  Original Elec energy = " << elec_energy << " kJ/mol";
        platform::log(LogLevel::DEBUG, ss.str());
    }
    
    // LJ energy is already capped in the LJ module, only cap elec_energy here
    elec_energy = std::min(elec_energy, static_cast<double>(MAX_SAFE_ENERGY));
    elec_energy = std::max(elec_energy, -static_cast<double>(MAX_SAFE_ENERGY));
    
    if (getEnergyDebugOutput() && (std::abs(vdw_energy) > MAX_SAFE_ENERGY || std::abs(elec_energy) > MAX_SAFE_ENERGY)) {
        std::stringstream ss;
        ss << "\nAfter individual capping:";
        ss << "\n  Capped VDW energy = " << vdw_energy << " kJ/mol";
        ss << "\n  Capped Elec energy = " << elec_energy << " kJ/mol";
        platform::log(LogLevel::DEBUG, ss.str());
    }
    
    // Cap total energy
    double total_energy = vdw_energy + elec_energy;
    double original_total = total_energy;
    float max_safe = MAX_SAFE_ENERGY;  // Use float version of safe value
    
    if (total_energy > max_safe) {
        double scale = max_safe / total_energy;
        vdw_energy *= scale;
        elec_energy *= scale;
        if (getEnergyDebugOutput()) {
            std::stringstream ss;
            ss << "\nTotal energy exceeded MAX_SAFE_ENERGY:";
            ss << "\n  Original total = " << original_total << " kJ/mol";
            ss << "\n  Scale factor = " << scale;
            ss << "\n  Final VDW = " << vdw_energy << " kJ/mol";
            ss << "\n  Final Elec = " << elec_energy << " kJ/mol";
            ss << "\n  Final total = " << (vdw_energy + elec_energy) << " kJ/mol";
            platform::log(LogLevel::DEBUG, ss.str());
        }
    } else if (total_energy < -max_safe) {
        double scale = -max_safe / total_energy;
        vdw_energy *= scale;
        elec_energy *= scale;
        if (getEnergyDebugOutput()) {
            std::stringstream ss;
            ss << "\nTotal energy below -MAX_SAFE_ENERGY:";
            ss << "\n  Original total = " << original_total << " kJ/mol";
            ss << "\n  Scale factor = " << scale;
            ss << "\n  Final VDW = " << vdw_energy << " kJ/mol";
            ss << "\n  Final Elec = " << elec_energy << " kJ/mol";
            ss << "\n  Final total = " << (vdw_energy + elec_energy) << " kJ/mol";
            platform::log(LogLevel::DEBUG, ss.str());
        }
    }

    if (getEnergyDebugOutput()) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(6);
        ss << "\n=== calcPairEnergy returning ===";
        ss << "\n  Final VDW energy = " << vdw_energy << " kJ/mol";
        ss << "\n  Final Electrostatic energy = " << elec_energy << " kJ/mol";
        ss << "\n  Total energy = " << (vdw_energy + elec_energy) << " kJ/mol";
        platform::log(LogLevel::DEBUG, ss.str());
    }
    
    return {vdw_energy, elec_energy};
}

} // namespace coulomb
} // namespace cpu
} // namespace platform
} // namespace pygcmc 