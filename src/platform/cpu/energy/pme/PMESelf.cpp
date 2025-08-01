#include "PMESelf.hpp"
#include "PMEGlobal.hpp"
#include "PMECore.hpp"
#include "platform/platform.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

// Use COULOMB constant from energyCommon.hpp

/**
 * @brief Compute self-energy term for PME
 * 
 * @param state MC state
 * @param movement_only Whether to compute only for moving atoms
 * @return double Self-energy
 */
double computeSelfEnergyPME(model::MCState& state, bool movement_only) {
    platform::log(LogLevel::INFO, "Computing self energy with alpha = ", pme_params.alpha);
    
    // Self-energy calculation same as Ewald
    double self_energy = 0.0;
    
    // Track sum of squared charges
    double sum_q2 = 0.0;
    int count = 0;
    
    if(movement_only) {
        for(const auto& movementInfo : state.movementResidues) {
            for(int i = movementInfo.startIndex; 
                i < movementInfo.startIndex + movementInfo.activeCount; i++) {
                if(!state.residues[i].active) continue;
                
                for(int j = state.residues[i].atomStart;
                    j < state.residues[i].atomStart + state.residues[i].atomCount; j++) {
                    if (j < state.activeAtomCount) { // Ensure index is valid
                        double charge = state.atoms[j].charge;
                        double q2 = charge * charge;
                        sum_q2 += q2;
                        
                        if (count < 5) {
                            platform::log(LogLevel::INFO, "Atom ", j, 
                                        " charge = ", charge, 
                                        ", q² = ", q2);
                            count++;
                        }
                    }
                }
            }
        }
    }
    else {
        for(int i = 0; i < state.activeAtomCount; i++) {
            // No longer check active member
            double charge = state.atoms[i].charge;
            double q2 = charge * charge;
            sum_q2 += q2;
            
            if (i < 5) {
                platform::log(LogLevel::INFO, "Atom ", i, 
                            " charge = ", charge, 
                            ", q² = ", q2);
            }
        }
    }
    
    platform::log(LogLevel::INFO, "Sum of q² = ", sum_q2);
    
    // Self-energy formula from pme.cpp: -ONE_4PI_EPS0 * alpha / sqrt(M_PI) * sum_q2
    // Ensure we use exactly the same formula, including the COULOMB constant
    double prefactor = -COULOMB * pme_params.alpha / sqrt(M_PI);
    self_energy = prefactor * sum_q2;
    
    platform::log(LogLevel::INFO, "Self energy prefactor = ", prefactor, 
                 ", resulting self energy = ", self_energy);
    
    // Self energy already includes COULOMB constant, no need to multiply elsewhere
    return self_energy;
}

// Note: computeSelfEnergy is implemented in energyEwald.cpp for general use
// PME module provides computeSelfEnergyPME for PME-specific calculations

/**
 * @brief Calculate self-energy for a specific charge
 */
double calculateSelfEnergyForCharge(double charge, double alpha) {
    if (std::abs(charge) < 1e-6) return 0.0;
    
    double prefactor = -COULOMB * alpha / sqrt(M_PI);
    return prefactor * charge * charge;
}

/**
 * @brief Calculate total self-energy for a set of charges
 */
double calculateTotalSelfEnergy(const std::vector<double>& charges, double alpha) {
    double sum_q2 = 0.0;
    
    for (double charge : charges) {
        sum_q2 += charge * charge;
    }
    
    double prefactor = -COULOMB * alpha / sqrt(M_PI);
    return prefactor * sum_q2;
}

/**
 * @brief Validate self-energy calculation parameters
 */
bool validateSelfEnergyParameters(double alpha) {
    if (alpha <= 0.0) {
        platform::log(LogLevel::ERROR, "Invalid alpha parameter for self-energy: ", alpha);
        return false;
    }
    
    return true;
}

/**
 * @brief Get self-energy contribution statistics
 */
void getSelfEnergyStatistics(const model::MCState& state, 
                           double& totalSelfEnergy,
                           double& maxAtomContribution,
                           int& chargedAtomCount) {
    totalSelfEnergy = 0.0;
    maxAtomContribution = 0.0;
    chargedAtomCount = 0;
    
    double prefactor = -COULOMB * pme_params.alpha / sqrt(M_PI);
    
    for (int i = 0; i < state.activeAtomCount; i++) {
        double charge = state.atoms[i].charge;
        if (std::abs(charge) > 1e-6) {
            chargedAtomCount++;
            double contribution = prefactor * charge * charge;
            totalSelfEnergy += contribution;
            maxAtomContribution = std::max(maxAtomContribution, std::abs(contribution));
        }
    }
}

/**
 * @brief Calculate self-energy difference for charge modifications
 */
double calculateSelfEnergyDifference(double oldCharge, double newCharge, double alpha) {
    double prefactor = -COULOMB * alpha / sqrt(M_PI);
    return prefactor * (newCharge * newCharge - oldCharge * oldCharge);
}

// <agent-hook:self_implementation>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 
