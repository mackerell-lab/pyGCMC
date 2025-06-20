#include "EwaldSelf.hpp"
#include "platform/platform.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Calculate self-energy correction for Ewald summation
 */
double computeSelfEnergy(model::MCState& state, bool movement_only) {
    double self_energy = 0.0;
    
    if(movement_only) {
        // Calculate self-energy only for moving residues
        for(const auto& movementInfo : state.movementResidues) {
            for(int i = movementInfo.startIndex; 
                i < movementInfo.startIndex + movementInfo.activeCount; i++) {
                if(!state.residues[i].active) continue;
                
                for(int j = state.residues[i].atomStart;
                    j < state.residues[i].atomStart + state.residues[i].atomCount; j++) {
                    double charge = state.atoms[j].charge;
                    self_energy += charge * charge;
                }
            }
        }
    } else {
        // Calculate self-energy for all active atoms
        for(int i = 0; i < state.activeAtomCount; i++) {
            double charge = state.atoms[i].charge;
            self_energy += charge * charge;
        }
    }
    
    // Apply the self-energy formula: -COULOMB * alpha / sqrt(π) * Σ q_i²
    self_energy = -COULOMB * ewald_params.alpha / SQRT_PI * self_energy;
    
    return self_energy;
}

/**
 * @brief Calculate self-energy for a specific set of charges
 */
double calculateSelfEnergyForCharges(const std::vector<double>& charges, double alpha) {
    double sum_q2 = 0.0;
    for(double charge : charges) {
        sum_q2 += charge * charge;
    }
    
    return -COULOMB * alpha / std::sqrt(M_PI) * sum_q2;
}

/**
 * @brief Get self-energy breakdown by residue
 */
void getSelfEnergyBreakdown(const model::MCState& state, 
                           std::vector<double>& selfEnergies) {
    selfEnergies.clear();
    selfEnergies.reserve(state.activeResidueCount);
    
    const double prefactor = -COULOMB * ewald_params.alpha / SQRT_PI;
    
    for(int i = 0; i < state.activeResidueCount; i++) {
        if(!state.residues[i].active) {
            selfEnergies.push_back(0.0);
            continue;
        }
        
        double residue_self_energy = 0.0;
        
        // Sum squared charges for this residue
        for(int j = state.residues[i].atomStart;
            j < state.residues[i].atomStart + state.residues[i].atomCount; j++) {
            double charge = state.atoms[j].charge;
            residue_self_energy += charge * charge;
        }
        
        // Apply self-energy formula
        residue_self_energy *= prefactor;
        selfEnergies.push_back(residue_self_energy);
    }
}

/**
 * @brief Validate self-energy calculation parameters
 */
bool validateSelfEnergyParameters(const model::MCState& state) {
    // Check if Ewald parameters are initialized
    if (!ewald_params.initialized) {
        platform::log(LogLevel::ERROR, "Ewald parameters not initialized for self-energy calculation");
        return false;
    }
    
    // Check alpha parameter
    if (ewald_params.alpha <= 0.0) {
        platform::log(LogLevel::ERROR, "Invalid alpha parameter for self-energy calculation: ", ewald_params.alpha);
        return false;
    }
    
    // Check if there are any charges
    bool hasCharges = false;
    for(int i = 0; i < state.activeAtomCount; i++) {
        if (std::abs(state.atoms[i].charge) > 1e-10) {
            hasCharges = true;
            break;
        }
    }
    
    if (!hasCharges) {
        platform::log(LogLevel::INFO, "No charges found in system - self-energy will be zero");
    }
    
    return true;
}

// <agent-hook:ewald_self_impl>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 