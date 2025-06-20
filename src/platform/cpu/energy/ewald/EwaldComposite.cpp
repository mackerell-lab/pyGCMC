#include "EwaldComposite.hpp"
#include "platform/cpu/energyLJ.hpp"
#include "platform/platform.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Initialize Ewald with automatic parameter optimization
 */
void EwaldComposite::initialize(double cutoff, 
                              const double box[3], 
                              double alpha,
                              double tolerance) {
    if (alpha <= 0.0) {
        autoAdjustParameters(tolerance, cutoff, box);
    } else {
        int kmax[3] = {15, 15, 15}; // Default value
        setEwaldParameters(alpha, kmax, tolerance);
        ewald_params.initializeTables(cutoff);
    }
    
    platform::log(LogLevel::INFO, "EwaldComposite initialized: alpha = ", ewald_params.alpha,
                 ", kmax = [", ewald_params.kmax[0], ",", ewald_params.kmax[1], ",", ewald_params.kmax[2], "]",
                 ", cutoff = ", ewald_params.cutoff);
}

/**
 * @brief Compute total system energy using Ewald summation
 */
void EwaldComposite::computeSystemEnergy(model::MCState& state) {
    // Validate setup
    if (!validateSetup(state)) {
        throw std::runtime_error("Ewald setup validation failed");
    }
    
    // Output Coulomb constant for debugging
    platform::log(LogLevel::INFO, "COULOMB constant in EwaldComposite = ", COULOMB);
    
    // Reset Ewald energy components
    state.ewald_energy.real_space = 0.0;
    state.ewald_energy.reciprocal = 0.0;
    state.ewald_energy.self = 0.0;
    state.ewald_energy.total = 0.0;
    
    // Clear electrostatic energy in residues
    for(auto& residue : state.residues) {
        if(residue.active) {
            residue.energy_elec = 0.0f;
        }
    }
    
    // Real-space part calculation
    computeRealSpaceEwald(state, false, true);
    
    // Apply COULOMB constant to real-space energy
    double real_space_total = state.ewald_energy.real_space * COULOMB;
    state.ewald_energy.real_space = real_space_total;
    
    // Apply COULOMB constant to energies in residues
    for(auto& residue : state.residues) {
        if(residue.active) {
            residue.energy_elec *= COULOMB;
        }
    }
    
    // VDW energy calculation
    computeSystemVdwEnergyDirect(state, true, true);
    
    // Reciprocal space part
    double recip_energy = computeReciprocalEnergy(state, false);
    state.ewald_energy.reciprocal = recip_energy;
    
    // Self-energy correction
    double self_energy = computeSelfEnergy(state, false);
    state.ewald_energy.self = self_energy;
    
    // Calculate total energy
    double residue_total = 0.0;
    for (const auto& residue : state.residues) {
        if (residue.active) {
            residue_total += residue.energy_vdw + residue.energy_elec;
        }
    }
    state.ewald_energy.total = residue_total + state.ewald_energy.reciprocal + state.ewald_energy.self;
    
    // Log energy components
    platform::log(LogLevel::INFO, "\n========== Ewald Energy Components ==========");
    platform::log(LogLevel::INFO, "Real Space Energy:       ", state.ewald_energy.real_space, " kJ/mol");
    platform::log(LogLevel::INFO, "Reciprocal Space Energy: ", state.ewald_energy.reciprocal, " kJ/mol");
    platform::log(LogLevel::INFO, "Self Energy:             ", state.ewald_energy.self, " kJ/mol");
    platform::log(LogLevel::INFO, "Total Ewald Energy:      ", state.ewald_energy.total, " kJ/mol");
    platform::log(LogLevel::INFO, "============================================");
}

/**
 * @brief Compute energy for moving residues only
 */
void EwaldComposite::computeMovementEnergy(model::MCState& state) {
    // Validate setup
    if (!validateSetup(state)) {
        throw std::runtime_error("Ewald setup validation failed");
    }
    
    // Reset Ewald energy components
    state.ewald_energy.real_space = 0.0;
    state.ewald_energy.reciprocal = 0.0;
    state.ewald_energy.self = 0.0;
    state.ewald_energy.total = 0.0;
    
    // Clear electrostatic energy for relevant residues
    for(const auto& movementInfo : state.movementResidues) {
        for(int i = movementInfo.startIndex;
            i < movementInfo.startIndex + movementInfo.activeCount; i++) {
            if(state.residues[i].active) {
                state.residues[i].energy_elec = 0.0f;
            }
        }
    }
    
    // Real-space part for moving residues
    computeRealSpaceEwald(state, true, true);
    
    // Apply COULOMB constant to real-space energy
    double real_space_total = state.ewald_energy.real_space * COULOMB;
    state.ewald_energy.real_space = real_space_total;
    
    // Apply COULOMB constant to energies in residues
    for(const auto& movementInfo : state.movementResidues) {
        for(int i = movementInfo.startIndex;
            i < movementInfo.startIndex + movementInfo.activeCount; i++) {
            if(state.residues[i].active) {
                state.residues[i].energy_elec *= COULOMB;
            }
        }
    }
    
    // VDW energy calculation
    computeSystemVdwEnergyDirect(state, true, true);
    
    // Reciprocal space part for moving residues
    double recip_energy = computeReciprocalEnergy(state, true);
    state.ewald_energy.reciprocal = recip_energy;
    
    // Self-energy correction for moving residues
    double self_energy = computeSelfEnergy(state, true);
    state.ewald_energy.self = self_energy;
    
    // Calculate total energy
    double residue_total = 0.0;
    for(const auto& movementInfo : state.movementResidues) {
        for(int i = movementInfo.startIndex;
            i < movementInfo.startIndex + movementInfo.activeCount; i++) {
            if(state.residues[i].active) {
                residue_total += state.residues[i].energy_vdw + state.residues[i].energy_elec;
            }
        }
    }
    state.ewald_energy.total = residue_total + state.ewald_energy.reciprocal + state.ewald_energy.self;
    
    // Log energy components
    platform::log(LogLevel::INFO, "\n========== Movement Ewald Energy Components ==========");
    platform::log(LogLevel::INFO, "Real Space Energy:       ", state.ewald_energy.real_space, " kJ/mol");
    platform::log(LogLevel::INFO, "Reciprocal Space Energy: ", state.ewald_energy.reciprocal, " kJ/mol");
    platform::log(LogLevel::INFO, "Self Energy:             ", state.ewald_energy.self, " kJ/mol");
    platform::log(LogLevel::INFO, "Total Ewald Energy:      ", state.ewald_energy.total, " kJ/mol");
    platform::log(LogLevel::INFO, "====================================================");
}

/**
 * @brief Validate Ewald setup and parameters
 */
bool EwaldComposite::validateSetup(const model::MCState& state) {
    if (!ewald_params.initialized) {
        platform::log(LogLevel::ERROR, "Ewald parameters not initialized");
        return false;
    }
    
    // Validate individual components
    if (!validateRealSpaceParameters(state)) {
        return false;
    }
    
    if (!validateReciprocalSpaceParameters(state)) {
        return false;
    }
    
    if (!validateSelfEnergyParameters(state)) {
        return false;
    }
    
    return true;
}

/**
 * @brief Get Ewald energy components breakdown
 */
void EwaldComposite::getEnergyBreakdown(const model::MCState& state,
                                      double& realSpace,
                                      double& reciprocal,
                                      double& selfEnergy,
                                      double& vdw,
                                      double& total) {
    realSpace = state.ewald_energy.real_space;
    reciprocal = state.ewald_energy.reciprocal;
    selfEnergy = state.ewald_energy.self;
    
    // Calculate VdW energy from residues
    vdw = 0.0;
    for (const auto& residue : state.residues) {
        if (residue.active) {
            vdw += residue.energy_vdw;
        }
    }
    
    total = state.ewald_energy.total;
}

/**
 * @brief Check if Ewald is properly initialized
 */
bool EwaldComposite::isInitialized() {
    return ewald_params.initialized;
}

/**
 * @brief Reset Ewald state for new calculation
 */
void EwaldComposite::reset() {
    ewald_params.initialized = false;
    ewald_params.erfcTable.clear();
    ewald_params.ewaldScaleTable.clear();
    ewald_params.expIkrTable.clear();
    ewald_params.expIkrXY.clear();
    
    platform::log(LogLevel::INFO, "Ewald state reset");
}

/**
 * @brief Get convergence information for all components
 */
void EwaldComposite::getConvergenceInfo(const model::MCState& state,
                                      double& realSpaceError,
                                      double& reciprocalError,
                                      double& totalError) {
    const double* box = reinterpret_cast<const double*>(state.info.box);
    
    realSpaceError = ewald_params.estimateRealSpaceError();
    reciprocalError = ewald_params.estimateReciprocalSpaceError(box);
    totalError = ewald_params.estimateTotalError(box);
    
    platform::log(LogLevel::INFO, 
        "Ewald convergence info: real space error = ", realSpaceError,
        ", reciprocal space error = ", reciprocalError,
        ", total error = ", totalError);
}

/**
 * @brief Validate system properties
 */
void EwaldComposite::validateSystemProperties(const model::MCState& state) {
    // Check charge neutrality
    double totalCharge = 0.0;
    for(int i = 0; i < state.activeAtomCount; i++) {
        totalCharge += state.atoms[i].charge;
    }
    
    if (std::abs(totalCharge) > 1e-6) {
        platform::log(LogLevel::WARNING, 
            "System charge (", totalCharge, ") deviates from neutrality. ",
            "This may affect Ewald summation accuracy.");
    }
}

/**
 * @brief Check parameter consistency
 */
void EwaldComposite::checkParameterConsistency() {
    if (ewald_params.alpha <= 0.0) {
        throw std::runtime_error("Invalid alpha parameter: must be positive");
    }
    
    for(int i = 0; i < 3; i++) {
        if (ewald_params.kmax[i] <= 0) {
            throw std::runtime_error("Invalid kmax parameter: must be positive");
        }
    }
}

// <agent-hook:ewald_composite_impl>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 