#include "PMEComposite.hpp"
#include "PMEGrid.hpp"
#include "PMERecip.hpp"
#include "PMEReal.hpp"
#include "PMESelf.hpp"
#include "../common/EnergyDirectCore.hpp"
#include "platform/platform.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

// Use COULOMB constant from energyCommon.hpp

// Global backup variables for FFT debugging (referenced from original code)
std::vector<std::complex<double>> fftGridBackup;

// Implementation of PMEComposite static methods

void PMEComposite::initialize(double cutoff, const double box[3], 
                            double tolerance, double alpha,
                            const int* meshSize, int splineOrder) {
    // Set box dimensions - ensure B-spline initialization uses correct volume
    pme_params.setBox(box);
    
    // If alpha is not specified, calculate the optimal value
    if (alpha <= 0.0) {
        // Auto-adjust parameters includes calling initializeTables and initializeBsplines
        autoAdjustPMEParameters(tolerance, cutoff, box);
    } else {
        // Use the specified parameters
        int mSize[3] = {64, 64, 64}; // Default value
        
        // If mesh size is provided, use it
        if (meshSize != nullptr) {
            mSize[0] = meshSize[0];
            mSize[1] = meshSize[1];
            mSize[2] = meshSize[2];
        }
        
        // Set parameters and initialize tables
        setPMEParameters(alpha, mSize, splineOrder, tolerance);
        pme_params.initializeTables(cutoff);
        pme_params.initializeBsplines();
    }
    
    // Log the final parameters for debugging
    platform::log(LogLevel::INFO, "PME parameters initialized: alpha = ", pme_params.alpha,
                 ", mesh size = [", pme_params.meshSize[0], ",", pme_params.meshSize[1], ",", pme_params.meshSize[2], "]",
                 ", spline order = ", pme_params.splineOrder,
                 ", initialized = ", pme_params.initialized);
}

void PMEComposite::computeSystemEnergy(model::MCState& state) {
    if (!pme_params.initialized) {
        throw std::runtime_error("PME parameters not initialized. Call initializePMEParameters() first.");
    }
    
    // Calculate energy for the entire system
    // 1. First calculate real space part - needs to be multiplied by COULOMB factor
    computeRealSpacePME(state, false, true);
    
    // 2. Then calculate reciprocal space part - computeReciprocalPME already includes COULOMB factor
    state.ewald_energy.reciprocal = computeReciprocalPME(state);
    
    // 3. Finally calculate self energy part - computeSelfEnergyPME already includes COULOMB factor
    state.ewald_energy.self = computeSelfEnergyPME(state, false);
    
    // 4. Calculate VDW energy - consistent with Ewald, use Direct method
    computeSystemVdwEnergyDirect(state, true, true);
    
    // Only multiply real space energy by COULOMB coefficient
    state.ewald_energy.real_space *= COULOMB;
    
    // Apply COULOMB constant to energies in residues - consistent with Ewald
    for(auto& residue : state.residues) {
        if(residue.active) {
            residue.energy_elec *= COULOMB;
        }
    }
    
    // Calculate total energy - consistent with Ewald, get total energy from residues (including VDW and real-space electrostatics)
    double residue_total = 0.0;
    for (const auto& residue : state.residues) {
        if (residue.active) {
            residue_total += residue.energy_vdw + residue.energy_elec;
        }
    }
    
    // Total energy = residue total energy (including VDW and real-space electrostatics) + Reciprocal + Self
    state.ewald_energy.total = residue_total + 
                             state.ewald_energy.reciprocal + 
                             state.ewald_energy.self;
    
    platform::log(LogLevel::INFO, "PME system energy components: real_space=", state.ewald_energy.real_space,
                 " reciprocal=", state.ewald_energy.reciprocal,
                 " self=", state.ewald_energy.self,
                 " residue_total=", residue_total,
                 " total=", state.ewald_energy.total);
}

void PMEComposite::computeMovementEnergy(model::MCState& state) {
    if (!pme_params.initialized) {
        throw std::runtime_error("PME parameters not initialized. Call initializePMEParameters() first.");
    }
    
    // Only calculate for residues that moved
    computeRealSpacePME(state, true, true);
    state.ewald_energy.reciprocal = computeReciprocalPME(state);
    state.ewald_energy.self = computeSelfEnergyPME(state, true);
    
    // Add VDW energy calculation for movement residues - consistent with Ewald
    computeSystemVdwEnergyDirect(state, true, true);
    
    // Apply Coulomb factor only to real-space component
    // Note: computeReciprocalPME and computeSelfEnergyPME already include COULOMB factor
    state.ewald_energy.real_space *= COULOMB;
    
    // Apply COULOMB constant to energies in movement residues - consistent with Ewald
    for(const auto& movementInfo : state.movementResidues) {
        for(int i = movementInfo.startIndex;
            i < movementInfo.startIndex + movementInfo.activeCount; i++) {
            if(state.residues[i].active) {
                state.residues[i].energy_elec *= COULOMB;
            }
        }
    }
    
    // Calculate total energy for movement residues - consistent with Ewald
    double residue_total = 0.0;
    for(const auto& movementInfo : state.movementResidues) {
        for(int i = movementInfo.startIndex;
            i < movementInfo.startIndex + movementInfo.activeCount; i++) {
            if(state.residues[i].active) {
                residue_total += state.residues[i].energy_vdw + state.residues[i].energy_elec;
            }
        }
    }
    
    // Total energy = residue total energy (including VDW and real-space electrostatics) + Reciprocal + Self
    state.ewald_energy.total = residue_total + 
                             state.ewald_energy.reciprocal + 
                             state.ewald_energy.self;
    
    platform::log(LogLevel::INFO, "PME movement energy components: real_space=", state.ewald_energy.real_space,
                 " reciprocal=", state.ewald_energy.reciprocal,
                 " self=", state.ewald_energy.self,
                 " residue_total=", residue_total,
                 " total=", state.ewald_energy.total);
}

bool PMEComposite::validateSetup(const model::MCState& /* state */) {
    return pme_params.initialized;
}

void PMEComposite::getEnergyBreakdown(const model::MCState& state,
                                    double& realSpace,
                                    double& reciprocal,
                                    double& selfEnergy,
                                    double& vdw,
                                    double& total) {
    realSpace = state.ewald_energy.real_space;
    reciprocal = state.ewald_energy.reciprocal;
    selfEnergy = state.ewald_energy.self;
    
    // Calculate VdW total
    vdw = 0.0;
    for (const auto& residue : state.residues) {
        if (residue.active) {
            vdw += residue.energy_vdw;
        }
    }
    
    total = state.ewald_energy.total;
}

bool PMEComposite::isInitialized() {
    return pme_params.initialized;
}

void PMEComposite::reset() {
    // Clear grid if allocated
    if (!pme_params.pmeGrid.empty()) {
        std::fill(pme_params.pmeGrid.begin(), pme_params.pmeGrid.end(), std::complex<double>(0.0, 0.0));
    }
}

void PMEComposite::setDebugMode(bool enable) {
    platform::log(LogLevel::INFO, "PME debug mode ", (enable ? "enabled" : "disabled"));
}

// Original main function implementations (not inline to avoid redefinition)

// Note: computeSystemEnergyPME and computeMovementEnergyPME are implemented 
// as inline functions in PMEComposite.hpp to avoid redefinition errors

// Note: setPMEParameters and autoAdjustPMEParameters are implemented in PMECore.cpp

// Note: setEnergyDebugOutput is implemented as inline function in energyCommon.hpp

// <agent-hook:composite_implementation>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 