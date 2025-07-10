// src/simulation/simulation.cpp
#include "simulation.hpp"
#include "../platform/cpu/energy/EnergyModule.hpp"
#include "../platform/cpu/energy/pme/PMEComposite.hpp"
#include "../platform/cpu/energy/pgp/PGPCore.hpp"
#include "../platform/cpu/energy/pgp/PGPComplete.hpp"
#include <cmath>

namespace pygcmc {
namespace simulation {

void Simulation::computeMovementEnergy(model::MCState& state) {
    if (is_debug_enabled()) {
        log(LogLevel::DEBUG, "Computing nonbonded energy for movement residues");
    }
    platform::cpu::computeMovementEnergy(state, platform::cpu::EnergyMethod::DIRECT, false, false);
}

void Simulation::computeMovementEnergyCutoff(model::MCState& state) {
    if (is_debug_enabled()) {
        log(LogLevel::DEBUG, "Computing nonbonded energy for movement residues with cutoff");
    }
    platform::cpu::computeMovementEnergyCutoff(state);
}

void Simulation::computeSystemEnergy(model::MCState& state) {
    if (is_debug_enabled()) {
        log(LogLevel::DEBUG, "Computing nonbonded energy for all active residues");
    }
    
    platform::cpu::computeSystemEnergy(state, platform::cpu::EnergyMethod::DIRECT, false, false);
    
    // Only log total energy in debug mode
    if (is_debug_enabled()) {
        float total_vdw = 0.0f;
        float total_elec = 0.0f;
        for (int i = 0; i < state.activeResidueCount; ++i) {
            if (state.residues[i].active) {
                total_vdw += state.residues[i].energy_vdw;
                total_elec += state.residues[i].energy_elec;
            }
        }

        total_vdw /= 2.0f;
        total_elec /= 2.0f; 
        
        log(LogLevel::DEBUG, "Total system energy: vdw=", total_vdw, 
            ", elec=", total_elec, 
            ", total=", (total_vdw + total_elec));
    }
}

void Simulation::computeSystemEnergyCutoff(model::MCState& state) {
    if (is_debug_enabled()) {
        log(LogLevel::DEBUG, "Computing cutoff nonbonded energy for all active residues");
    }
    
    platform::cpu::computeSystemEnergyCutoff(state);
    
    // Only log total energy in debug mode
    if (is_debug_enabled()) {
        float total_vdw = 0.0f;
        float total_elec = 0.0f;
        for (int i = 0; i < state.activeResidueCount; ++i) {
            if (state.residues[i].active) {
                total_vdw += state.residues[i].energy_vdw;
                total_elec += state.residues[i].energy_elec;
            }
        }
        
        total_vdw /= 2.0f;
        total_elec /= 2.0f;
        
        log(LogLevel::DEBUG, "Total system energy: vdw=", total_vdw, 
            ", elec=", total_elec, 
            ", total=", (total_vdw + total_elec));
    }
}

void Simulation::computeSystemEnergyPBC(model::MCState& state) {
    // Validate box dimensions before proceeding
    if (state.info.box[0] <= 0.0f || state.info.box[1] <= 0.0f || state.info.box[2] <= 0.0f) {
        throw std::runtime_error("Invalid box dimensions for PBC calculation");
    }
    
    if (is_debug_enabled()) {
        log(LogLevel::DEBUG, "Computing PBC nonbonded energy (no cutoff): box=", 
            state.info.box[0], "x", state.info.box[1], "x", state.info.box[2], " nm");
    }
    
    platform::cpu::computeSystemEnergyPBC(state);
    
    // Only log total energy in debug mode
    if (is_debug_enabled()) {
        float total_vdw = 0.0f;
        float total_elec = 0.0f;
        for (int i = 0; i < state.activeResidueCount; ++i) {
            if (state.residues[i].active) {
                total_vdw += state.residues[i].energy_vdw;
                total_elec += state.residues[i].energy_elec;
            }
        }
        
        total_vdw /= 2.0f;
        total_elec /= 2.0f;
        
        log(LogLevel::DEBUG, "Total system energy with PBC (no cutoff): vdw=", total_vdw, 
            ", elec=", total_elec, 
            ", total=", (total_vdw + total_elec));
    }
}

void Simulation::computeSystemEnergyPBCCutoff(model::MCState& state) {
    if (is_debug_enabled()) {
        log(LogLevel::DEBUG, "Computing nonbonded energy for all active residues with PBC and cutoff");
    }
    platform::cpu::computeSystemEnergyPBCCutoff(state);
}

void Simulation::computeSystemVdwEnergyCutoff(model::MCState& state) {
    if (is_debug_enabled()) {
        log(LogLevel::DEBUG, "Computing VDW energies with cutoff");
    }
    platform::cpu::computeSystemVdwEnergyCutoff(state);
}

void Simulation::setEwaldParameters(float alpha, const int kmax[3], float tolerance) {
    platform::cpu::setEwaldParameters(alpha, kmax, tolerance);
}

void Simulation::initializeEwaldParameters(float cutoff, const float box[3], 
                                        float alpha, float tolerance) {
    // Convert float parameters to double
    double cutoff_d = static_cast<double>(cutoff);
    double box_d[3] = {
        static_cast<double>(box[0]),
        static_cast<double>(box[1]),
        static_cast<double>(box[2])
    };
    double alpha_d = static_cast<double>(alpha);
    double tolerance_d = static_cast<double>(tolerance);
    
    platform::cpu::initializeEwaldParameters(cutoff_d, box_d, alpha_d, tolerance_d);
}

void Simulation::computeSystemEnergyEwald(model::MCState& state) {
    if (is_debug_enabled()) {
        log(LogLevel::DEBUG, "Computing Ewald energy for all active residues");
    }
    platform::cpu::computeSystemEnergyEwald(state);
}

void Simulation::computeMovementEnergyEwald(model::MCState& state) {
    if (is_debug_enabled()) {
        log(LogLevel::DEBUG, "Computing Ewald energy for movement residues");
    }
    platform::cpu::computeMovementEnergyEwald(state);
}

void Simulation::setEnergyDebugOutput(bool enable) {
    // Set both simulation-level and platform-level debug output
    if (enable) {
        set_verbose(true);
        set_log_level(LogLevel::DEBUG);
    }
    platform::cpu::setEnergyDebugOutput(enable);
}

// Implementation of PME-related methods
void Simulation::setPMEParameters(float alpha, const int meshSize[3], int splineOrder, float tolerance) {
    platform::cpu::setPMEParameters(alpha, meshSize, splineOrder, tolerance);
}

void Simulation::initializePMEParameters(float cutoff, const float box[3], 
                                      float alpha, const int* meshSize,
                                      int splineOrder, float tolerance) {
    // Convert float parameters to double
    double cutoff_d = static_cast<double>(cutoff);
    double box_d[3] = {
        static_cast<double>(box[0]),
        static_cast<double>(box[1]),
        static_cast<double>(box[2])
    };
    double alpha_d = static_cast<double>(alpha);
    double tolerance_d = static_cast<double>(tolerance);
    
    // Handle meshSize conversion
    int meshSize_d[3] = {0, 0, 0};
    if (meshSize != nullptr) {
        meshSize_d[0] = meshSize[0];
        meshSize_d[1] = meshSize[1];
        meshSize_d[2] = meshSize[2];
    }
    
    platform::cpu::initializePMEParameters(cutoff_d, box_d, alpha_d, 
                                        meshSize != nullptr ? meshSize_d : nullptr, 
                                        splineOrder, tolerance_d);
}

void Simulation::computeSystemEnergyPME(model::MCState& state) {
    if (is_debug_enabled()) {
        log(LogLevel::DEBUG, "Computing PME energy for all active residues");
    }
    platform::cpu::computeSystemEnergyPME(state);
}

void Simulation::computeMovementEnergyPME(model::MCState& state) {
    if (is_debug_enabled()) {
        log(LogLevel::DEBUG, "Computing PME energy for movement residues");
    }
    platform::cpu::computeMovementEnergyPME(state);
}

void Simulation::clearPMEEngine() {
    if (is_debug_enabled()) {
        log(LogLevel::DEBUG, "Clearing PME engine state");
    }
    platform::cpu::clearPMEState();
}

void Simulation::computeSystemEnergyPMEComplete(model::MCState& state) {
    if (is_debug_enabled()) {
        log(LogLevel::DEBUG, "Computing complete PME system energy with intramolecular interactions");
    }
    platform::cpu::computeSystemEnergyPMEComplete(state);
}

void Simulation::computeSystemEnergyCutoffComplete(model::MCState& state) {
    if (is_debug_enabled()) {
        log(LogLevel::DEBUG, "Computing complete cutoff system energy with intramolecular interactions");
    }
    platform::cpu::computeSystemEnergyCutoffComplete(state);
}

// Implementation of PGP-related methods
void Simulation::setPGPParameters(float alpha, const int meshSize[3], float potential_cutoff, 
                               const int potentialGridSize[3], int splineOrder, float tolerance) {
    // Convert float parameters to double
    double alpha_d = static_cast<double>(alpha);
    double potential_cutoff_d = static_cast<double>(potential_cutoff);
    double tolerance_d = static_cast<double>(tolerance);
    
    // Convert int arrays
    int meshSize_d[3] = {meshSize[0], meshSize[1], meshSize[2]};
    int potentialGridSize_d[3] = {potentialGridSize[0], potentialGridSize[1], potentialGridSize[2]};
    
    // Call the CPU platform implementation
    platform::cpu::setPGPParameters(alpha_d, meshSize_d, potential_cutoff_d, 
                                 potentialGridSize_d, splineOrder, tolerance_d);
}

void Simulation::precomputeGridPotential(model::MCState& state, bool fixed_only) {
    if (is_debug_enabled()) {
        log(LogLevel::DEBUG, "Precomputing grid potential for PGP-PME: fixed_only=", fixed_only);
    }
    platform::cpu::precomputeGridPotential(state, fixed_only);
}

void Simulation::interpolateMoleculeEnergy(model::MCState& state, double& energy) {
    platform::cpu::interpolateMoleculeEnergy(state, energy);
}

double Simulation::calculateMoleculeEnergy(model::MCState& state) {
    if (is_debug_enabled()) {
        log(LogLevel::DEBUG, "Calculating molecule energy using PGP interpolation");
    }
    return platform::cpu::calculateMoleculeEnergy(state);
}

void Simulation::computeSystemEnergyPGP(model::MCState& state) {
    if (is_debug_enabled()) {
        log(LogLevel::DEBUG, "Computing PGP energy for all active residues");
    }
    platform::cpu::computeSystemEnergyPGP(state);
}

void Simulation::computeMovementEnergyPGP(model::MCState& state) {
    if (is_debug_enabled()) {
        log(LogLevel::DEBUG, "Computing PGP energy for movement residues");
    }
    platform::cpu::computeMovementEnergyPGP(state);
}

void Simulation::computeSystemEnergyPGPFixed(model::MCState& state) {
    if (is_debug_enabled()) {
        log(LogLevel::DEBUG, "Computing PGP-Fixed energy for all active residues");
    }
    platform::cpu::computeSystemEnergyPGPFixed(state);
}

void Simulation::computeMovementEnergyPGPFixed(model::MCState& state) {
    if (is_debug_enabled()) {
        log(LogLevel::DEBUG, "Computing PGP-Fixed energy for movement residues");
    }
    platform::cpu::computeMovementEnergyPGPFixed(state);
}

std::pair<double, double> Simulation::getTotalEnergyComponents(const model::MCState& state) {
    double total_elec = 0.0;
    double total_vdw = 0.0;
    
    // Sum energy components from all active residues
    for (int i = 0; i < state.activeResidueCount; ++i) {
        if (state.residues[i].active) {
            total_elec += state.residues[i].energy_elec;
            total_vdw += state.residues[i].energy_vdw;
        }
    }
    
    // For direct energy calculations, divide by 2 to account for double counting
    // Note: For Ewald/PME methods, the reciprocal and self energy terms should be added separately
    total_elec /= 2.0;
    total_vdw /= 2.0;
    
    return std::make_pair(total_elec, total_vdw);
}

void Simulation::computeSystemEnergyPGPComplete(model::MCState& state) {
    if (is_debug_enabled()) {
        log(LogLevel::DEBUG, "Computing PGP Complete energy for all active residues");
    }
    platform::cpu::computeSystemEnergyPGPComplete(state);
}

void Simulation::computeMovementEnergyPGPComplete(model::MCState& state) {
    if (is_debug_enabled()) {
        log(LogLevel::DEBUG, "Computing PGP Complete energy for movement residues");
    }
    platform::cpu::computeMovementEnergyPGPComplete(state);
}

} // namespace simulation
} // namespace pygcmc