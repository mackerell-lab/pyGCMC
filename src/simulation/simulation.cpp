// src/simulation/simulation.cpp
#include "simulation.hpp"
#include "../platform/cpu/energy.hpp"
#include "../platform/cpu/energyPME.hpp"
#include "../platform/cpu/energyPGP.hpp"
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

// Implementation of PGP-related methods
void Simulation::setPGPParameters(float alpha, const int meshSize[3], float pair_cutoff, 
                                  const int pairGridSize[3], int splineOrder, float tolerance) {
    // Convert float parameters to double
    double alpha_d = static_cast<double>(alpha);
    double pair_cutoff_d = static_cast<double>(pair_cutoff);
    double tolerance_d = static_cast<double>(tolerance);
    
    // Call the platform implementation
    platform::cpu::setPGPParameters(alpha_d, meshSize, pair_cutoff_d, pairGridSize, splineOrder, tolerance_d);
    
    log(LogLevel::INFO, "PGP parameters set: alpha=", alpha, 
        ", meshSize=[", meshSize[0], ",", meshSize[1], ",", meshSize[2], "]",
        ", pair_cutoff=", pair_cutoff,
        ", pairGridSize=[", pairGridSize[0], ",", pairGridSize[1], ",", pairGridSize[2], "]",
        ", splineOrder=", splineOrder,
        ", tolerance=", tolerance);
}

void Simulation::precomputeGridPotential(model::MCState& state, bool fixed_only) {
    platform::cpu::precomputeGridPotential(state, fixed_only);
}

void Simulation::interpolateMoleculeEnergy(model::MCState& state, double& energy) {
    platform::cpu::interpolateMoleculeEnergy(state, energy);
}

double Simulation::calculateMoleculeEnergy(model::MCState& state) {
    return platform::cpu::calculateMoleculeEnergy(state);
}

} // namespace simulation
} // namespace pygcmc