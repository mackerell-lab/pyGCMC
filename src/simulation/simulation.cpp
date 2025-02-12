// src/simulation/simulation.cpp
#include "simulation.hpp"
#include "../platform/cpu/energy.hpp"
#include <cmath>

namespace pygcmc {
namespace simulation {

void Simulation::computeMovementEnergy(model::MCState& state) {
    if (is_debug_enabled()) {
        log(LogLevel::DEBUG, "Computing nonbonded energy for movement residues");
    }
    platform::cpu::computeMovementEnergy(state);
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
    
    platform::cpu::computeSystemEnergy(state);
    
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
        log(LogLevel::DEBUG, "Computing PBC nonbonded energy: box=", 
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
        
        log(LogLevel::DEBUG, "Total system energy with PBC: vdw=", total_vdw, 
            ", elec=", total_elec, 
            ", total=", (total_vdw + total_elec));
    }
}

void Simulation::setEnergyDebugOutput(bool enable) {
    // Set both simulation-level and platform-level debug output
    if (enable) {
        set_verbose(true);
        set_log_level(LogLevel::DEBUG);
    }
    platform::cpu::setEnergyDebugOutput(enable);
}

void Simulation::setEwaldParameters(float alpha, const int kmax[3], float tolerance) {
    platform::cpu::setEwaldParameters(alpha, kmax, tolerance);
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

} // namespace simulation
} // namespace pygcmc