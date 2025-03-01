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
    // 转换 float 参数为 double
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

} // namespace simulation
} // namespace pygcmc