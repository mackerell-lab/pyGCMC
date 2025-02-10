// src/simulation/simulation.cpp
#include "simulation.hpp"
#include "../platform/cpu/energy.hpp"
#include <cmath>

namespace pygcmc {
namespace simulation {

void Simulation::computeMovementResiduesEnergy(model::MCState& state) {
    log(LogLevel::DEBUG, "Computing nonbonded energy for movement residues");
    
    platform::cpu::computeMovementResiduesEnergy(state);
    
    // Log the energy of the first movement residue
    const auto& firstMovementInfo = state.movementResidues[0];
    const auto& residue = state.residues[firstMovementInfo.startIndex];
    
    log(LogLevel::DEBUG, "First movement residue energy: vdw=", residue.energy_vdw, 
        ", elec=", residue.energy_elec, 
        ", total=", (residue.energy_vdw + residue.energy_elec));
}

void Simulation::computeFullSystemEnergy(model::MCState& state) {
    log(LogLevel::DEBUG, "Computing nonbonded energy for all active residues");
    
    platform::cpu::computeFullSystemEnergy(state);
    
    // Log total system energy
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

void Simulation::computeFullSystemCutoffEnergy(model::MCState& state) {
    log(LogLevel::DEBUG, "Computing cutoff nonbonded energy for all active residues");
    
    platform::cpu::computeFullSystemCutoffEnergy(state);
    
    // Log total system energy
    float total_vdw = 0.0f;
    float total_elec = 0.0f;
    for (int i = 0; i < state.activeResidueCount; ++i) {
        if (state.residues[i].active) {
            total_vdw += state.residues[i].energy_vdw;
            total_elec += state.residues[i].energy_elec;
        }
    }
    
    // Total energy divided by 2 (since each interaction is counted twice)
    total_vdw /= 2.0f;
    total_elec /= 2.0f;
    
    log(LogLevel::DEBUG, "Total system energy: vdw=", total_vdw, 
        ", elec=", total_elec, 
        ", total=", (total_vdw + total_elec));
}

} // namespace simulation
} // namespace pygcmc