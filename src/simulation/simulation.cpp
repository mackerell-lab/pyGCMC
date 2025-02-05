// src/simulation/simulation.cpp
#include "simulation.hpp"
#include "../platform/cpu/energy.hpp"
#include <cmath>

namespace pygcmc {
namespace simulation {

float Simulation::computeNaiveNonbondedEnergy(model::MCState& state) {
    log(LogLevel::DEBUG, "Computing naive nonbonded energy for first movement residue");
    
    platform::cpu::computeNaiveNonbondedEnergy(state);
    
    // Return the energy of the first movement residue
    const auto& firstMovementInfo = state.movementResidues[0];
    const auto& residue = state.residues[firstMovementInfo.startIndex];
    float energy = residue.energy_vdw + residue.energy_elec;
    
    log(LogLevel::DEBUG, "First movement residue energy: vdw=", residue.energy_vdw, 
        ", elec=", residue.energy_elec, ", total=", energy);
    
    return energy;
}

} // namespace simulation
} // namespace pygcmc