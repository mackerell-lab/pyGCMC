// src/simulation/simulation.cpp
#include "simulation.hpp"
#include "../platform/cpu/energy.hpp"
#include <cmath>

namespace pygcmc {
namespace simulation {

void Simulation::computeNaiveNonbondedEnergy(model::MCState& state) {
    log(LogLevel::DEBUG, "Computing naive nonbonded energy for movement residues");
    
    platform::cpu::computeNaiveNonbondedEnergy(state);
    
    // Log the energy of the first movement residue
    const auto& firstMovementInfo = state.movementResidues[0];
    const auto& residue = state.residues[firstMovementInfo.startIndex];
    
    log(LogLevel::DEBUG, "First movement residue energy: vdw=", residue.energy_vdw, 
        ", elec=", residue.energy_elec, 
        ", total=", (residue.energy_vdw + residue.energy_elec));
}

} // namespace simulation
} // namespace pygcmc