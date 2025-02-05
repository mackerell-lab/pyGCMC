// src/simulation/simulation.cpp
#include "simulation.hpp"
#include "../platform/cpu/energy.hpp"
#include <cmath>

namespace pygcmc {
namespace simulation {

float Simulation::computeNaiveNonbondedEnergy(model::MCState& state) {
    platform::cpu::computeNaiveNonbondedEnergy(state);
    // Return the energy of the first movement residue
    const auto& firstMovementInfo = state.movementResidues[0];
    const auto& residue = state.residues[firstMovementInfo.startIndex];
    return residue.energy_vdw + residue.energy_elec;
}

} // namespace simulation
} // namespace pygcmc