// src/simulation/simulation.cpp
#include "simulation.hpp"
#include "../platform/cpu/energy.hpp"
#include <cmath>

namespace pygcmc {
namespace simulation {

float Simulation::computeNaiveNonbondedEnergy(const model::MCState& state) {
    return platform::cpu::computeNaiveNonbondedEnergy(state);
}

} // namespace simulation
} // namespace pygcmc