// src/simulation/simulation.hpp

#pragma once

#include "../platform/platform.hpp"
#include "../model/montecarlo.hpp"
#include <memory>

namespace pygcmc {
namespace simulation {

class Simulation {
public:
    // Construct the simulation with a platform (e.g. CPU or CUDA)
    explicit Simulation(std::unique_ptr<platform::IPlatform> platform)
      : platform_(std::move(platform)) {}

    // Upload the initial state to the platform
    void initialize(const model::MCState& initialState) {
        platform_->initialize(initialState);
    }

    // Run the simulation for a number of steps by attempting moves
    void run(int steps) {
        for (int i = 0; i < steps; ++i) {
            // For simplicity, using "translate" as the move type
            platform_->attemptMove("translate");
        }
    }

    // Compute and return the total energy from the platform
    float computeEnergy() {
        return platform_->computeTotalEnergy();
    }

    // Compute naive nonbonded energy between active movement residues and all other active residues
    // Note: This function modifies the energy_vdw and energy_elec parameters of residues in the state
    static float computeNaiveNonbondedEnergy(model::MCState& state);

    // Finalize the simulation and (optionally) download final data
    void finalize() {
        platform_->finalize();
    }

private:
    std::unique_ptr<platform::IPlatform> platform_;
};

} // namespace simulation
} // namespace pygcmc