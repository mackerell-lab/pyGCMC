#pragma once

#include "../model/montecarlo.hpp"
#include <string>

namespace pygcmc {
namespace platform {

class IPlatform {
public:
    virtual ~IPlatform() = default;

    // Initializes the platform with the initial state (e.g. uploading data to GPU or storing for CPU computation)
    virtual void initialize(const model::MCState& state) = 0;

    // Finalizes the platform by cleaning up resources or optionally downloading final data to CPU
    virtual void finalize() = 0;

    // Computes the total non-bonded energy of the current system state
    virtual float computeTotalEnergy() = 0;

    // Computes the energy contribution for a specific residue
    virtual float computeResidueEnergy(int residueIndex) = 0;

    // Attempts a Monte Carlo move (e.g., translate, insert, delete) and updates the state accordingly
    virtual bool attemptMove(const std::string& moveType) = 0;
};

} // namespace platform
} // namespace pygcmc 