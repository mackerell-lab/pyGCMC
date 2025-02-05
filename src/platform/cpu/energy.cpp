// src/platform/cpu/energy.cpp
#include "energy.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

float computeNaiveNonbondedEnergy(const model::MCState& state) {
    float energy = 0.0f;
    const auto& residues = state.residues;
    const auto& forcefield = state.forcefield;
    const float r = 1.0f;  // Fixed distance for naive implementation

    // Iterate through all active movement molecules
    for (const auto& movementInfo : state.movementResidues) {
        int start = movementInfo.startIndex;
        int count = movementInfo.activeCount;
        
        // For each active movement molecule
        for (int i = start; i < start + count; ++i) {
            if (!residues[i].active) continue;
            int moveType = residues[i].type;
            
            // Compute interaction with all other active molecules
            for (size_t j = 0; j < residues.size(); ++j) {
                if (!residues[j].active || static_cast<int>(j) == i) continue;
                
                int resType = residues[j].type;
                // For movement type 0 and fixed type 1, index should be 0 * maxTypes + 1 = 1
                // So we need to use the correct index in the force field arrays
                size_t index = static_cast<size_t>(moveType * forcefield.maxTypes + resType);
                if (index >= forcefield.ljEps.size()) continue;  // Safety check
                
                float eps = forcefield.ljEps[index];
                float sigma = forcefield.ljSigma[index];
                
                // Debug output
                printf("Computing energy between residues %d (type %d) and %zu (type %d)\n", i, moveType, j, resType);
                printf("Using force field parameters: eps = %f, sigma = %f at index %zu\n", eps, sigma, index);
                
                float term6 = std::pow(sigma / r, 6);
                float term12 = term6 * term6;
                
                // V = eps * [(sigma/r)^12 - 2*(sigma/r)^6]
                float pairEnergy = eps * (term12 - 2.0f * term6);
                printf("Computed pair energy: %f\n", pairEnergy);
                
                energy += pairEnergy;
            }
        }
    }
    
    printf("Total energy: %f\n", energy);
    return energy;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc