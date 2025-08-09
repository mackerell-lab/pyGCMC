#include "DrudeNBTholeBuilder.hpp"
#include "DrudeSCFOM.hpp"
#include <cmath>
#include <array>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace exp {

void buildNBTholePairs(const std::vector<DrudeParticle>& particles,
                       double thole,
                       double /* cutoff_nm */,
                       std::vector<ScreenedPair>& outPairs) {
    // Simple version: add all pairs without distance checking
    // This is suitable when we don't have access to state/positions
    // or when we want all possible pairs
    
    outPairs.reserve(outPairs.size() + particles.size() * (particles.size() - 1) / 2);
    
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            ScreenedPair pair;
            pair.dipole1 = static_cast<int>(i);
            pair.dipole2 = static_cast<int>(j);
            pair.thole = thole;
            outPairs.push_back(pair);
        }
    }
}

void buildNBTholePairsWithCutoff(const std::vector<DrudeParticle>& particles,
                                  const model::MCState& state,
                                  double thole,
                                  double cutoff_nm,
                                  std::vector<ScreenedPair>& outPairs) {
    // Advanced version: only add pairs within cutoff distance
    // based on parent-parent distances
    
    outPairs.reserve(outPairs.size() + particles.size() * 10);  // Estimate
    
    std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
    double cutoff2 = cutoff_nm * cutoff_nm;
    
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle1 = particles[i];
        if (particle1.parentIndex < 0 || particle1.parentIndex >= state.activeAtomCount) {
            continue;
        }
        const auto& parent1 = state.atoms[particle1.parentIndex];
        
        for (size_t j = i + 1; j < particles.size(); ++j) {
            const auto& particle2 = particles[j];
            if (particle2.parentIndex < 0 || particle2.parentIndex >= state.activeAtomCount) {
                continue;
            }
            const auto& parent2 = state.atoms[particle2.parentIndex];
            
            // Calculate parent-parent distance
            double dx = parent2.x - parent1.x;
            double dy = parent2.y - parent1.y;
            double dz = parent2.z - parent1.z;
            
            // Apply PBC
            DrudeSCFOM::applyPBC(dx, dy, dz, box);
            
            double r2 = dx*dx + dy*dy + dz*dz;
            
            // Check cutoff
            if (r2 < cutoff2) {
                ScreenedPair pair;
                pair.dipole1 = static_cast<int>(i);
                pair.dipole2 = static_cast<int>(j);
                pair.thole = thole;
                outPairs.push_back(pair);
            }
        }
    }
}

} // namespace exp
} // namespace cpu
} // namespace platform
} // namespace pygcmc