#pragma once
#include "../DrudeStructures.hpp"
#include "model/ModelModule.hpp"
#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace exp {

// Build NBTHOLE pairs for all Drude particle pairs
void buildNBTholePairs(const std::vector<DrudeParticle>& particles,
                       double thole,
                       double cutoff_nm,
                       std::vector<ScreenedPair>& outPairs);

// Build NBTHOLE pairs based on distance cutoff (requires state)
void buildNBTholePairsWithCutoff(const std::vector<DrudeParticle>& particles,
                                  const model::MCState& state,
                                  double thole,
                                  double cutoff_nm,
                                  std::vector<ScreenedPair>& outPairs);

} // namespace exp
} // namespace cpu
} // namespace platform
} // namespace pygcmc