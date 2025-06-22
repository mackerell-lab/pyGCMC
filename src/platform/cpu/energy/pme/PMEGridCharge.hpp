#ifndef PMEGRIDCHARGE_HPP
#define PMEGRIDCHARGE_HPP

#include "model/ModelModule.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Spread charges onto the PME grid
 * 
 * @param state MC state
 * @param fixed_only If true, only consider charges from atoms in fixed residues
 */
void spreadChargesOntoGrid(model::MCState& state, bool fixed_only);

/**
 * @brief Calculate fractional coordinates in the grid
 */
void calculateFractionalCoordinates(const float position[3], 
                                  const float box[3], 
                                  const int meshSize[3], 
                                  double fractional[3]);

} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PMEGRIDCHARGE_HPP 