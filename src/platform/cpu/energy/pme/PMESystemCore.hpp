#ifndef PMEENERGYCALC_HPP
#define PMEENERGYCALC_HPP

#include "model/ModelModule.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Compute energy from the PME grid after FFT
 *
 * @param energy Output energy
 * @param box Box dimensions
 */
void computeEnergyFromGrid(double& energy, const double box[3]);

/**
 * @brief Compute reciprocal space energy using PME
 *
 * @param state MC state
 * @return Reciprocal space energy
 */
double computeReciprocalPME(model::MCState& state);

} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PMEENERGYCALC_HPP
