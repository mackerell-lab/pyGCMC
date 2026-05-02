#ifndef EWALDSYSTEMENERGY_HPP
#define EWALDSYSTEMENERGY_HPP

#include "model/ModelModule.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Compute total system energy using Ewald summation
 *
 * @param state MC state
 */
void computeSystemEwaldEnergy(model::MCState& state);

/**
 * @brief Validate Ewald setup and parameters
 *
 * @param state MC state
 * @return true if setup is valid
 */
bool validateEwaldSetup(const model::MCState& state);

/**
 * @brief Get Ewald energy components breakdown
 *
 * @param state MC state
 * @param realSpace Real space energy
 * @param reciprocal Reciprocal space energy
 * @param selfEnergy Self energy
 * @param vdw Van der Waals energy
 * @param total Total energy
 */
void getEwaldEnergyBreakdown(const model::MCState& state,
                            double& realSpace,
                            double& reciprocal,
                            double& selfEnergy,
                            double& vdw,
                            double& total);

/**
 * @brief Validate system properties
 *
 * @param state MC state
 */
void validateEwaldSystemProperties(const model::MCState& state);

} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // EWALDSYSTEMENERGY_HPP
