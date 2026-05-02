#ifndef PMEGRIDMAPPING_HPP
#define PMEGRIDMAPPING_HPP

#include "model/ModelModule.hpp"
#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Calculate grid indices and fractional coordinates for atoms
 *
 * @param state MC state
 * @param atomsToProcess List of atom indices to process
 * @param recipBoxVectors Reciprocal lattice vectors
 * @param gridIndices Output grid indices for each atom
 * @param gridFractions Output fractional coordinates for each atom
 */
void calculateGridIndicesAndFractions(const model::MCState& state,
                                     const std::vector<int>& atomsToProcess,
                                     const double recipBoxVectors[3][3],
                                     std::vector<std::vector<int>>& gridIndices,
                                     std::vector<std::vector<double>>& gridFractions);

/**
 * @brief Calculate B-spline coefficients for all atoms
 *
 * @param atomsToProcess List of atom indices to process
 * @param gridFractions Fractional coordinates for each atom
 * @param bsplines_theta Output B-spline coefficients
 */
void calculateBSplineCoefficients(const std::vector<int>& atomsToProcess,
                                 const std::vector<std::vector<double>>& gridFractions,
                                 std::vector<std::vector<double>>& bsplines_theta);

/**
 * @brief Distribute charges to the grid
 *
 * @param state MC state
 * @param atomsToProcess List of atom indices to process
 * @param gridIndices Grid indices for each atom
 * @param bsplines_theta B-spline coefficients
 * @return Total charge distributed to grid
 */
double distributeChargesToGrid(const model::MCState& state,
                              const std::vector<int>& atomsToProcess,
                              const std::vector<std::vector<int>>& gridIndices,
                              const std::vector<std::vector<double>>& bsplines_theta);

} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PMEGRIDMAPPING_HPP
