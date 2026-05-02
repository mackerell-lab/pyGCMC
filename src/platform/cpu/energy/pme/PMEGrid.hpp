#pragma once

#include "model/ModelModule.hpp"
#include "PMECore.hpp"
#include <complex>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Grid operations for PME calculations
 *
 * This module handles the charge spreading and grid manipulation operations
 * required for PME. It includes functions for distributing particle charges
 * onto the grid and managing grid data structures.
 */

/**
 * @brief Spread particle charges onto the PME grid
 *
 * This function distributes the charges of particles onto the PME grid using
 * B-spline interpolation. This is a critical step in the PME algorithm.
 *
 * @param state MC state containing particle positions and charges
 * @param fixed_only If true, only process charges from fixed residues
 */
void spreadChargesOntoGrid(model::MCState& state, bool fixed_only = false);

/**
 * @brief Initialize and clear the PME grid
 *
 * @param params PME parameters containing grid dimensions
 */
void initializePMEGrid(PMEParams& params);

/**
 * @brief Clear all values in the PME grid
 *
 * @param grid Grid to clear
 */
void clearGrid(std::vector<std::complex<double>>& grid);

/**
 * @brief Validate grid indices for bounds checking
 *
 * @param x, y, z Grid indices
 * @param nx, ny, nz Grid dimensions
 * @return true if indices are valid
 */
bool validateGridIndices(int x, int y, int z, int nx, int ny, int nz);

/**
 * @brief Convert 3D grid coordinates to linear index
 *
 * @param x, y, z Grid coordinates
 * @param ny, nz Grid dimensions in Y and Z
 * @return Linear grid index
 */
inline int gridIndex3D(int x, int y, int z, int ny, int nz) {
    return x * ny * nz + y * nz + z;
}

/**
 * @brief Calculate fractional coordinates in the grid
 *
 * @param position Real space position
 * @param box Box dimensions
 * @param meshSize Grid dimensions
 * @param fractional Output fractional coordinates
 */
void calculateFractionalCoordinates(const float position[3],
                                  const float box[3],
                                  const int meshSize[3],
                                  double fractional[3]);

/**
 * @brief Get grid statistics for debugging
 *
 * @param grid PME grid
 * @param nonZeroCount Output: number of non-zero grid points
 * @param maxValue Output: maximum grid value magnitude
 * @param totalCharge Output: total charge on grid
 */
void getGridStatistics(const std::vector<std::complex<double>>& grid,
                      int& nonZeroCount,
                      double& maxValue,
                      double& totalCharge);

// <agent-hook:pme_grid>

} // namespace cpu
} // namespace platform
} // namespace pygcmc
