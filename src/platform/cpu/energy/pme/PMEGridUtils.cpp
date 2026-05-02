#include "PMEGridUtils.hpp"
#include "PMECore.hpp"
#include <algorithm>
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Initialize and clear the PME grid
 */
void initializePMEGrid(PMEParams& params) {
    int totalGridPoints = params.meshSize[0] * params.meshSize[1] * params.meshSize[2];
    params.pmeGrid.resize(totalGridPoints);
    clearGrid(params.pmeGrid);
}

/**
 * @brief Clear all values in the PME grid
 */
void clearGrid(std::vector<std::complex<double>>& grid) {
    std::fill(grid.begin(), grid.end(), std::complex<double>(0.0, 0.0));
}

/**
 * @brief Validate grid indices for bounds checking
 */
bool validateGridIndices(int x, int y, int z, int nx, int ny, int nz) {
    return (x >= 0 && x < nx && y >= 0 && y < ny && z >= 0 && z < nz);
}

/**
 * @brief Get grid statistics for debugging
 */
void getGridStatistics(const std::vector<std::complex<double>>& grid,
                      int& nonZeroCount,
                      double& maxValue,
                      double& totalCharge) {
    nonZeroCount = 0;
    maxValue = 0.0;
    totalCharge = 0.0;

    for (const auto& point : grid) {
        double magnitude = std::abs(point);
        if (magnitude > 1e-10) {
            nonZeroCount++;
            maxValue = std::max(maxValue, magnitude);
        }
        totalCharge += point.real();
    }
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
