#ifndef PMEGRIDUTILS_HPP
#define PMEGRIDUTILS_HPP

#include <vector>
#include <complex>

namespace pygcmc {
namespace platform {
namespace cpu {

// Forward declaration
struct PMEParams;

/**
 * @brief Initialize and clear the PME grid
 */
void initializePMEGrid(PMEParams& params);

/**
 * @brief Clear all values in the PME grid
 */
void clearGrid(std::vector<std::complex<double>>& grid);

/**
 * @brief Validate grid indices for bounds checking
 */
bool validateGridIndices(int x, int y, int z, int nx, int ny, int nz);

/**
 * @brief Get grid statistics for debugging
 */
void getGridStatistics(const std::vector<std::complex<double>>& grid,
                      int& nonZeroCount, 
                      double& maxValue, 
                      double& totalCharge);

} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PMEGRIDUTILS_HPP 