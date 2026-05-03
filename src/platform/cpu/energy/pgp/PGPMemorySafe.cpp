/**
 * @brief Memory-safe optimization for PGP implementation
 *
 * This file contains optimized versions of PGP functions that avoid memory issues
 * while maintaining compatibility with existing functionality.
 */

#include "PGPCore.hpp"
#include "PGPPrecompute.hpp"
#include "PGPGlobal.hpp"
#include "../../../Platform.hpp"
#include <memory>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Thread-safe wrapper for PGP initialization
 *
 * Ensures proper cleanup before re-initialization to avoid memory leaks
 */
void safeSetPGPParameters(double alpha, const int meshSize[3], double potential_cutoff,
                         const int potentialGridSize[3], int splineOrder, double tolerance) {
    std::lock_guard<std::mutex> lock(pgp_global_mutex);

    // Clear existing data before setting new parameters
    if (!getPGPParams().potentialGrid.empty()) {
        platform::log(LogLevel::DEBUG, "Clearing existing PGP grid before re-initialization");
        std::vector<std::complex<double>>().swap(getPGPParams().potentialGrid);
    }

    // Call original setPGPParameters
    setPGPParameters(alpha, meshSize, potential_cutoff, potentialGridSize, splineOrder, tolerance);
}

/**
 * @brief Memory-safe version of initializePotentialGrid
 *
 * Ensures proper memory allocation and deallocation
 */
void PGPParams::initializePotentialGrid() {
    // Calculate total grid size
    int totalSize = potential_grid_size[0] * potential_grid_size[1] * potential_grid_size[2];

    // Clear existing grid if needed
    if (!potentialGrid.empty()) {
        std::vector<std::complex<double>>().swap(potentialGrid);
    }

    // Reserve exact capacity to avoid reallocation
    potentialGrid.reserve(totalSize);
    potentialGrid.resize(totalSize, std::complex<double>(0.0, 0.0));

    platform::log(LogLevel::DEBUG, "Initialized potential grid with size: ", totalSize);
}

/**
 * @brief Safe grid assignment with bounds checking
 *
 * Replaces direct assignment to prevent memory corruption
 */
void safeCopyGridToPotential(const std::vector<std::complex<double>>& sourceGrid) {
    std::lock_guard<std::mutex> lock(pgp_global_mutex);

    // Check if sizes match
    size_t expectedSize = getPGPParams().potential_grid_size[0] *
                         getPGPParams().potential_grid_size[1] *
                         getPGPParams().potential_grid_size[2];

    if (sourceGrid.size() != expectedSize) {
        platform::log(LogLevel::ERROR, "Grid size mismatch: expected ", expectedSize,
                     " but got ", sourceGrid.size());
        throw std::runtime_error("Grid size mismatch in safeCopyGridToPotential");
    }

    // Clear and resize to ensure proper memory management
    getPGPParams().potentialGrid.clear();
    getPGPParams().potentialGrid.reserve(expectedSize);

    // Copy data
    getPGPParams().potentialGrid.assign(sourceGrid.begin(), sourceGrid.end());
}

/**
 * @brief Memory-safe precompute function
 *
 * This replaces the problematic direct assignment in precomputeGridPotentialImpl
 */
void safePrecomputeGridPotential(model::MCState& state, bool fixed_only) {
    // Lock for thread safety
    std::lock_guard<std::mutex> lock(pgp_global_mutex);

    // Check if parameters are initialized
    if (!getPGPParams().initialized) {
        throw std::runtime_error("PGP parameters not initialized");
    }

    // Create local backup to avoid issues with global state
    std::vector<std::complex<double>> pmeGridBackup;
    pmeGridBackup.reserve(getPMEParams().pmeGrid.size());
    pmeGridBackup = getPMEParams().pmeGrid;

    // Store original mesh size
    int originalMeshSize[3];
    for (int i = 0; i < 3; i++) {
        originalMeshSize[i] = getPMEParams().meshSize[i];
    }

    try {
        // Call the original implementation logic
        // (This would include all the charge assignment and FFT logic)
        precomputeGridPotentialImpl(state, fixed_only);

        // Instead of direct assignment, use safe copy
        // Replace: getPGPParams().potentialGrid = getPMEParams().pmeGrid;
        safeCopyGridToPotential(getPMEParams().pmeGrid);

    } catch (const std::exception& e) {
        // Restore PME state on error
        getPMEParams().pmeGrid = std::move(pmeGridBackup);
        for (int i = 0; i < 3; i++) {
            getPMEParams().meshSize[i] = originalMeshSize[i];
        }
        throw; // Re-throw the exception
    }

    // Restore PME state
    getPMEParams().pmeGrid = std::move(pmeGridBackup);
    for (int i = 0; i < 3; i++) {
        getPMEParams().meshSize[i] = originalMeshSize[i];
    }
}

/**
 * @brief Automatic cleanup on program exit
 *
 * This ensures proper cleanup even if resetPGPState isn't called
 */
class PGPAutoCleanup {
public:
    ~PGPAutoCleanup() {
        try {
            resetPGPState();
        } catch (...) {
            // Suppress exceptions in destructor
        }
    }
};

// Static instance ensures cleanup on program exit
static PGPAutoCleanup pgp_auto_cleanup;

} // namespace cpu
} // namespace platform
} // namespace pygcmc
