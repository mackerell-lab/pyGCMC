#include "PGPGrid.hpp"
#include <algorithm>
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

// Note: pgp_params is defined in PGPCore.cpp

/**
 * @brief Initialize the three-dimensional grid for precomputed potential
 * 
 * This is a fundamental step in the PGP-PME algorithm, responsible for creating and initializing the 3D grid used to store precomputed potentials.
 * This function allocates grid memory based on potential_grid_size parameters and calculates appropriate grid spacing.
 */
void PGPParams::initializePotentialGrid() {
    // Calculate total grid size and allocate memory
    int totalSize = potential_grid_size[0] * potential_grid_size[1] * potential_grid_size[2];
    potentialGrid.resize(totalSize);
    
    // Set grid spacing, taking the minimum value of the three dimensions
    grid_spacing = std::min({
        box[0] / potential_grid_size[0],
        box[1] / potential_grid_size[1],
        box[2] / potential_grid_size[2]
    });
    
    // Only output debug information in debug mode
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "PGP potential grid initialized with size: ", 
                    potential_grid_size[0], "x", potential_grid_size[1], "x", potential_grid_size[2],
                    ", grid spacing: ", grid_spacing);
    }
}

void initializePotentialGridImpl() {
    // Implementation wrapper for external access
    pgp_params.initializePotentialGrid();
}

// <agent-hook:pgp_grid_impl>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 