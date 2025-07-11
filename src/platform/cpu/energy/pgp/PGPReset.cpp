/**
 * Reset function for PGP state - PROPOSED FIX
 * 
 * This function should be added to properly clean up global state
 * between test runs or when reinitializing PGP parameters.
 */

#include "PGPCore.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

void resetPGPState() {
    // Clear the potential grid
    pgp_params.potentialGrid.clear();
    pgp_params.potentialGrid.shrink_to_fit();  // Actually release memory
    
    // Reset grid parameters
    pgp_params.potential_cutoff = 0.0;
    for (int i = 0; i < 3; i++) {
        pgp_params.potential_grid_size[i] = 0;
    }
    pgp_params.grid_spacing = 0.0;
    
    // Clear inherited PME data
    pgp_params.pmeGrid.clear();
    pgp_params.pmeCharge.clear();
    
    // Reset all other parameters to defaults
    pgp_params.initialized = false;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc