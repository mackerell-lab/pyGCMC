#include "PGPCore.hpp"
#include "../../../platform.hpp"
#include <thread>
#include <chrono>

namespace pygcmc {
namespace platform {
namespace cpu {

// Define the global mutex
std::mutex pgp_mutex;

void resetPGPState() {
    // Lock the mutex to ensure thread safety
    std::lock_guard<std::mutex> lock(pgp_mutex);
    
    platform::log(LogLevel::INFO, "Resetting PGP global state");
    
    // Clear the potential grid using swap trick to ensure memory is freed
    std::vector<std::complex<double>>().swap(pgp_params.potentialGrid);
    
    // Reset grid parameters
    pgp_params.potential_cutoff = 0.0;
    pgp_params.grid_spacing = 0.0;
    for (int i = 0; i < 3; i++) {
        pgp_params.potential_grid_size[i] = 0;
    }
    
    // Clear inherited PME data structures
    pgp_params.pmeGrid.clear();
    pgp_params.pmeGrid.shrink_to_fit();
    
    pgp_params.pmeCharge.clear();
    pgp_params.pmeCharge.shrink_to_fit();
    
    // Clear B-spline moduli
    for (int i = 0; i < 3; i++) {
        pgp_params.bsplineModuli[i].clear();
        pgp_params.bsplineModuli[i].shrink_to_fit();
    }
    
    // Clear lookup tables
    pgp_params.erfcTable.clear();
    pgp_params.erfcTable.shrink_to_fit();
    
    pgp_params.ewaldScaleTable.clear(); 
    pgp_params.ewaldScaleTable.shrink_to_fit();
    
    // Reset scalar parameters
    pgp_params.initialized = false;
    pgp_params.alpha = 0.0;
    pgp_params.tolerance = 1e-5f;
    pgp_params.epsilon_r = 1.0;
    pgp_params.splineOrder = 4;
    pgp_params.ewaldDX = 0.0;
    pgp_params.ewaldDXInv = 0.0;
    pgp_params.erfcDXInv = 0.0;
    pgp_params.cutoff = 0.0;
    pgp_params.debug_mode = true;
    
    // Reset arrays
    for (int i = 0; i < 3; i++) {
        pgp_params.box[i] = 0.0;
        pgp_params.meshSize[i] = 0;
    }
    
    platform::log(LogLevel::INFO, "PGP state reset complete");
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc