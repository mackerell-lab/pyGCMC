#include "PGPCore.hpp"
#include "platform/cpu/energy/pme/PMEComposite.hpp"
#include "platform/cpu/energy/pme/PMESetup.hpp"
#include "platform/platform.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

// Global PGP parameters instance
PGPParams pgp_params;

// Mathematical constants
const double TWO_PI = 2.0 * M_PI;
const double SQRT_PI = sqrt(M_PI);

void setPGPParameters(double alpha, const int meshSize[3], double potential_cutoff, 
                        const int potentialGridSize[3], int splineOrder, double tolerance) {
    // First set standard PME parameters
    setPMEParameters(alpha, meshSize, splineOrder, tolerance);
    
    // Copy standard PME parameters to PGP parameter structure
    pgp_params.alpha = pme_params.alpha;
    pgp_params.tolerance = pme_params.tolerance;
    pgp_params.initialized = pme_params.initialized;
    pgp_params.cutoff = pme_params.cutoff;
    pgp_params.epsilon_r = pme_params.epsilon_r;
    pgp_params.splineOrder = pme_params.splineOrder;
    
    // Copy box size and grid size
    for (int i = 0; i < 3; i++) {
        pgp_params.box[i] = pme_params.box[i];
        pgp_params.meshSize[i] = pme_params.meshSize[i];
    }
    
    // Copy PME lookup tables
    pgp_params.erfcTable = pme_params.erfcTable;
    pgp_params.ewaldScaleTable = pme_params.ewaldScaleTable;
    pgp_params.ewaldDX = pme_params.ewaldDX;
    pgp_params.ewaldDXInv = pme_params.ewaldDXInv;
    pgp_params.erfcDXInv = pme_params.erfcDXInv;
    
    // Copy B-spline moduli
    for (int i = 0; i < 3; i++) {
        pgp_params.bsplineModuli[i] = pme_params.bsplineModuli[i];
    }
    
    // Copy PME grid
    pgp_params.pmeGrid = pme_params.pmeGrid;
    pgp_params.pmeCharge = pme_params.pmeCharge;
    
    // Set PGP-specific parameters
    pgp_params.potential_cutoff = potential_cutoff;
    for (int i = 0; i < 3; i++) {
        pgp_params.potential_grid_size[i] = potentialGridSize[i];
    }
    
    // Initialize grid for precomputed potential
    pgp_params.initializePotentialGrid();
    
    // Mark as initialized
    pgp_params.initialized = true;
    
    // Output parameter setting information
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "PGP parameters set: alpha=", alpha, 
                    ", potential_cutoff=", potential_cutoff, 
                    ", potentialGrid=[", potentialGridSize[0], ",", potentialGridSize[1], ",", potentialGridSize[2], "]");
    }
}

void initializePGPParameters(double cutoff, const double box[3], 
                           double alpha,
                           const int meshSize[3], 
                           double potentialCutoff,
                           const int potentialGridSize[3],
                           int splineOrder,
                           double tolerance) {
    // Set box dimensions
    pgp_params.box[0] = box[0];
    pgp_params.box[1] = box[1];
    pgp_params.box[2] = box[2];
    pgp_params.cutoff = cutoff;
    
    // Set PGP parameters
    setPGPParameters(alpha, meshSize, potentialCutoff, potentialGridSize, splineOrder, tolerance);
    
    // Initialize PME parameters (PGP inherits from PME)
    pgp_params.setBox(box);
    pgp_params.initializeTables(cutoff);
    pgp_params.initializeBsplines();
    
    platform::log(LogLevel::INFO, "PGP parameters initialized successfully");
}

// Note: The following functions are implemented in PGPEvaluator.cpp:
// - computeSystemEnergyPGP
// - computeMovementEnergyPGP  
// - computeRealSpacePGP
// - computeSelfEnergyPGP

// <agent-hook:pgp_core_impl>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 