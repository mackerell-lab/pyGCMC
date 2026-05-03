#include "PGPCore.hpp"
#include "PGPGlobal.hpp"
#include "platform/cpu/energy/pme/PMEComposite.hpp"
#include "platform/cpu/energy/pme/PMEGlobal.hpp"
#include "platform/cpu/energy/pme/PMESetup.hpp"
#include "platform/cpu/energy/common/MemorySafetyChecks.hpp"
#include "platform/Platform.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

// Note: getPGPParams() is now defined via macro in PGPGlobal.hpp

// Mathematical constants
const double TWO_PI = 2.0 * M_PI;
const double SQRT_PI = sqrt(M_PI);

void setPGPParameters(double alpha, const int meshSize[3], double potential_cutoff,
                        const int potentialGridSize[3], int splineOrder, double tolerance) {
    platform::log(LogLevel::DEBUG, "setPGPParameters called");

    // Mark as initialized
    SafetyChecks::markInitialized();

    // Lock for thread safety
    platform::log(LogLevel::DEBUG, "Acquiring lock...");
    std::lock_guard<std::mutex> lock(pgp_global_mutex);
    platform::log(LogLevel::DEBUG, "Lock acquired");

    // Ensure pgp_params_ptr exists before using it
    if (!pgp_params_ptr) {
        pgp_params_ptr = std::make_unique<PGPParams>();
    }

    // Now we can safely use the pointer directly (we already have the lock)
    // Note: We use getPGPParams() throughout for consistency

    // Clear existing potential grid before setting new parameters
    if (!getPGPParams().potentialGrid.empty()) {
        platform::log(LogLevel::DEBUG, "Clearing existing potential grid");
        std::vector<std::complex<double>>().swap(getPGPParams().potentialGrid);
    }

    // First set standard PME parameters (this will set all inherited fields)
    setPMEParameters(alpha, meshSize, splineOrder, tolerance);

    // Copy PME parameters to PGP instance (since they are separate global instances)
    getPGPParams().alpha = getPMEParams().alpha;
    getPGPParams().tolerance = getPMEParams().tolerance;
    getPGPParams().epsilon_r = getPMEParams().epsilon_r;
    getPGPParams().splineOrder = getPMEParams().splineOrder;

    // Copy arrays
    for (int i = 0; i < 3; i++) {
        getPGPParams().box[i] = getPMEParams().box[i];
        getPGPParams().meshSize[i] = getPMEParams().meshSize[i];
    }

    // Copy lookup tables and grid data
    getPGPParams().erfcTable = getPMEParams().erfcTable;
    getPGPParams().ewaldScaleTable = getPMEParams().ewaldScaleTable;
    getPGPParams().ewaldDX = getPMEParams().ewaldDX;
    getPGPParams().ewaldDXInv = getPMEParams().ewaldDXInv;
    getPGPParams().erfcDXInv = getPMEParams().erfcDXInv;

    // Copy B-spline moduli
    for (int i = 0; i < 3; i++) {
        getPGPParams().bsplineModuli[i] = getPMEParams().bsplineModuli[i];
    }

    // Copy PME grids
    getPGPParams().pmeGrid = getPMEParams().pmeGrid;
    getPGPParams().pmeCharge = getPMEParams().pmeCharge;

    // Set PGP-specific parameters
    getPGPParams().potential_cutoff = potential_cutoff;
    for (int i = 0; i < 3; i++) {
        getPGPParams().potential_grid_size[i] = potentialGridSize[i];
    }

    // Initialize grid for precomputed potential
    getPGPParams().initializePotentialGrid();

    // Mark PGP as initialized
    getPGPParams().initialized = true;

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
    // Set PGP parameters (this will also set PME parameters via setPMEParameters)
    setPGPParameters(alpha, meshSize, potentialCutoff, potentialGridSize, splineOrder, tolerance);

    // Initialize PME parameters (PGP inherits from PME)
    // Use getPGPParams() here since we're not holding the lock
    getPGPParams().setBox(box);
    getPGPParams().cutoff = cutoff;
    getPGPParams().initializeTables(cutoff);
    getPGPParams().initializeBsplines();

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
