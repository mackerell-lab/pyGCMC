#include "PGPComposite.hpp"
#include "platform/platform.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

void PGPComposite::initialize(double cutoff, 
                            const double box[3],
                            double alpha,
                            const int meshSize[3],
                            double potentialCutoff,
                            const int potentialGridSize[3],
                            int splineOrder,
                            double tolerance) {
    
    // Log initialization start
    platform::log(LogLevel::INFO, "Initializing PGP parameters...");
    
    // Set PGP parameters (this also sets PME parameters via inheritance)
    setPGPParameters(alpha, meshSize, potentialCutoff, potentialGridSize, splineOrder, tolerance);
    
    // Initialize PME parameters (PGP inherits from PME)
    pgp_params.setBox(box);
    pgp_params.cutoff = cutoff;
    pgp_params.initializeTables(cutoff);
    pgp_params.initializeBsplines();
    
    // Initialize PGP-specific grids (if not already done by setPGPParameters)
    if (pgp_params.potentialGrid.empty()) {
        pgp_params.initializePotentialGrid();
    }
    
    // Log final parameters
    platform::log(LogLevel::INFO, "PGP parameters initialized: alpha = ", pgp_params.alpha,
                 ", mesh size = [", pgp_params.meshSize[0], ",", pgp_params.meshSize[1], ",", pgp_params.meshSize[2], "]",
                 ", potential grid size = [", pgp_params.potential_grid_size[0], ",", 
                 pgp_params.potential_grid_size[1], ",", pgp_params.potential_grid_size[2], "]",
                 ", spline order = ", pgp_params.splineOrder);
}

void PGPComposite::computeSystemEnergy(model::MCState& state) {
    if (!isInitialized()) {
        platform::log(LogLevel::ERROR, "PGP not initialized before energy calculation");
        return;
    }
    
    // Delegate to the existing implementation
    computeSystemEnergyPGPImpl(state);
}

void PGPComposite::computeMovementEnergy(model::MCState& state) {
    if (!isInitialized()) {
        platform::log(LogLevel::ERROR, "PGP not initialized before movement energy calculation");
        return;
    }
    
    // Delegate to the existing implementation
    computeMovementEnergyPGPImpl(state);
}

bool PGPComposite::validateSetup(const model::MCState& state) {
    if (!isInitialized()) {
        platform::log(LogLevel::WARNING, "PGP parameters not initialized");
        return false;
    }
    
    // Validate system properties
    validateSystemProperties(state);
    
    // Check grid compatibility
    checkGridCompatibility(pgp_params.box);
    
    return true;
}

void PGPComposite::getEnergyBreakdown(const model::MCState& state,
                                    double& realSpace,
                                    double& reciprocal,
                                    double& selfEnergy,
                                    double& vdw,
                                    double& total) {
    // Initialize all components
    realSpace = 0.0;
    reciprocal = 0.0;
    selfEnergy = 0.0;
    vdw = 0.0;
    
    // Get current energy values from state
    realSpace = state.ewald_energy.real_space;
    reciprocal = state.ewald_energy.reciprocal;  // Grid energy is stored in reciprocal field
    selfEnergy = state.ewald_energy.self;
    
    // Calculate VdW energy from active residues
    for (const auto& residue : state.residues) {
        if (residue.active) {
            vdw += residue.energy_vdw;
        }
    }
    
    total = realSpace + reciprocal + selfEnergy + vdw;
}

bool PGPComposite::isInitialized() {
    return pgp_params.initialized;
}

void PGPComposite::reset() {
    pgp_params.initialized = false;
    pgp_params.potentialGrid.clear();
    pgp_params.pmeGrid.clear();
    pgp_params.pmeCharge.clear();
    
    platform::log(LogLevel::INFO, "PGP state reset");
}

void PGPComposite::setDebugMode(bool enable) {
    pgp_params.debug_mode = enable;
    platform::log(LogLevel::INFO, "PGP debug mode ", enable ? "enabled" : "disabled");
}

void PGPComposite::precomputeGrids(model::MCState& state, bool fixedOnly) {
    if (!isInitialized()) {
        platform::log(LogLevel::ERROR, "PGP not initialized before grid precomputation");
        return;
    }
    
    // Delegate to existing precomputation function
    precomputeGridPotential(state, fixedOnly);
}

// Private helper functions
void PGPComposite::validateSystemProperties(const model::MCState& state) {
    // Check for essential system properties
    if (state.residues.empty()) {
        platform::log(LogLevel::WARNING, "No residues in system for PGP calculation");
    }
    
    // Check box dimensions
    for (int i = 0; i < 3; i++) {
        if (pgp_params.box[i] <= 0.0) {
            platform::log(LogLevel::ERROR, "Invalid box dimension ", i, ": ", pgp_params.box[i]);
        }
    }
}

void PGPComposite::checkGridCompatibility(const double box[3]) {
    // Check if grid dimensions are compatible with box size
    for (int i = 0; i < 3; i++) {
        double gridSpacing = box[i] / pgp_params.meshSize[i];
        if (gridSpacing > pgp_params.cutoff / 2.0) {
            platform::log(LogLevel::WARNING, "Grid spacing ", gridSpacing, 
                         " may be too large compared to cutoff ", pgp_params.cutoff);
        }
    }
}

void PGPComposite::initializeGrids() {
    // Initialize potential and PME grids
    pgp_params.initializePotentialGrid();
    
    // Resize PME grids based on mesh size
    size_t gridSize = pgp_params.meshSize[0] * pgp_params.meshSize[1] * pgp_params.meshSize[2];
    pgp_params.pmeGrid.resize(gridSize);
    pgp_params.pmeCharge.resize(gridSize);
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc 