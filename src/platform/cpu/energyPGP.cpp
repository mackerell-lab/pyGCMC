#include "energyPGP.hpp"
#include "energyPME.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace pygcmc {
namespace platform {
namespace cpu {

// Initialize global PGP parameters
PGPParams pgp_params;

// Initialize pair grid
void PGPParams::initializePairGrid() {
    // Calculate grid size and allocate memory for pair grid
    int totalSize = pair_grid_size[0] * pair_grid_size[1] * pair_grid_size[2];
    pairGrid.resize(totalSize);
    
    // Set grid spacing based on box dimensions and grid size
    grid_spacing = std::min({
        box[0] / pair_grid_size[0],
        box[1] / pair_grid_size[1],
        box[2] / pair_grid_size[2]
    });
    
    platform::log(LogLevel::DEBUG, "PGP pair grid initialized with size: ", 
                 pair_grid_size[0], "x", pair_grid_size[1], "x", pair_grid_size[2],
                 ", grid spacing: ", grid_spacing);
}

// Set PGP parameters
void setPGPParameters(double alpha, const int meshSize[3], double pair_cutoff, 
                        const int pairGridSize[3], int splineOrder, double tolerance) {
    // Set PME parameters (reuses PME functionality)
    setPMEParameters(alpha, meshSize, splineOrder, tolerance);
    
    // Copy PME parameters to PGP params
    pgp_params.alpha = pme_params.alpha;
    pgp_params.tolerance = pme_params.tolerance;
    pgp_params.initialized = pme_params.initialized;
    pgp_params.cutoff = pme_params.cutoff;
    pgp_params.epsilon_r = pme_params.epsilon_r;
    pgp_params.splineOrder = pme_params.splineOrder;
    
    // Copy box dimensions
    for (int i = 0; i < 3; i++) {
        pgp_params.box[i] = pme_params.box[i];
        pgp_params.meshSize[i] = pme_params.meshSize[i];
    }
    
    // Copy lookup tables
    pgp_params.erfcTable = pme_params.erfcTable;
    pgp_params.ewaldScaleTable = pme_params.ewaldScaleTable;
    pgp_params.ewaldDX = pme_params.ewaldDX;
    pgp_params.ewaldDXInv = pme_params.ewaldDXInv;
    pgp_params.erfcDXInv = pme_params.erfcDXInv;
    
    // Copy B-spline moduli
    for (int i = 0; i < 3; i++) {
        pgp_params.bsplineModuli[i] = pme_params.bsplineModuli[i];
    }
    
    // Copy reciprocal space grid
    pgp_params.pmeGrid = pme_params.pmeGrid;
    pgp_params.pmeCharge = pme_params.pmeCharge;
    
    // Set PGP-specific parameters
    pgp_params.pair_cutoff = pair_cutoff;
    for (int i = 0; i < 3; i++) {
        pgp_params.pair_grid_size[i] = pairGridSize[i];
    }
    
    platform::log(LogLevel::INFO, "PGP parameters set: alpha=", alpha, 
                 ", pair_cutoff=", pair_cutoff, 
                 ", pairGrid=[", pairGridSize[0], ",", pairGridSize[1], ",", pairGridSize[2], "]");
}

// Auto-adjust PGP parameters
void autoAdjustPGPParameters(double error_tolerance, double cutoff_distance, 
                               double pair_cutoff, const double box[3]) {
    // First, auto-adjust the PME parameters
    autoAdjustPMEParameters(error_tolerance, cutoff_distance, box);
    
    // Copy PME parameters to PGP params (same as in setPGPParameters)
    pgp_params.alpha = pme_params.alpha;
    pgp_params.tolerance = pme_params.tolerance;
    pgp_params.initialized = pme_params.initialized;
    pgp_params.cutoff = pme_params.cutoff;
    pgp_params.epsilon_r = pme_params.epsilon_r;
    pgp_params.splineOrder = pme_params.splineOrder;
    
    // Copy box dimensions
    for (int i = 0; i < 3; i++) {
        pgp_params.box[i] = pme_params.box[i];
        pgp_params.meshSize[i] = pme_params.meshSize[i];
    }
    
    // Copy lookup tables
    pgp_params.erfcTable = pme_params.erfcTable;
    pgp_params.ewaldScaleTable = pme_params.ewaldScaleTable;
    pgp_params.ewaldDX = pme_params.ewaldDX;
    pgp_params.ewaldDXInv = pme_params.ewaldDXInv;
    pgp_params.erfcDXInv = pme_params.erfcDXInv;
    
    // Copy B-spline moduli
    for (int i = 0; i < 3; i++) {
        pgp_params.bsplineModuli[i] = pme_params.bsplineModuli[i];
    }
    
    // Copy reciprocal space grid
    pgp_params.pmeGrid = pme_params.pmeGrid;
    pgp_params.pmeCharge = pme_params.pmeCharge;
    
    // Set PGP-specific parameters
    pgp_params.pair_cutoff = pair_cutoff;
    
    // Calculate appropriate pair grid size based on box and pair cutoff
    for (int i = 0; i < 3; i++) {
        // Simple heuristic: ratio of box size to cutoff, with minimum size
        pgp_params.pair_grid_size[i] = std::max(32, static_cast<int>(box[i] / pair_cutoff * 2.0));
        // Ensure it's a power of 2 for FFT efficiency
        pgp_params.pair_grid_size[i] = 1 << static_cast<int>(std::ceil(std::log2(pgp_params.pair_grid_size[i])));
    }
    
    // Initialize the pair grid
    pgp_params.initializePairGrid();
    
    platform::log(LogLevel::INFO, "PGP parameters auto-adjusted: ", 
                 "alpha=", pgp_params.alpha, 
                 ", pair_cutoff=", pgp_params.pair_cutoff, 
                 ", pairGrid=[", pgp_params.pair_grid_size[0], ",", 
                 pgp_params.pair_grid_size[1], ",", pgp_params.pair_grid_size[2], "]");
}

// Spread pairs onto grid for PP part
void spreadPairsOntoGrid([[maybe_unused]] model::MCState& state, [[maybe_unused]] bool movement_only) {
    // Reset the pair grid
    std::fill(pgp_params.pairGrid.begin(), pgp_params.pairGrid.end(), std::complex<double>(0.0, 0.0));
    
    // Implementation will depend on the specific PGP algorithm
    platform::log(LogLevel::DEBUG, "Spreading pairs onto grid. Movement only: ", movement_only);
    
    // TODO: Implement the pair grid spreading algorithm
}

// Compute energy from pair grid
void computeEnergyFromPairGrid([[maybe_unused]] double& energy) {
    // Implementation will depend on the specific PGP algorithm
    platform::log(LogLevel::DEBUG, "Computing energy from pair grid");
    
    // TODO: Implement energy calculation from pair grid
}

// Compute reciprocal space energy using PGP
double computeReciprocalPGP(model::MCState& state, bool movement_only) {
    // Reuse the PME reciprocal calculation for the long-range part
    double reciprocal_energy = computeReciprocalPME(state, movement_only);
    
    // Add the pair-grid contribution
    spreadPairsOntoGrid(state, movement_only);
    double pair_grid_energy = 0.0;
    computeEnergyFromPairGrid(pair_grid_energy);
    
    platform::log(LogLevel::DEBUG, "PGP reciprocal energy: PME=", reciprocal_energy, 
                 ", pair-grid=", pair_grid_energy, 
                 ", total=", reciprocal_energy + pair_grid_energy);
    
    return reciprocal_energy + pair_grid_energy;
}

// Compute self energy (mostly the same as PME)
double computeSelfEnergyPGP(model::MCState& state, bool movement_only) {
    // The self energy calculation is the same as in PME
    return computeSelfEnergyPME(state, movement_only);
}

// Compute real space energy using the pair grid approach
void computeRealSpacePGP(model::MCState& state, bool movement_only, bool store_in_residues) {
    // This is a modified version of the real space calculation
    // TODO: Implement the PGP real space calculation
    
    // For now, just use the PME real space calculation
    computeRealSpacePME(state, movement_only, store_in_residues);
    
    platform::log(LogLevel::DEBUG, "PGP real space energy calculated");
}

// Compute pair grid interactions
void computePairGridPGP([[maybe_unused]] model::MCState& state, [[maybe_unused]] bool movement_only) {
    // This is a new function specific to PGP
    // TODO: Implement the pair grid calculations
    
    platform::log(LogLevel::DEBUG, "PGP pair grid calculation");
}

// Main system energy calculation function
void computeSystemEnergyPGP(model::MCState& state) {
    if (!pgp_params.initialized) {
        throw std::runtime_error("PGP parameters not initialized");
    }
    
    // Check system neutrality
    checkSystemNeutrality(state);
    
    // Reset the ewald_energy structure
    state.ewald_energy = {0.0, 0.0, 0.0};
    
    // Compute real space energy (includes vdw)
    computeRealSpacePGP(state, false, true);
    
    // Compute reciprocal space energy
    state.ewald_energy.reciprocal = computeReciprocalPGP(state, false);
    
    // Compute self energy
    state.ewald_energy.self = computeSelfEnergyPGP(state, false);
    
    platform::log(LogLevel::DEBUG, "PGP system energy: real=", state.ewald_energy.real_space,
                 ", reciprocal=", state.ewald_energy.reciprocal,
                 ", self=", state.ewald_energy.self,
                 ", total=", state.ewald_energy.real_space + state.ewald_energy.reciprocal + state.ewald_energy.self);
}

// Movement energy calculation (only affected residues)
void computeMovementEnergyPGP(model::MCState& state) {
    if (!pgp_params.initialized) {
        throw std::runtime_error("PGP parameters not initialized");
    }
    
    // Check system neutrality
    checkSystemNeutrality(state);
    
    // Reset the ewald_energy structure
    state.ewald_energy = {0.0, 0.0, 0.0};
    
    // Compute real space energy for movement residues
    computeRealSpacePGP(state, true, true);
    
    // Compute reciprocal space energy for movement
    state.ewald_energy.reciprocal = computeReciprocalPGP(state, true);
    
    // Compute self energy for movement
    state.ewald_energy.self = computeSelfEnergyPGP(state, true);
    
    platform::log(LogLevel::DEBUG, "PGP movement energy: real=", state.ewald_energy.real_space,
                 ", reciprocal=", state.ewald_energy.reciprocal,
                 ", self=", state.ewald_energy.self,
                 ", total=", state.ewald_energy.real_space + state.ewald_energy.reciprocal + state.ewald_energy.self);
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc 