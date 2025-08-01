#pragma once

#include "model/ModelModule.hpp"
#include "platform/platform.hpp"
#include "../common/EnergyConstants.hpp"
#include "../common/EnergyUtils.hpp"
#include "../pme/PMEComposite.hpp"
#include <array>
#include <vector>
#include <complex>
#include <mutex>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief PGP-PME algorithm parameter structure (Precomputed Grid-Potential Particle Mesh Ewald)
 * 
 * This structure contains all parameters and data structures required by the Precomputed Grid-Potential Particle Mesh Ewald algorithm.
 * PGP-PME is an optimized PME method that accelerates energy evaluation in Monte Carlo simulations by precomputing potential grids.
 */
struct PGPParams : public PMEParams {
    // PGP-specific parameters (not inherited from PMEParams)
    double potential_cutoff;              // Cutoff distance for potential calculation
    int potential_grid_size[3];           // Precomputed potential grid dimensions
    double grid_spacing;                  // Grid spacing
    std::vector<std::complex<double>> potentialGrid;  // Precomputed potential grid data
    
    // Debug flags
    bool debug_mode = true;  // Debug mode enabled by default

    /**
     * @brief Initialize the 3D grid for precomputed potential
     */
    void initializePotentialGrid();
};

// Note: Global PGP parameters are now managed by PGPGlobal.hpp
// Use getPGPParams() to access the parameters

// Core function declarations
void setPGPParameters(double alpha, const int meshSize[3], double potential_cutoff, 
                        const int potentialGridSize[3], int splineOrder, double tolerance);

void precomputeGridPotential(model::MCState& state, bool fixed_only = true);

void interpolateMoleculeEnergy(model::MCState& state, double& energy);

double calculateMoleculeEnergy(model::MCState& state);

// Reset function to clear global state - fixes memory corruption bug
void resetPGPState();

double computeMoleculeEnergyGlobal(model::MCState& state, const std::vector<int>& movementResidues, 
                                   const std::vector<int>& nearbyResidues, int threadIndex);

void computeRealSpacePGP(model::MCState& state, bool movement_only, bool store_in_residues = true);

double computeSelfEnergyPGP(model::MCState& state, bool movement_only);

void computeSystemEnergyPGP(model::MCState& state);

void computeMovementEnergyPGP(model::MCState& state);

// Fixed versions that correctly use pgp_params for erfc calculations
void computeSystemEnergyPGPFixed(model::MCState& state);

void computeMovementEnergyPGPFixed(model::MCState& state);

// <agent-hook:pgp_core>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 