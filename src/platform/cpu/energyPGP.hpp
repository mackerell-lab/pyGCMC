#pragma once

#include "model/montecarlo.hpp"
#include "platform/platform.hpp"
#include "energyCommon.hpp"
#include <vector>
#include <complex>

namespace pygcmc {
namespace platform {
namespace cpu {

// Structure for PGP (Pair Grid PME) parameters
struct PGPParams {
    // Basic PME parameters
    double alpha;               // Ewald splitting parameter
    double cutoff;              // Cutoff distance
    double tolerance;           // Error tolerance
    double epsilon_r;           // Relative permittivity
    int splineOrder;            // B-spline order
    double box[3];              // Box dimensions
    int meshSize[3];            // FFT mesh size
    bool initialized = false;   // Whether parameters are initialized
    
    // Additional pair grid parameters
    double pair_cutoff;         // Cutoff for pair interactions
    int pair_grid_size[3];      // Size of the pair grid
    double grid_spacing;        // Grid spacing for pair grid
    
    // Lookup tables (same as PME)
    std::vector<double> erfcTable;
    std::vector<double> ewaldScaleTable;
    double ewaldDX;
    double ewaldDXInv;
    double erfcDXInv;
    
    // B-spline moduli (same as PME)
    std::vector<std::vector<double>> bsplineModuli;
    
    // Reciprocal space grid (same as PME)
    std::vector<std::complex<double>> pmeGrid;
    std::vector<double> pmeCharge;
    
    // Pair grid for PP part
    std::vector<std::complex<double>> pairGrid;
    
    // Initialize the pair grid
    void initializePairGrid();
};

// Global PGP parameters
extern PGPParams pgp_params;

// Set PGP parameters
void setPGPParameters(double alpha, const int meshSize[3], double pair_cutoff, 
                        const int pairGridSize[3], int splineOrder = 4, double tolerance = 1e-5);

// Auto-adjust PGP parameters based on error tolerance and box size
void autoAdjustPGPParameters(double error_tolerance, double cutoff_distance, 
                               double pair_cutoff, const double box[3]);

// Spread pairs onto grid for PP part
void spreadPairsOntoGrid(model::MCState& state, bool movement_only);

// Compute energy from pair grid
void computeEnergyFromPairGrid(double& energy);

// Compute reciprocal space energy using PGP
double computeReciprocalPGP(model::MCState& state, bool movement_only);

// Compute self energy
double computeSelfEnergyPGP(model::MCState& state, bool movement_only);

// Compute real space energy using the pair grid approach
void computeRealSpacePGP(model::MCState& state, bool movement_only, bool store_in_residues);

// Compute pair grid interactions
void computePairGridPGP(model::MCState& state, bool movement_only);

// Main system energy calculation function
void computeSystemEnergyPGP(model::MCState& state);

// Movement energy calculation (only affected residues)
void computeMovementEnergyPGP(model::MCState& state);

} // namespace cpu
} // namespace platform
} // namespace pygcmc 