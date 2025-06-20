#include "PGPPrecompute.hpp"
#include "PGPGrid.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Set all parameters for the PGP-PME algorithm
 * 
 * This function is the entry point for the PGP-PME algorithm, used to configure all parameters needed for the algorithm to run.
 * It first sets standard PME parameters, then adds PGP-specific parameters, and finally initializes the potential grid.
 */
void setPGPParametersImpl(double alpha, const int meshSize[3], double potential_cutoff, 
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
    
    // Only output parameter setting information in debug mode
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "PGP parameters set: alpha=", alpha, 
                    ", potential_cutoff=", potential_cutoff, 
                    ", potentialGrid=[", potentialGridSize[0], ",", potentialGridSize[1], ",", potentialGridSize[2], "]");
    }
}

/**
 * @brief Precompute grid potential for fixed parts of the system
 * 
 * This is one of the core functions of the PGP-PME algorithm, responsible for precomputing the electrostatic potential field of the fixed parts of the system.
 * This function assigns charges from the fixed parts to the grid, computes the potential through FFT transformation, and stores the result for later use.
 */
void precomputeGridPotentialImpl(model::MCState& state, bool fixed_only) {
    // Check if parameters are initialized
    if (!pgp_params.initialized) {
        throw std::runtime_error("PGP parameters not initialized");
    }
    
    // Only output debug information in debug mode
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Starting grid potential precomputation");
        platform::log(LogLevel::DEBUG, "Processing: ", (fixed_only ? "fixed parts only" : "all parts"));
        platform::log(LogLevel::DEBUG, "Grid size: ", pgp_params.potential_grid_size[0], "x", 
                     pgp_params.potential_grid_size[1], "x", 
                     pgp_params.potential_grid_size[2]);
    }
    
    // Backup PME grid, will restore later
    std::vector<std::complex<double>> pmeGridBackup = pme_params.pmeGrid;
    
    // Reset PME grid, prepare for new calculation
    std::fill(pme_params.pmeGrid.begin(), pme_params.pmeGrid.end(), std::complex<double>(0.0, 0.0));
    
    // Statistics - only calculated in debug mode
    int fixed_residues_count = 0;
    
    // Calculate number of fixed residues - only computed in debug mode or when checking fixed_only validity
    if (platform::is_debug_mode() || fixed_only) {
        for (int i = 0; i < state.activeResidueCount; ++i) {
            const auto& res = state.residues[i];
            if (res.fixed && res.active) fixed_residues_count++;
        }
        
        if (platform::is_debug_mode()) {
            platform::log(LogLevel::DEBUG, "Number of fixed residues: ", fixed_residues_count);
        }
    }
    
    // If no fixed residues but requesting fixed_only, issue warning and switch automatically
    if (fixed_only && fixed_residues_count == 0) {
        platform::log(LogLevel::WARNING, "No fixed residues found, switching to process all atoms");
        fixed_only = false;
    }
    
    // Set pme_params grid size to match pgp_params grid size
    for (int i = 0; i < 3; i++) {
        pme_params.meshSize[i] = pgp_params.potential_grid_size[i];
    }
    
    // Adjust pme_params grid size to fit new grid dimensions
    int totalGridSize = pgp_params.potential_grid_size[0] * pgp_params.potential_grid_size[1] * pgp_params.potential_grid_size[2];
    pme_params.pmeGrid.resize(totalGridSize, std::complex<double>(0.0, 0.0));
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Calling PME charge spreading function (fixed_only=", fixed_only, ")");
    }
    
    spreadChargesOntoGrid(state, fixed_only);
    performFFTForward();
    
    // Apply Ewald factor and convert to potential
    int nx = pgp_params.potential_grid_size[0];
    int ny = pgp_params.potential_grid_size[1];
    int nz = pgp_params.potential_grid_size[2];
    double volume = pgp_params.box[0] * pgp_params.box[1] * pgp_params.box[2];
    
    // Calculate constants needed for Ewald factor application
    double alpha = pgp_params.alpha;
    double factor = 1.0/(4.0*alpha*alpha);
    
    // Get maximum k vector index
    int maxkx = (nx+1)/2;
    int maxky = (ny+1)/2;
    int maxkz = (nz+1)/2;
    
    // Calculate reciprocal lattice vectors
    double recipBoxVectors[3][3] = {{0}};
    recipBoxVectors[0][0] = 2.0 * M_PI / pgp_params.box[0]; 
    recipBoxVectors[1][1] = 2.0 * M_PI / pgp_params.box[1]; 
    recipBoxVectors[2][2] = 2.0 * M_PI / pgp_params.box[2];
    
    // Apply Ewald factor
    for (int kx = 0; kx < nx; kx++) {
        double mx = (kx < maxkx) ? kx : (kx-nx);
        double mhx = mx * recipBoxVectors[0][0];
        
        for (int ky = 0; ky < ny; ky++) {
            double my = (ky < maxky) ? ky : (ky-ny);
            double mhy = my * recipBoxVectors[1][1];
            
            for (int kz = 0; kz < nz; kz++) {
                // Skip zero frequency
                if (kx == 0 && ky == 0 && kz == 0) {
                    continue;
                }
                
                double mz = (kz < maxkz) ? kz : (kz-nz);
                double mhz = mz * recipBoxVectors[2][2];
                
                // Grid index
                int index = kx * ny * nz + ky * nz + kz;
                std::complex<double> structureFactor = pme_params.pmeGrid[index];
                
                // Calculate |k|^2
                double m2 = mhx * mhx + mhy * mhy + mhz * mhz;
                
                // Apply B-spline coefficients
                double bx = pgp_params.bsplineModuli[0][kx];
                double by = pgp_params.bsplineModuli[1][ky];
                double bz = pgp_params.bsplineModuli[2][kz];
                double denom = m2 * bx * by * bz;
                
                // Avoid division by zero problem
                if (denom < 1e-10) {
                    denom = 1e-10;
                }
                
                // Apply k-dependent factor: exp(-k²/(4α²))/(k² · B)
                double kDependentFactor = exp(-m2 * factor) / denom;
                pme_params.pmeGrid[index] = structureFactor * kDependentFactor;
            }
        }
    }
    
    performFFTBackward();
    
    // Apply normalization and conversion factors
    int totalFFTPoints = nx * ny * nz;
    double constantFactor = 4.0 * M_PI / volume;
    double ONE_4PI_EPS0 = 138.935456; // kJ·mol^-1·nm·e^-2
    double physicalUnitFactor = ONE_4PI_EPS0 / pgp_params.epsilon_r;
    double totalFactor = totalFFTPoints * constantFactor * physicalUnitFactor * 0.5;
    
    // Apply total correction factor to each grid point
    for (int i = 0; i < totalGridSize; i++) {
        pme_params.pmeGrid[i] *= totalFactor;
    }
    
    // Copy modified PME grid to PGP's potentialGrid
    pgp_params.potentialGrid = pme_params.pmeGrid;

    // Restore original PME grid
    pme_params.pmeGrid = pmeGridBackup;
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Potential precomputation completed");
    }
}

// Public interface wrappers
void setPGPParameters(double alpha, const int meshSize[3], double potential_cutoff, 
                        const int potentialGridSize[3], int splineOrder, double tolerance) {
    setPGPParametersImpl(alpha, meshSize, potential_cutoff, potentialGridSize, splineOrder, tolerance);
}

void precomputeGridPotential(model::MCState& state, bool fixed_only) {
    precomputeGridPotentialImpl(state, fixed_only);
}

// <agent-hook:pgp_precompute_impl>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 