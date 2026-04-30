#include "PGPPrecompute.hpp"
#include "PGPGrid.hpp"
#include "PGPCore.hpp"  // For resetPGPState
#include "PGPGlobal.hpp" // For getPGPParams
#include <cmath>
#include <algorithm>
#include <iostream>

namespace pygcmc {
namespace platform {
namespace cpu {

// Note: setPGPParameters has been moved to PGPCore.cpp

/**
 * @brief Precompute grid potential for fixed parts of the system
 * 
 * This is one of the core functions of the PGP-PME algorithm, responsible for precomputing the electrostatic potential field of the fixed parts of the system.
 * This function assigns charges from the fixed parts to the grid, computes the potential through FFT transformation, and stores the result for later use.
 */
void precomputeGridPotentialImpl(model::MCState& state, bool fixed_only) {
    // Lock for thread safety
    std::lock_guard<std::mutex> lock(pgp_global_mutex);
    
    // Check if parameters are initialized
    if (!getPGPParams().initialized) {
        throw std::runtime_error("PGP parameters not initialized");
    }
    
    // Only output debug information in debug mode
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Starting grid potential precomputation");
        platform::log(LogLevel::DEBUG, "Processing: ", (fixed_only ? "fixed parts only" : "all parts"));
        platform::log(LogLevel::DEBUG, "Grid size: ", getPGPParams().potential_grid_size[0], "x", 
                     getPGPParams().potential_grid_size[1], "x", 
                     getPGPParams().potential_grid_size[2]);
    }
    
    // Backup PME grid and its size, will restore later
    std::vector<std::complex<double>> pmeGridBackup = getPMEParams().pmeGrid;
    size_t originalGridSize = pmeGridBackup.size();
    
    // Store original PME mesh size FIRST before any modifications
    int originalMeshSize[3];
    for (int i = 0; i < 3; i++) {
        originalMeshSize[i] = getPMEParams().meshSize[i];
    }
    
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
    
    // Set PME mesh size to PGP dimensions
    for (int i = 0; i < 3; i++) {
        getPMEParams().meshSize[i] = getPGPParams().potential_grid_size[i];
    }
    
    // Adjust pme_params grid size to fit new grid dimensions
    int totalGridSize = getPGPParams().potential_grid_size[0] * getPGPParams().potential_grid_size[1] * getPGPParams().potential_grid_size[2];
    getPMEParams().pmeGrid.clear();
    getPMEParams().pmeGrid.resize(totalGridSize, std::complex<double>(0.0, 0.0));
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Calling PME charge spreading function (fixed_only=", fixed_only, ")");
    }
    
    spreadChargesOntoGrid(state, fixed_only);
    performFFTForward();
    
    // Apply Ewald factor and convert to potential
    int nx = getPGPParams().potential_grid_size[0];
    int ny = getPGPParams().potential_grid_size[1];
    int nz = getPGPParams().potential_grid_size[2];
    double volume = getPGPParams().box[0] * getPGPParams().box[1] * getPGPParams().box[2];
    
    // Calculate constants needed for Ewald factor application
    double alpha = getPGPParams().alpha;
    double factor = 1.0/(4.0*alpha*alpha);
    
    // Get maximum k vector index
    int maxkx = (nx+1)/2;
    int maxky = (ny+1)/2;
    int maxkz = (nz+1)/2;
    
    // Calculate reciprocal lattice vectors
    double recipBoxVectors[3][3] = {{0}};
    recipBoxVectors[0][0] = 2.0 * M_PI / getPGPParams().box[0]; 
    recipBoxVectors[1][1] = 2.0 * M_PI / getPGPParams().box[1]; 
    recipBoxVectors[2][2] = 2.0 * M_PI / getPGPParams().box[2];
    
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
                std::complex<double> structureFactor = getPMEParams().pmeGrid[index];
                
                // Calculate |k|^2
                double m2 = mhx * mhx + mhy * mhy + mhz * mhz;
                
                // Apply B-spline coefficients
                double bx = getPGPParams().bsplineModuli[0][kx];
                double by = getPGPParams().bsplineModuli[1][ky];
                double bz = getPGPParams().bsplineModuli[2][kz];
                double denom = m2 * bx * by * bz;
                
                // Avoid division by zero problem
                if (denom < 1e-10) {
                    denom = 1e-10;
                }
                
                // Apply k-dependent factor: exp(-k²/(4α²))/(k² · B)
                double kDependentFactor = exp(-m2 * factor) / denom;
                getPMEParams().pmeGrid[index] = structureFactor * kDependentFactor;
            }
        }
    }

    // Remove k=0 mode explicitly before inverse FFT.
    // Keeping the DC component injects a uniform potential offset proportional to total charge.
    if (!getPMEParams().pmeGrid.empty()) {
        getPMEParams().pmeGrid[0] = std::complex<double>(0.0, 0.0);
    }
    
    performFFTBackward();
    
    // Apply normalization and conversion factors
    int totalFFTPoints = nx * ny * nz;
    double constantFactor = 4.0 * M_PI / volume;
    double ONE_4PI_EPS0 = 138.935456; // kJ·mol^-1·nm·e^-2
    double physicalUnitFactor = ONE_4PI_EPS0 / getPGPParams().epsilon_r;
    double totalFactor = totalFFTPoints * constantFactor * physicalUnitFactor;
    
    // Apply total correction factor to each grid point
    for (int i = 0; i < totalGridSize; i++) {
        getPMEParams().pmeGrid[i] *= totalFactor;
    }
    
    // Safe copy of PME grid to PGP's potentialGrid
    size_t gridSize = getPMEParams().meshSize[0] * getPMEParams().meshSize[1] * getPMEParams().meshSize[2];
    
    // Clear existing grid safely
    getPGPParams().potentialGrid.clear();
    getPGPParams().potentialGrid.shrink_to_fit();
    
    // Reserve and copy
    getPGPParams().potentialGrid.reserve(gridSize);
    getPGPParams().potentialGrid.assign(getPMEParams().pmeGrid.begin(), getPMEParams().pmeGrid.end());
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Copied ", gridSize, " grid points to potential grid");
    }

    // Restore original PME mesh size FIRST
    for (int i = 0; i < 3; i++) {
        getPMEParams().meshSize[i] = originalMeshSize[i];
    }
    
    // Then restore original PME grid with correct size
    getPMEParams().pmeGrid = std::move(pmeGridBackup);
    
    // Verify grid size is restored correctly
    if (getPMEParams().pmeGrid.size() != originalGridSize) {
        platform::log(LogLevel::ERROR, "PME grid size mismatch after restore: expected ", 
                     originalGridSize, " but got ", getPMEParams().pmeGrid.size());
        // Force correct size
        getPMEParams().pmeGrid.resize(originalGridSize, std::complex<double>(0.0, 0.0));
    }
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Potential precomputation completed");
    }
}

// Public interface wrappers
void precomputeGridPotential(model::MCState& state, bool fixed_only) {
    // Removed static call count to avoid thread safety and memory issues
    // The cleanup mechanism should be handled by the user calling _cleanup() explicitly
    
    precomputeGridPotentialImpl(state, fixed_only);
}

// <agent-hook:pgp_precompute_impl>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 
