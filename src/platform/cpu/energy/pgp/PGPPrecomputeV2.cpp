#include "PGPPrecompute.hpp"
#include "PGPCore.hpp"
#include "PGPGlobal.hpp"
#include "platform/cpu/energy/pme/PMEComposite.hpp"
#include "platform/cpu/energy/pme/PMEGlobal.hpp"
#include "platform/cpu/energy/pme/PMEGrid.hpp"
#include "platform/cpu/energy/pme/PMEFFTCore.hpp"
#include "platform/platform.hpp"
#include <cmath>
#include <algorithm>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Thread-local workspace for PGP calculations
 */
struct PGPWorkspace {
    std::vector<std::complex<double>> tempGrid;
    std::vector<double> tempBsplineModuli[3];
    int currentGridSize[3] = {0, 0, 0};
    
    void ensureSize(int nx, int ny, int nz) {
        int newSize = nx * ny * nz;
        if (tempGrid.size() < static_cast<size_t>(newSize)) {
            tempGrid.resize(newSize);
        }
        // Clear the grid
        std::fill(tempGrid.begin(), tempGrid.begin() + newSize, std::complex<double>(0.0, 0.0));
        
        currentGridSize[0] = nx;
        currentGridSize[1] = ny;
        currentGridSize[2] = nz;
    }
};

// Thread-local workspace to avoid global state modifications
thread_local PGPWorkspace workspace;

/**
 * @brief Spread charges onto a local grid without modifying global state
 */
void spreadChargesOntoLocalGrid(const model::MCState& state, bool fixed_only,
                               std::vector<std::complex<double>>& grid,
                               const int gridSize[3]) {
    // Zero out the grid
    std::fill(grid.begin(), grid.end(), std::complex<double>(0.0, 0.0));
    
    // Process atoms based on fixed_only flag
    for (int res_idx = 0; res_idx < state.activeResidueCount; ++res_idx) {
        const auto& res = state.residues[res_idx];
        
        // Skip if not active or if we only want fixed and this isn't fixed
        if (!res.active || (fixed_only && !res.fixed)) {
            continue;
        }
        
        // Process atoms in this residue
        for (int i = 0; i < res.atomCount; ++i) {
            int atom_idx = res.atomStart + i;
            const auto& atom = state.atoms[atom_idx];
            
            // Calculate fractional coordinates
            double fx = atom.x / getPGPParams().box[0];
            double fy = atom.y / getPGPParams().box[1];
            double fz = atom.z / getPGPParams().box[2];
            
            // Wrap to [0,1)
            fx = fx - std::floor(fx);
            fy = fy - std::floor(fy);
            fz = fz - std::floor(fz);
            
            // Convert to grid coordinates
            double gx = fx * gridSize[0];
            double gy = fy * gridSize[1];
            double gz = fz * gridSize[2];
            
            // Find base grid point
            int ix0 = static_cast<int>(std::floor(gx));
            int iy0 = static_cast<int>(std::floor(gy));
            int iz0 = static_cast<int>(std::floor(gz));
            
            // Calculate B-spline weights
            double dx = gx - ix0;
            double dy = gy - iy0;
            double dz = gz - iz0;
            
            // Simple 4th order B-spline weights (simplified)
            double wx[4], wy[4], wz[4];
            
            // Calculate weights for each dimension
            for (int order = 0; order < 4; order++) {
                double t = dx - order + 1.5;
                wx[order] = (t > 0 && t < 1) ? (1 - t) : 0;
                
                t = dy - order + 1.5;
                wy[order] = (t > 0 && t < 1) ? (1 - t) : 0;
                
                t = dz - order + 1.5;
                wz[order] = (t > 0 && t < 1) ? (1 - t) : 0;
            }
            
            // Spread charge to grid
            for (int ix = 0; ix < 4; ix++) {
                int gix = (ix0 + ix - 1 + gridSize[0]) % gridSize[0];
                for (int iy = 0; iy < 4; iy++) {
                    int giy = (iy0 + iy - 1 + gridSize[1]) % gridSize[1];
                    for (int iz = 0; iz < 4; iz++) {
                        int giz = (iz0 + iz - 1 + gridSize[2]) % gridSize[2];
                        
                        int idx = gix * gridSize[1] * gridSize[2] + giy * gridSize[2] + giz;
                        double weight = wx[ix] * wy[iy] * wz[iz];
                        grid[idx] += std::complex<double>(atom.charge * weight, 0.0);
                    }
                }
            }
        }
    }
}

/**
 * @brief Safe precompute that doesn't modify global state
 */
void precomputeGridPotentialSafe(model::MCState& state, bool fixed_only) {
    // Lock for thread safety
    static std::mutex local_mutex;
    std::lock_guard<std::mutex> lock(local_mutex);
    
    // Check if parameters are initialized
    if (!getPGPParams().initialized) {
        throw std::runtime_error("PGP parameters not initialized");
    }
    
    // Get grid dimensions
    int nx = getPGPParams().potential_grid_size[0];
    int ny = getPGPParams().potential_grid_size[1];
    int nz = getPGPParams().potential_grid_size[2];
    int totalSize = nx * ny * nz;
    
    // Ensure workspace is large enough
    workspace.ensureSize(nx, ny, nz);
    
    // Spread charges onto local grid
    spreadChargesOntoLocalGrid(state, fixed_only, workspace.tempGrid, 
                              getPGPParams().potential_grid_size);
    
    // Perform FFT on local grid
    CustomFFT::padded_fft(workspace.tempGrid.data(), totalSize, false);
    
    // Apply Ewald factor
    double alpha = getPGPParams().alpha;
    double factor = 1.0/(4.0*alpha*alpha);
    double volume = getPGPParams().box[0] * getPGPParams().box[1] * getPGPParams().box[2];
    
    int maxkx = (nx+1)/2;
    int maxky = (ny+1)/2;
    int maxkz = (nz+1)/2;
    
    double recipBoxVectors[3][3] = {{0}};
    recipBoxVectors[0][0] = 2.0 * M_PI / getPGPParams().box[0];
    recipBoxVectors[1][1] = 2.0 * M_PI / getPGPParams().box[1];
    recipBoxVectors[2][2] = 2.0 * M_PI / getPGPParams().box[2];
    
    // Apply Ewald factor in Fourier space
    for (int kx = 0; kx < nx; kx++) {
        double mx = (kx < maxkx) ? kx : (kx-nx);
        double mhx = mx * recipBoxVectors[0][0];
        
        for (int ky = 0; ky < ny; ky++) {
            double my = (ky < maxky) ? ky : (ky-ny);
            double mhy = my * recipBoxVectors[1][1];
            
            for (int kz = 0; kz < nz; kz++) {
                if (kx == 0 && ky == 0 && kz == 0) {
                    continue;
                }
                
                double mz = (kz < maxkz) ? kz : (kz-nz);
                double mhz = mz * recipBoxVectors[2][2];
                
                int index = kx * ny * nz + ky * nz + kz;
                double m2 = mhx * mhx + mhy * mhy + mhz * mhz;
                
                // Simple Ewald factor without B-spline correction
                double ewaldFactor = exp(-m2 * factor) / m2;
                workspace.tempGrid[index] *= ewaldFactor;
            }
        }
    }
    
    // Perform inverse FFT
    CustomFFT::padded_fft(workspace.tempGrid.data(), totalSize, true);
    
    // Apply normalization
    double ONE_4PI_EPS0 = 138.935456;
    double constantFactor = 4.0 * M_PI / volume;
    double physicalUnitFactor = ONE_4PI_EPS0 / getPGPParams().epsilon_r;
    double totalFactor = totalSize * constantFactor * physicalUnitFactor;
    
    for (int i = 0; i < totalSize; i++) {
        workspace.tempGrid[i] *= totalFactor;
    }
    
    // Copy to PGP potential grid with proper locking
    getPGPParams().potentialGrid.clear();
    getPGPParams().potentialGrid.reserve(totalSize);
    getPGPParams().potentialGrid.assign(workspace.tempGrid.begin(), 
                                       workspace.tempGrid.begin() + totalSize);
    
    platform::log(LogLevel::DEBUG, "Grid potential precomputation completed safely");
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc