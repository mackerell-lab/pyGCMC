#include "PGPPrecompute.hpp"
#include "PGPCore.hpp"
#include "PGPGlobal.hpp"
#include "platform/cpu/energy/pme/PMEComposite.hpp"
#include "platform/cpu/energy/pme/PMEGlobal.hpp"
#include "platform/cpu/energy/pme/PMEGrid.hpp"
#include "platform/cpu/energy/pme/PMEFFTCore.hpp"
#include "platform/platform.hpp"
#include <cmath>
#include <mutex>

namespace pygcmc {
namespace platform {
namespace cpu {

// Thread-local storage for temporary PME state
thread_local struct {
    std::vector<std::complex<double>> tempPmeGrid;
    int tempMeshSize[3] = {0, 0, 0};
    bool inUse = false;
} tls_pme_state;

// Mutex for protecting global state modifications
static std::mutex precompute_mutex;

/**
 * @brief Safely precompute grid potential without modifying global PME state
 */
void safePrecomputeGridPotentialImpl(model::MCState& state, bool fixed_only) {
    std::lock_guard<std::mutex> lock(precompute_mutex);
    
    // Check if we're already in a precompute operation (prevent recursion)
    if (tls_pme_state.inUse) {
        platform::log(LogLevel::ERROR, "Recursive precompute detected, aborting");
        return;
    }
    
    // Mark thread-local state as in use
    tls_pme_state.inUse = true;
    
    // Store original PME state
    int originalMeshSize[3];
    for (int i = 0; i < 3; i++) {
        originalMeshSize[i] = getPMEParams().meshSize[i];
        tls_pme_state.tempMeshSize[i] = getPGPParams().potential_grid_size[i];
    }
    
    // Save original PME grid
    auto originalPmeGrid = std::move(getPMEParams().pmeGrid);
    
    try {
        // Temporarily set PME parameters for PGP calculation
        for (int i = 0; i < 3; i++) {
            getPMEParams().meshSize[i] = getPGPParams().potential_grid_size[i];
        }
        
        // Create new PME grid with PGP dimensions
        int totalGridSize = getPGPParams().potential_grid_size[0] * 
                           getPGPParams().potential_grid_size[1] * 
                           getPGPParams().potential_grid_size[2];
        getPMEParams().pmeGrid.clear();
        getPMEParams().pmeGrid.resize(totalGridSize, std::complex<double>(0.0, 0.0));
        
        // Perform PME calculation
        spreadChargesOntoGrid(state, fixed_only);
        performFFTForward();
        
        // Apply Ewald factor and convert to potential
        int nx = getPGPParams().potential_grid_size[0];
        int ny = getPGPParams().potential_grid_size[1];
        int nz = getPGPParams().potential_grid_size[2];
        double volume = getPGPParams().box[0] * getPGPParams().box[1] * getPGPParams().box[2];
        
        double alpha = getPGPParams().alpha;
        double factor = 1.0/(4.0*alpha*alpha);
        
        int maxkx = (nx+1)/2;
        int maxky = (ny+1)/2;
        int maxkz = (nz+1)/2;
        
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
                    if (kx == 0 && ky == 0 && kz == 0) {
                        continue;
                    }
                    
                    double mz = (kz < maxkz) ? kz : (kz-nz);
                    double mhz = mz * recipBoxVectors[2][2];
                    
                    int index = kx * ny * nz + ky * nz + kz;
                    std::complex<double> structureFactor = getPMEParams().pmeGrid[index];
                    
                    double m2 = mhx * mhx + mhy * mhy + mhz * mhz;
                    double expterm = std::exp(-factor * m2) / m2;
                    
                    getPMEParams().pmeGrid[index] = structureFactor * 
                        std::complex<double>(4.0 * M_PI * expterm / volume, 0.0);
                }
            }
        }
        
        // Perform inverse FFT
        performFFTBackward();
        
        // Copy result to PGP potential grid (with proper synchronization)
        getPGPParams().potentialGrid.clear();
        getPGPParams().potentialGrid.reserve(totalGridSize);
        getPGPParams().potentialGrid.assign(getPMEParams().pmeGrid.begin(), 
                                           getPMEParams().pmeGrid.end());
        
        // Restore original PME state
        getPMEParams().pmeGrid = std::move(originalPmeGrid);
        for (int i = 0; i < 3; i++) {
            getPMEParams().meshSize[i] = originalMeshSize[i];
        }
        
    } catch (...) {
        // Ensure we restore PME state even if exception occurs
        getPMEParams().pmeGrid = std::move(originalPmeGrid);
        for (int i = 0; i < 3; i++) {
            getPMEParams().meshSize[i] = originalMeshSize[i];
        }
        tls_pme_state.inUse = false;
        throw;
    }
    
    // Clear thread-local state
    tls_pme_state.inUse = false;
    
    platform::log(LogLevel::DEBUG, "Grid potential precomputation completed safely");
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc