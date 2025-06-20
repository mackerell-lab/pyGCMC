#include "PMEEnergyCalc.hpp"
#include "PMECore.hpp"
#include "PMEFFT.hpp"
#include "PMEGridCharge.hpp"
#include "PMEReciprocal.hpp"
#include "platform/platform.hpp"
#include <cmath>
#include <algorithm>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Compute energy from the PME grid after FFT
 * 
 * @param energy Output energy
 */
void computeEnergyFromGrid(double& energy, const double box[3]) {
    // Log beginning of processing
    platform::log(LogLevel::DEBUG, "Computing energy from PME grid");
    
    // Use exactly the same calculation method as pme.cpp
    double volume = box[0] * box[1] * box[2];
    // Use exactly the same constants as pme.cpp
    double one_4pi_eps = 138.935456/pme_params.epsilon_r; // Ensure using the same Coulomb constant as pme.cpp
    double factor = M_PI*M_PI/(pme_params.alpha*pme_params.alpha);
    // Calculate boxfactor: exactly like pme.cpp
    double boxfactor = M_PI * volume;
    
    platform::log(LogLevel::DEBUG, "Computing energy from grid with box = [" + std::to_string(box[0]) + "," + 
                 std::to_string(box[1]) + "," + std::to_string(box[2]) + "], alpha = " + 
                 std::to_string(pme_params.alpha) + ", volume = " + std::to_string(volume));
    
    platform::log(LogLevel::DEBUG, "Energy parameters: one_4pi_eps = " + std::to_string(one_4pi_eps) + 
                 ", factor = " + std::to_string(factor) + ", boxfactor = " + 
                 std::to_string(boxfactor));
    
    // Get grid dimensions
    int nx = pme_params.meshSize[0];
    int ny = pme_params.meshSize[1];
    int nz = pme_params.meshSize[2];
    
    // Calculate reciprocal lattice vectors
    double recipBoxVectors[3][3] = {{0}};
    
    // Process box vectors correctly - maintain complete consistency with pme.cpp
    double periodicBoxVectors[3][3] = {
        {box[0], 0.0, 0.0},
        {0.0, box[1], 0.0},
        {0.0, 0.0, box[2]}
    };
    
    // Diagonal box check
    bool isDiagonalBox = true;  // Force diagonal box
    
    // Calculate reciprocal vectors exactly as in pme.cpp
    if (isDiagonalBox) {
        recipBoxVectors[0][0] = 1.0 / box[0]; 
        recipBoxVectors[1][1] = 1.0 / box[1]; 
        recipBoxVectors[2][2] = 1.0 / box[2]; 
    } else {
        // Non-diagonal box calculation
        double det = periodicBoxVectors[0][0] * (periodicBoxVectors[1][1] * periodicBoxVectors[2][2] - periodicBoxVectors[1][2] * periodicBoxVectors[2][1]) -
                     periodicBoxVectors[0][1] * (periodicBoxVectors[1][0] * periodicBoxVectors[2][2] - periodicBoxVectors[1][2] * periodicBoxVectors[2][0]) +
                     periodicBoxVectors[0][2] * (periodicBoxVectors[1][0] * periodicBoxVectors[2][1] - periodicBoxVectors[1][1] * periodicBoxVectors[2][0]);
        
        // Calculate reciprocal lattice vectors
        recipBoxVectors[0][0] = (periodicBoxVectors[1][1] * periodicBoxVectors[2][2] - periodicBoxVectors[1][2] * periodicBoxVectors[2][1]) / det;
        recipBoxVectors[0][1] = (periodicBoxVectors[0][2] * periodicBoxVectors[2][1] - periodicBoxVectors[0][1] * periodicBoxVectors[2][2]) / det;
        recipBoxVectors[0][2] = (periodicBoxVectors[0][1] * periodicBoxVectors[1][2] - periodicBoxVectors[0][2] * periodicBoxVectors[1][1]) / det;
        
        recipBoxVectors[1][0] = (periodicBoxVectors[1][2] * periodicBoxVectors[2][0] - periodicBoxVectors[1][0] * periodicBoxVectors[2][2]) / det;
        recipBoxVectors[1][1] = (periodicBoxVectors[0][0] * periodicBoxVectors[2][2] - periodicBoxVectors[0][2] * periodicBoxVectors[2][0]) / det;
        recipBoxVectors[1][2] = (periodicBoxVectors[0][2] * periodicBoxVectors[1][0] - periodicBoxVectors[0][0] * periodicBoxVectors[1][2]) / det;
        
        recipBoxVectors[2][0] = (periodicBoxVectors[1][0] * periodicBoxVectors[2][1] - periodicBoxVectors[1][1] * periodicBoxVectors[2][0]) / det;
        recipBoxVectors[2][1] = (periodicBoxVectors[0][1] * periodicBoxVectors[2][0] - periodicBoxVectors[0][0] * periodicBoxVectors[2][1]) / det;
        recipBoxVectors[2][2] = (periodicBoxVectors[0][0] * periodicBoxVectors[1][1] - periodicBoxVectors[0][1] * periodicBoxVectors[1][0]) / det;
    }
    
    // Initialize energy
    energy = 0.0;
    
    int maxkx = (nx+1)/2;
    int maxky = (ny+1)/2;
    int maxkz = (nz+1)/2;
    
    // Calculate energy exactly as in pme.cpp
    for (int kx = 0; kx < nx; kx++) {
        // Calculate frequency
        double mx = (kx < maxkx) ? kx : (kx-nx);
        double mhx = mx * recipBoxVectors[0][0];
        
        // Key modification: ensure completely consistent calculation with pme.cpp
        // Exactly reproduce the B-spline moduli application method from pme.cpp
        double bx = boxfactor * pme_params.bsplineModuli[0][kx];
        
        for (int ky = 0; ky < ny; ky++) {
            double my = (ky < maxky) ? ky : (ky-ny);
            // Consistent with pme.cpp, considering non-diagonal terms
            double mhy = mx*recipBoxVectors[1][0] + my*recipBoxVectors[1][1];
            
            // Note: Don't apply boxfactor to by, consistent with pme.cpp
            double by = pme_params.bsplineModuli[1][ky];
            
            for (int kz = 0; kz < nz; kz++) {
                // Skip zero frequency term - for neutral systems
                if (kx == 0 && ky == 0 && kz == 0) {
                    continue;
                }
                
                double mz = (kz < maxkz) ? kz : (kz-nz);
                // Consistent with pme.cpp, considering non-diagonal terms
                double mhz = mx*recipBoxVectors[2][0] + my*recipBoxVectors[2][1] + mz*recipBoxVectors[2][2];
                
                // Get grid data
                int index = kx * ny * nz + ky * nz + kz;
                double d1 = pme_params.pmeGrid[index].real();
                double d2 = pme_params.pmeGrid[index].imag();
                
                // Calculate convolution
                double m2 = mhx * mhx + mhy * mhy + mhz * mhz;
                
                // Don't apply boxfactor to bz, consistent with pme.cpp
                double bz = pme_params.bsplineModuli[2][kz];
                
                // Calculate denom exactly as in pme.cpp
                double denom = m2 * bx * by * bz;
                
                // Improve numerical stability
                if (denom < 1e-10) {
                    denom = 1e-10;
                }
                
                // Ensure energy calculation formula is exactly the same as in pme.cpp
                double eterm = one_4pi_eps * exp(-factor * m2) / denom;
                double struct2 = d1*d1 + d2*d2;
                double energyContrib = eterm * struct2;
                
                // Update grid value - exactly reproduce pme.cpp method
                std::complex<double> updatedValue(d1 * eterm, d2 * eterm);
                pme_params.pmeGrid[index] = updatedValue;
                
                // Accumulate energy
                energy += energyContrib;
            }
        }
    }
    
    // Consistent with pme.cpp: multiply by 0.5
    energy *= 0.5;
    
    // Output important statistics
    platform::log(LogLevel::INFO, "PME reciprocal energy = ", energy);
}

/**
 * @brief Compute reciprocal space energy using PME
 */
double computeReciprocalPME(model::MCState& state) {
    platform::log(LogLevel::INFO, "Computing PME reciprocal space energy");
    
    // Simplified system information output
    platform::log(LogLevel::INFO, "Box: [", state.info.box[0], ", ", 
                 state.info.box[1], ", ", state.info.box[2], "], Alpha: ", pme_params.alpha);
    
    // Simplified system charge check
    double totalCharge = 0.0;
    for (int i = 0; i < state.activeAtomCount; i++) {
        totalCharge += state.atoms[i].charge;
    }
    
    if (std::abs(totalCharge) > 1e-6) {
        platform::log(LogLevel::WARNING, "System is not neutral! Total charge = ", totalCharge);
    }
    
    // PME grid check
    if (pme_params.pmeGrid.empty()) {
        platform::log(LogLevel::ERROR, "PME grid not initialized!");
        return 0.0;
    }
    
    // Initialize B-spline functions
    if (pme_params.bsplineModuli[0].empty() || 
        pme_params.bsplineModuli[1].empty() || 
        pme_params.bsplineModuli[2].empty()) {
        platform::log(LogLevel::INFO, "Initializing B-splines for PME calculation...");
        pme_params.initializeBsplines();
    }
    
    // Reset grid
    std::fill(pme_params.pmeGrid.begin(), pme_params.pmeGrid.end(), std::complex<double>(0.0, 0.0));
    
    // Execute PME calculation steps
    spreadChargesOntoGrid(state, false);
    performFFTForward();
    
    // Convert box to double array
    double box[3];
    for (int i = 0; i < 3; i++) {
        box[i] = static_cast<double>(state.info.box[i]);
    }
    
    // Update box size in PME parameters, ensure B-splines and energy calculation use the same volume
    pme_params.setBox(box);
    
    // Calculate energy
    double energy = 0.0;
    computeEnergyFromGrid(energy, box);
    
    // Energy already includes Coulomb constant, no need to multiply again
    double reciprocal_energy = energy;
    
    platform::log(LogLevel::INFO, "Reciprocal space energy = ", reciprocal_energy);
    platform::log(LogLevel::DEBUG, "Reciprocal space energy = " + std::to_string(reciprocal_energy));
    
    return reciprocal_energy;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc 