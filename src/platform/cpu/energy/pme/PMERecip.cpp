#include "PMERecip.hpp"
#include "PMECore.hpp"
#include "PMEFFT.hpp"
#include "PMESystemCore.hpp"
#include "PMEGridCharge.hpp"
#include "platform/platform.hpp"
#include <cmath>
#include <algorithm>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Perform forward FFT on the grid
 * 
 * Uses custom FFT implementation
 */
void performFFTForward() {
    // Console output - same as pme.cpp
    platform::log(LogLevel::DEBUG, "Performing forward FFT on PME grid");
    
    // For debugging backup data, only execute in debug mode
    if (platform::is_debug_mode()) {
        // In debug mode, backup grid data for comparison
        extern std::vector<std::complex<double>> fftGridBackup;
        fftGridBackup.resize(pme_params.pmeGrid.size());
        std::copy(pme_params.pmeGrid.begin(), pme_params.pmeGrid.end(), fftGridBackup.begin());
    }
    
    // Only execute test code in debug mode
    if (platform::is_debug_mode()) {
        // Calculate non-zero points before FFT
        int nonZeroBeforeFFT = 0;
        for (size_t i = 0; i < pme_params.pmeGrid.size(); i++) {
            if (std::abs(pme_params.pmeGrid[i]) > 1e-10) nonZeroBeforeFFT++;
        }
        platform::log(LogLevel::DEBUG, "Grid before FFT: non-zero points = " + std::to_string(nonZeroBeforeFFT));
    }
    
    // Get grid dimensions
    int nx = pme_params.meshSize[0];
    int ny = pme_params.meshSize[1];
    int nz = pme_params.meshSize[2];
    
    // Use custom FFT implementation
    CustomFFT::fft3D_forward(pme_params.pmeGrid.data(), nx, ny, nz);
    
    // Non-zero point counting - only perform detailed counting in debug mode
    if (platform::is_debug_mode()) {
        // Calculate non-zero points after FFT
        int nonZeroAfterFFT = 0;
        for (size_t i = 0; i < pme_params.pmeGrid.size(); i++) {
            if (std::abs(pme_params.pmeGrid[i]) > 1e-10) nonZeroAfterFFT++;
        }
        platform::log(LogLevel::INFO, "Grid after FFT: non-zero points = ", nonZeroAfterFFT);
        platform::log(LogLevel::DEBUG, "Grid after FFT: non-zero points = " + std::to_string(nonZeroAfterFFT));
    }
}

/**
 * @brief Perform backward FFT on the grid
 * 
 * Uses custom FFT implementation, matches pme.cpp implementation
 */
void performFFTBackward() {
    // Get grid dimensions
    int nx = pme_params.meshSize[0];
    int ny = pme_params.meshSize[1];
    int nz = pme_params.meshSize[2];
    
    // Only execute test code in debug mode
    if (platform::is_debug_mode()) {
        // Calculate non-zero points before backward FFT
        int nonZeroBeforeFFT = 0;
        for (size_t i = 0; i < pme_params.pmeGrid.size(); i++) {
            if (std::abs(pme_params.pmeGrid[i]) > 1e-10) nonZeroBeforeFFT++;
        }
        platform::log(LogLevel::DEBUG, "Grid before backward FFT: non-zero points = " + std::to_string(nonZeroBeforeFFT));
    }
    
    // Use custom FFT implementation
    CustomFFT::fft3D_backward(pme_params.pmeGrid.data(), nx, ny, nz);
    
    // Only execute test code in debug mode
    if (platform::is_debug_mode()) {
        // Calculate non-zero points after backward FFT
        int nonZeroAfterFFT = 0;
        for (size_t i = 0; i < pme_params.pmeGrid.size(); i++) {
            if (std::abs(pme_params.pmeGrid[i]) > 1e-10) nonZeroAfterFFT++;
        }
        platform::log(LogLevel::DEBUG, "Grid after backward FFT: non-zero points = " + std::to_string(nonZeroAfterFFT));
    }
}

// Functions moved to PMEEnergyCalc.cpp

// <agent-hook:reciprocal_implementation>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 