#include "PMECore.hpp"
#include "PMEConfig.hpp"
#include "PMEBSpline.hpp"
#include "PMEFFTCore.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

// Global PME parameters instance
PMEParams pme_params;

// Implementation moved to PMEParams.cpp and PMEBSpline.cpp
// This file now only contains the global PME parameters instance

void clearPMEState() {
    // Reset to default-constructed state
    pme_params = PMEParams{};
    
    // Clear FFT static weights to force regeneration on next use
    CustomFFT::clearFFTWeights();
}

// <agent-hook:pme_core_impl>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 