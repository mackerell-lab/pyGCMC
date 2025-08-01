#include "PMECore.hpp"
#include "PMEGlobal.hpp"
#include "PMEConfig.hpp"
#include "PMEBSpline.hpp"
#include "PMEFFTCore.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

// Global PME parameters instance - moved to PMEGlobal.cpp
// PMEParams pme_params; // REMOVED - now using smart pointer

// Implementation moved to PMEParams.cpp and PMEBSpline.cpp
// This file now only contains the global PME parameters instance

void clearPMEState() {
    // Reset PME parameters using smart pointer
    resetPMEParamsPtr();
    
    // Clear FFT static weights to force regeneration on next use
    CustomFFT::clearFFTWeights();
}

// <agent-hook:pme_core_impl>

} // namespace cpu
} // namespace platform
} // namespace pygcmc 