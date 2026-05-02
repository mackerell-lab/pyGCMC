#include "MemoryGlobalCleanup.hpp"
#include "../pme/PMEGlobal.hpp"
#include "../pme/PMECore.hpp"
#include "../pme/PMEFFTCore.hpp"
#include "../pgp/PGPGlobal.hpp"
#include "../pgp/PGPCore.hpp"
#include "MemoryPool.hpp"
#include "MemorySafetyChecks.hpp"
#include "platform/platform.hpp"
#include <thread>
#include <chrono>

namespace pygcmc {
namespace platform {
namespace cpu {

void cleanupAllGlobalState() {
    platform::log(LogLevel::INFO, "Starting global state cleanup...");

    // IMPORTANT: Avoid double-free by not calling functions that reset the same pointers

    // Only clear FFT weights - this is safe as it's just a vector clear
    try {
        CustomFFT::clearFFTWeights();
    } catch (...) {
        // Ignore errors during cleanup
    }

    // Reset smart pointers ONCE
    try {
        resetPGPParamsPtr();
    } catch (...) {
        // Ignore errors
    }

    try {
        resetPMEParamsPtr();
    } catch (...) {
        // Ignore errors
    }

    // Clear memory pool to release all cached buffers
    try {
        getGridMemoryPool().clear();
        platform::log(LogLevel::DEBUG, "Memory pool cleared.");
    } catch (...) {
        // Ignore errors
    }

    // Mark that cleanup has been called
    SafetyChecks::markCleanupCalled();

    platform::log(LogLevel::INFO, "Global state cleanup completed.");
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
