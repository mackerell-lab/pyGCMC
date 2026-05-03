#include "MemoryGlobalCleanup.hpp"
#include "../pme/PMEGlobal.hpp"
#include "../pme/PMECore.hpp"
#include "../pme/PMEFFTCore.hpp"
#include "../pgp/PGPGlobal.hpp"
#include "../pgp/PGPCore.hpp"
#include "MemoryPool.hpp"
#include "platform/Platform.hpp"
#include <thread>
#include <chrono>
#include <new>

// For malloc_trim on GNU/Linux systems
#ifdef __linux__
#include <malloc.h>
#endif

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Force complete reset of all global state
 *
 * This is a more aggressive cleanup that completely resets all global state
 * to prevent any memory corruption issues in extreme cases.
 */
void forceCompleteReset() {
    platform::log(LogLevel::INFO, "Starting forced complete reset...");

    // Step 1: Clear FFT weights
    try {
        CustomFFT::clearFFTWeights();
    } catch (...) {
        // Ignore errors
    }

    // Step 2: Clear memory pool
    try {
        getGridMemoryPool().clear();
    } catch (...) {
        // Ignore errors
    }

    // Step 3: Force reset of global parameters
    // This is more aggressive than just reset - we destroy and recreate
    try {
        // Lock both mutexes
        std::lock_guard<std::mutex> pme_lock(pme_global_mutex);
        std::lock_guard<std::mutex> pgp_lock(pgp_global_mutex);

        // Destroy existing instances
        if (pme_params_ptr) {
            pme_params_ptr.reset();
            // Force deallocation
            pme_params_ptr = nullptr;
        }

        if (pgp_params_ptr) {
            pgp_params_ptr.reset();
            // Force deallocation
            pgp_params_ptr = nullptr;
        }

        // Give OS time to reclaim memory
        std::this_thread::sleep_for(std::chrono::milliseconds(10));

        // Create new instances
        pme_params_ptr = std::make_unique<PMEParams>();
        pgp_params_ptr = std::make_unique<PGPParams>();

    } catch (...) {
        platform::log(LogLevel::ERROR, "Error during forced reset");
    }

    // Step 4: Force garbage collection hint
    // This is a hint to the memory allocator to release unused memory
#ifdef __linux__
    // malloc_trim is a GNU extension, only available on Linux
    malloc_trim(0);
#endif

    platform::log(LogLevel::INFO, "Forced complete reset completed.");
}

/**
 * @brief Safe cleanup that prevents memory corruption
 *
 * This version is safer for extreme cases where cleanup is called
 * multiple times or after many iterations.
 */
void safeCleanupAllGlobalState() {
    static std::mutex cleanup_mutex;
    static bool cleanup_in_progress = false;

    // Prevent concurrent cleanup
    std::lock_guard<std::mutex> lock(cleanup_mutex);

    if (cleanup_in_progress) {
        platform::log(LogLevel::WARNING, "Cleanup already in progress, skipping...");
        return;
    }

    cleanup_in_progress = true;

    try {
        // Use the more aggressive reset
        forceCompleteReset();
    } catch (...) {
        platform::log(LogLevel::ERROR, "Exception during cleanup, continuing anyway");
    }

    cleanup_in_progress = false;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
