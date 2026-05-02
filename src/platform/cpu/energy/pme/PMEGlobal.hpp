#pragma once

#include <memory>
#include <mutex>
#include "PMECore.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

// Global PME parameters management
extern std::unique_ptr<PMEParams> pme_params_ptr;
extern std::mutex pme_global_mutex;

// Get PME parameters (thread-safe lazy initialization)
inline PMEParams& getPMEParams() {
    if (!pme_params_ptr) {
        std::lock_guard<std::mutex> lock(pme_global_mutex);
        if (!pme_params_ptr) {
            pme_params_ptr = std::make_unique<PMEParams>();
        }
    }
    return *pme_params_ptr;
}

// Reset PME parameters (for testing)
inline void resetPMEParamsPtr() {
    std::lock_guard<std::mutex> lock(pme_global_mutex);
    pme_params_ptr.reset();
}

// Compatibility macro for easier migration
#define pme_params (getPMEParams())

} // namespace cpu
} // namespace platform
} // namespace pygcmc
