#pragma once

#include "PGPCore.hpp"
#include <memory>
#include <mutex>

namespace pygcmc {
namespace platform {
namespace cpu {

// Use smart pointer to manage global PGP parameters
extern std::unique_ptr<PGPParams> pgp_params_ptr;
extern std::mutex pgp_global_mutex;

// Thread-safe function to get PGP parameters
inline PGPParams& getPGPParams() {
    if (!pgp_params_ptr) {
        std::lock_guard<std::mutex> lock(pgp_global_mutex);
        if (!pgp_params_ptr) {
            pgp_params_ptr = std::make_unique<PGPParams>();
        }
    }
    return *pgp_params_ptr;
}

// Reset PGP parameters
inline void resetPGPParamsPtr() {
    std::lock_guard<std::mutex> lock(pgp_global_mutex);
    pgp_params_ptr.reset();
    // New instance will be created on next call to getPGPParams()
}

// Not using macro to avoid recursive call issues

} // namespace cpu
} // namespace platform
} // namespace pygcmc