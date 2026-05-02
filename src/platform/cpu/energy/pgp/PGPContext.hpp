#pragma once

#include <memory>
#include <vector>
#include <complex>
#include <mutex>
#include "PGPCore.hpp"
#include "../pme/PMECore.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Context manager for PGP calculations
 *
 * This class provides isolated context for PGP calculations to avoid
 * global state corruption during concurrent or repeated operations
 */
class PGPContext {
private:
    // Local copies of parameters
    PMEParams local_pme_params;
    PGPParams local_pgp_params;

    // Backup of global state
    PMEParams* global_pme_backup;
    PGPParams* global_pgp_backup;

    // Thread safety
    static std::mutex context_mutex;
    bool active = false;

public:
    PGPContext() : global_pme_backup(nullptr), global_pgp_backup(nullptr) {}

    /**
     * @brief Enter the context - save global state and setup local
     */
    void enter() {
        std::lock_guard<std::mutex> lock(context_mutex);

        if (active) {
            throw std::runtime_error("PGPContext already active");
        }

        // Create deep copies of global parameters
        local_pme_params = getPMEParams();
        local_pgp_params = getPGPParams();

        // Save pointers to globals (we'll swap them)
        global_pme_backup = &getPMEParams();
        global_pgp_backup = &getPGPParams();

        active = true;
    }

    /**
     * @brief Exit the context - restore global state
     */
    void exit() {
        if (!active) return;

        std::lock_guard<std::mutex> lock(context_mutex);

        // No need to restore anything since we used local copies
        active = false;
    }

    /**
     * @brief Get local PME parameters
     */
    PMEParams& getLocalPMEParams() {
        if (!active) {
            throw std::runtime_error("PGPContext not active");
        }
        return local_pme_params;
    }

    /**
     * @brief Get local PGP parameters
     */
    PGPParams& getLocalPGPParams() {
        if (!active) {
            throw std::runtime_error("PGPContext not active");
        }
        return local_pgp_params;
    }

    ~PGPContext() {
        if (active) {
            exit();
        }
    }
};

// Initialize static member
inline std::mutex PGPContext::context_mutex;

/**
 * @brief RAII helper for PGPContext
 */
class PGPContextGuard {
private:
    PGPContext& ctx;

public:
    explicit PGPContextGuard(PGPContext& context) : ctx(context) {
        ctx.enter();
    }

    ~PGPContextGuard() {
        ctx.exit();
    }

    // Delete copy operations
    PGPContextGuard(const PGPContextGuard&) = delete;
    PGPContextGuard& operator=(const PGPContextGuard&) = delete;
};

} // namespace cpu
} // namespace platform
} // namespace pygcmc
