#include "PGPCore.hpp"
#include "PGPGlobal.hpp"
#include "../../../platform.hpp"
#include <thread>
#include <chrono>

namespace pygcmc {
namespace platform {
namespace cpu {

// Note: pgp_mutex is now pgp_global_mutex defined in PGPGlobal.cpp

void resetPGPState() {
    platform::log(LogLevel::INFO, "Resetting PGP state completely");

    // Simply reset the smart pointer - this will automatically:
    // 1. Call destructor of PGPParams which cleans up all vectors
    // 2. Free all memory
    // 3. Next call to getPGPParams() will create a fresh instance
    resetPGPParamsPtr();

    platform::log(LogLevel::INFO, "PGP state reset complete");
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
