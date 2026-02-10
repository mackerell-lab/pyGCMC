#pragma once

#include "SystemInterface.hpp"

namespace pygcmc {
namespace system {
namespace common {

/**
 * @brief Shared logging configuration (SYSTEM + PLATFORM components)
 *
 * This header centralizes logging state across translation units using C++17
 * inline variables. It keeps the hot path fast: emitters should always
 * "check first, format second".
 */
struct LoggingState {
    // SYSTEM component (SystemLogger + LogMain)
    inline static bool system_enabled = false;
    inline static LogLevel system_level = LogLevel::INFO;
    inline static bool system_debug_mode = false;

    // PLATFORM component (platform::log / platform::is_debug_mode)
    inline static bool platform_enabled = false;
    // Keep in sync with pygcmc::platform::LogLevel ordering:
    // DEBUG=0, INFO=1, WARNING=2, ERROR=3.
    inline static int platform_level_int = 1; // INFO
    inline static bool platform_debug_mode = false;

    static inline bool should_log_system(LogLevel level) noexcept {
        return system_enabled && level >= system_level;
    }

    static inline bool should_log_platform_int(int level_int) noexcept {
        return platform_enabled && level_int >= platform_level_int;
    }
};

} // namespace common
} // namespace system
} // namespace pygcmc

