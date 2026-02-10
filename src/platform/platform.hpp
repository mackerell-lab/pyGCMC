#pragma once

#include "../model/ModelModule.hpp"
#include "../system/common/LoggingState.hpp"
#include <string>
#include <sstream>
#include <iostream>
#include <vector>

namespace pygcmc {
namespace platform {

// Log level enum
enum class LogLevel {
    DEBUG,
    INFO,
    WARNING,
    ERROR
};

// Logging functions
inline void set_verbose(bool verbose) { ::pygcmc::system::common::LoggingState::platform_enabled = verbose; }
inline void set_log_level(LogLevel level) {
    ::pygcmc::system::common::LoggingState::platform_level_int = static_cast<int>(level);
}
inline void set_debug_mode(bool debug_mode) { ::pygcmc::system::common::LoggingState::platform_debug_mode = debug_mode; }

// Helper to check if debug mode is enabled (for test output)
inline bool is_debug_mode() { return ::pygcmc::system::common::LoggingState::platform_debug_mode; }

template<typename... Args>
inline void log(LogLevel level, Args... args) {
    if (!::pygcmc::system::common::LoggingState::should_log_platform_int(static_cast<int>(level))) return;
    
    std::stringstream ss;
    (ss << ... << args);
    
    switch (level) {
        case LogLevel::DEBUG:
            std::cout << "[PLATFORM DEBUG] ";
            break;
        case LogLevel::INFO:
            std::cout << "[PLATFORM INFO] ";
            break;
        case LogLevel::WARNING:
            std::cout << "[PLATFORM WARNING] ";
            break;
        case LogLevel::ERROR:
            std::cout << "[PLATFORM ERROR] ";
            break;
    }
    std::cout << ss.str() << std::endl;
}

class IPlatform {
public:
    virtual ~IPlatform() = default;

    // Initializes the platform with the initial state (e.g. uploading data to GPU or storing for CPU computation)
    virtual void initialize(const model::MCState& state) = 0;

    // Finalizes the platform by cleaning up resources or optionally downloading final data to CPU
    virtual void finalize() = 0;

    // Computes the total non-bonded energy of the current system state
    virtual float computeTotalEnergy() = 0;

    // Computes the energy contribution for a specific residue
    virtual float computeResidueEnergy(int residueIndex) = 0;

    // Attempts a Monte Carlo move (e.g., translate, insert, delete) and updates the state accordingly
    virtual bool attemptMove(const std::string& moveType) = 0;
};

} // namespace platform
} // namespace pygcmc 
