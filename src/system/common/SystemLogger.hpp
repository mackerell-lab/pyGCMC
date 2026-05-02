#pragma once

#include "LoggingState.hpp"

#include <iostream>
#include <sstream>
#include <string>

namespace pygcmc {
namespace system {
namespace common {

/**
 * @brief System-wide logger for PyGCMC
 *
 * This logger provides basic logging functionality with different log levels.
 * It was extracted from simulation.hpp to be reusable across modules.
 */
class SystemLogger {
public:
    // Configuration methods
    static void setVerbose(bool verbose) { LoggingState::system_enabled = verbose; }
    static void setLogLevel(LogLevel level) { LoggingState::system_level = level; }
    static void setDebugMode(bool debug) { LoggingState::system_debug_mode = debug; }

    static bool isVerbose() { return LoggingState::system_enabled; }
    static bool isDebugEnabled() {
        return LoggingState::system_enabled && LoggingState::system_level <= LogLevel::DEBUG;
    }
    static bool isDebugMode() { return LoggingState::system_debug_mode; }

    // Logging method with variadic templates
    template<typename... Args>
    static void log(LogLevel level, Args... args) {
        if (!LoggingState::should_log_system(level)) return;

        std::stringstream ss;
        (ss << ... << args);

        switch (level) {
            case LogLevel::DEBUG:
                std::cout << "[DEBUG] ";
                break;
            case LogLevel::INFO:
                std::cout << "[INFO] ";
                break;
            case LogLevel::WARNING:
                std::cout << "[WARN] ";
                break;
            case LogLevel::ERROR:
                std::cout << "[ERROR] ";
                break;
        }
        std::cout << ss.str() << std::endl;
    }

    // Convenience methods
    template<typename... Args>
    static void debug(Args... args) { log(LogLevel::DEBUG, args...); }

    template<typename... Args>
    static void info(Args... args) { log(LogLevel::INFO, args...); }

    template<typename... Args>
    static void warn(Args... args) { log(LogLevel::WARNING, args...); }

    template<typename... Args>
    static void error(Args... args) { log(LogLevel::ERROR, args...); }
};

} // namespace common
} // namespace system
} // namespace pygcmc
