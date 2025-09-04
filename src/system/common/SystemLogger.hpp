#pragma once

#include <iostream>
#include <sstream>
#include <string>

namespace pygcmc {
namespace system {
namespace common {

/**
 * @brief Log level enumeration
 */
enum class LogLevel {
    DEBUG,
    INFO,
    WARNING,
    ERROR
};

/**
 * @brief System-wide logger for PyGCMC
 * 
 * This logger provides basic logging functionality with different log levels.
 * It was extracted from simulation.hpp to be reusable across modules.
 */
class SystemLogger {
private:
    static bool verbose_;
    static LogLevel log_level_;
    static bool debug_mode_;

public:
    // Configuration methods
    static void setVerbose(bool verbose) { verbose_ = verbose; }
    static void setLogLevel(LogLevel level) { log_level_ = level; }
    static void setDebugMode(bool debug) { debug_mode_ = debug; }
    
    static bool isVerbose() { return verbose_; }
    static bool isDebugEnabled() { return verbose_ && log_level_ <= LogLevel::DEBUG; }
    static bool isDebugMode() { return debug_mode_; }
    
    // Logging method with variadic templates
    template<typename... Args>
    static void log(LogLevel level, Args... args) {
        if (!verbose_ || level < log_level_) return;
        
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

// Static member definitions
inline bool SystemLogger::verbose_ = false;
inline LogLevel SystemLogger::log_level_ = LogLevel::WARNING;
inline bool SystemLogger::debug_mode_ = false;

} // namespace common
} // namespace system
} // namespace pygcmc