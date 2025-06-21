#pragma once
#ifndef PYGCMC_SYSTEM_LOG_LOGMAIN_HPP
#define PYGCMC_SYSTEM_LOG_LOGMAIN_HPP

#include "../common/SystemInterface.hpp"
#include <sstream>
#include <iostream>

namespace pygcmc {
namespace system {
namespace log {

/**
 * @brief Main logging implementation
 */
class LogMain : public common::ILogger {
public:
    // ILogger interface implementation
    void setVerbose(bool enable) override;
    void setLogLevel(common::LogLevel level) override;
    common::LogLevel getLogLevel() const override;
    bool isVerbose() const override;
    
    // Static interface for backward compatibility
    static void set_verbose(bool verbose);
    static void set_log_level(common::LogLevel level);
    static bool get_verbose();
    static common::LogLevel get_log_level();
    
    // Template logging function (must be in header for template instantiation)
    template<typename... Args>
    static void log(common::LogLevel level, Args... args) {
        if (!verbose_ || level < log_level_) return;
        
        std::stringstream ss;
        (ss << ... << args);
        
        switch (level) {
            case common::LogLevel::DEBUG:
                std::cout << "[DEBUG] ";
                break;
            case common::LogLevel::INFO:
                std::cout << "[INFO] ";
                break;
            case common::LogLevel::WARNING:
                std::cout << "[WARNING] ";
                break;
            case common::LogLevel::ERROR:
                std::cout << "[ERROR] ";
                break;
        }
        std::cout << ss.str() << std::endl;
    }
    
    // Convenience logging functions
    template<typename... Args>
    static void debug(Args... args) {
        log(common::LogLevel::DEBUG, args...);
    }
    
    template<typename... Args>
    static void info(Args... args) {
        log(common::LogLevel::INFO, args...);
    }
    
    template<typename... Args>
    static void warning(Args... args) {
        log(common::LogLevel::WARNING, args...);
    }
    
    template<typename... Args>
    static void error(Args... args) {
        log(common::LogLevel::ERROR, args...);
    }

private:
    // Static members for global logging state
    static bool verbose_;
    static common::LogLevel log_level_;
    
    // Instance members for object-based logging
    bool instance_verbose_;
    common::LogLevel instance_log_level_;
};

/**
 * @brief Initialize logging system with specified settings
 * 
 * @param verbose Enable verbose logging output
 * @param level Minimum logging level to display
 */
inline void initializeLogging(bool verbose = false, common::LogLevel level = common::LogLevel::INFO) {
    LogMain::set_verbose(verbose);
    LogMain::set_log_level(level);
}

/**
 * @brief Quick access to set verbose mode
 * 
 * @param verbose Enable verbose mode
 */
inline void setVerbose(bool verbose) {
    LogMain::set_verbose(verbose);
}

/**
 * @brief Quick access to set log level
 * 
 * @param level Log level to set
 */
inline void setLogLevel(common::LogLevel level) {
    LogMain::set_log_level(level);
}

} // namespace log
} // namespace system
} // namespace pygcmc

#endif // PYGCMC_SYSTEM_LOG_LOGMAIN_HPP 