#pragma once
#ifndef PYGCMC_SYSTEM_COMMON_SYSTEMMAIN_HPP
#define PYGCMC_SYSTEM_COMMON_SYSTEMMAIN_HPP

#include "SystemInterface.hpp"
#include "../log/LogMain.hpp"
#include "model/param.hpp"

namespace pygcmc {
namespace system {

/**
 * @brief Original System class for backward compatibility
 * 
 * This class maintains the exact same interface as the original System class,
 * but internally delegates to the new modular implementation.
 */
class System {
public:
    // Static log control (backward compatibility)
    static bool verbose_;
    static common::LogLevel log_level_;
    
    // Logging functionality (delegated to log::LogMain)
    static void set_verbose(bool verbose) {
        log::LogMain::set_verbose(verbose);
        verbose_ = verbose;  // Keep static member in sync
    }
    
    static void set_log_level(common::LogLevel level) {
        log::LogMain::set_log_level(level);
        log_level_ = level;  // Keep static member in sync
    }
    
    // Template logging function (delegated to log::LogMain)
    template<typename... Args>
    static void log(common::LogLevel level, Args... args) {
        log::LogMain::log(level, args...);
    }

    // Original functionality from system.cpp
    void initialize_parameters();
    void process_cavity_list();
    void initialize_mc_time_list();

private:
    model::Param params_;  // System parameters
};

// Static member initialization for backward compatibility
inline bool System::verbose_ = false;
inline common::LogLevel System::log_level_ = common::LogLevel::INFO;

} // namespace system
} // namespace pygcmc

#endif // PYGCMC_SYSTEM_COMMON_SYSTEMMAIN_HPP 