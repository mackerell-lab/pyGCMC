#pragma once

// Include all system main modules
#include "common/SystemMain.hpp"
#include "log/LogMain.hpp"
#include "molecular/MolecularMain.hpp"
#include "montecarlo/MCMain.hpp"

namespace pygcmc {
namespace system {

/**
 * @brief System Module - Unified Entry Point
 * 
 * This module provides comprehensive system management functionality including:
 * 1. System initialization and parameter management
 * 2. Logging system with configurable levels
 * 3. Molecular system building from PDB/PSF/topology files
 * 4. Monte Carlo GCMC simulation engine with energy calculations
 * 
 * Usage examples:
 * 
 * // Initialize logging
 * initializeLogging(true, LogLevel::DEBUG);
 * 
 * // Use molecular system
 * MolecularSystem molSys;
 * auto molecular = molSys.combine(structure, topology);
 * 
 * // Use Monte Carlo system
 * MonteCarloSystem mcSys;
 * mcSys.initializeFromMolecular(molecular);
 */

// Re-export main classes for backward compatibility
using MolecularSystem = molecular::MolecularMain;
using MonteCarloSystem = montecarlo::MCMain;

// Re-export common types
using LogLevel = common::LogLevel;
using MovementMolecularInfo = montecarlo::MovementMolecularInfo;

/**
 * @brief Initialize logging system
 * 
 * @param verbose Enable verbose logging output
 * @param level Minimum logging level to display
 */
inline void initializeLogging(bool verbose = false, LogLevel level = LogLevel::INFO) {
    log::LogMain::set_verbose(verbose);
    log::LogMain::set_log_level(level);
}

/**
 * @brief Set verbose logging mode
 * 
 * @param verbose Enable verbose mode
 */
inline void setVerbose(bool verbose) {
    log::LogMain::set_verbose(verbose);
}

/**
 * @brief Set minimum log level
 * 
 * @param level Log level to set
 */
inline void setLogLevel(LogLevel level) {
    log::LogMain::set_log_level(level);
}

/**
 * @brief Get system module version
 * 
 * @return Version string
 */
inline std::string getSystemVersion() {
    return "1.0.0-refactored";
}

} // namespace system
} // namespace pygcmc 