#pragma once

// Include all system sub-modules
#include "common/SystemInterface.hpp"
#include "common/SystemConstants.hpp"
#include "common/SystemUtils.hpp"
#include "common/SystemMain.hpp"
#include "common/SystemFactory.hpp"

#include "log/LogMain.hpp"

#include "molecular/MolecularComposite.hpp"
#include "molecular/MolecularMain.hpp"

#include "montecarlo/MCComposite.hpp"
#include "montecarlo/MCMain.hpp"

/**
 * @file SystemModule.hpp
 * @brief Unified System Module - Single Entry Point
 * 
 * This is the ONLY header file you need to include to access all system functionality.
 * It provides comprehensive system management capabilities including:
 * 
 * 1. **Common utilities**: Constants, utility functions, and interfaces
 * 2. **Logging system**: Configurable logging with multiple levels
 * 3. **Molecular system**: Building and managing molecular structures from PDB/PSF/topology files
 * 4. **Monte Carlo system**: GCMC simulation engine with energy calculations
 * 
 * **Module Organization:**
 * ```
 * system/
 * ├── common/         # Shared utilities, constants, and interfaces
 * ├── log/            # Logging system with configurable levels  
 * ├── molecular/      # Molecular system building and management
 * └── montecarlo/     # Monte Carlo simulation functionality
 * ```
 * 
 * **Usage Examples:**
 * 
 * ```cpp
 * // Include this single header for all system functionality
 * #include "system/SystemModule.hpp"
 * using namespace pygcmc::system;
 * 
 * // Initialize logging
 * initializeLogging(true, common::LogLevel::DEBUG);
 * 
 * // Use molecular system (backward compatible API)
 * MolecularSystem molSys;
 * auto molecular = molSys.combine(structure, topology);
 * 
 * // Use Monte Carlo system (backward compatible API)
 * MonteCarloSystem mcSys;
 * mcSys.initializeFromMolecular(molecular);
 * 
 * // Use modern modular API
 * auto molComposite = std::make_unique<molecular::MolecularComposite>();
 * auto mcComposite = std::make_unique<montecarlo::MCComposite>();
 * ```
 * 
 * **Backward Compatibility:**
 * All original APIs (System, MolecularSystem, MonteCarloSystem) are preserved
 * and work exactly as before. The implementation is now modular and maintainable.
 */

namespace pygcmc {
namespace system {

// Re-export major classes for easy access
// Note: System class is already defined in this namespace in common/SystemMain.hpp
using MolecularSystem = molecular::MolecularMain;
using MonteCarloSystem = montecarlo::MCMain;

// Re-export common types
using SystemKind = common::SystemKind;
using LogLevel = common::LogLevel;
using MovementMolecularInfo = montecarlo::MovementMolecularInfo;

/**
 * @brief Factory function for creating system instances
 * Delegates to common::createSystem()
 */
inline std::unique_ptr<common::ISystem> createSystem(SystemKind kind) {
    return common::createSystem(kind);
}

/**
 * @brief Initialize logging system with specified settings
 * Delegates to log::initializeLogging()
 */
inline void initializeLogging(bool verbose = false, LogLevel level = LogLevel::INFO) {
    return log::initializeLogging(verbose, level);
}

/**
 * @brief Get version information for the system module
 * Delegates to common::getSystemModuleVersion()
 */
inline std::string getSystemModuleVersion() {
    return common::getSystemModuleVersion();
}

/**
 * @brief Quick access to set verbose logging
 * Delegates to log::setVerbose()
 */
inline void setVerbose(bool verbose) {
    return log::setVerbose(verbose);
}

/**
 * @brief Quick access to set log level
 * Delegates to log::setLogLevel()
 */
inline void setLogLevel(LogLevel level) {
    return log::setLogLevel(level);
}

} // namespace system
} // namespace pygcmc 