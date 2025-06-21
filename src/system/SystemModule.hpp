#pragma once

#include <memory>
#include <string>
#include <stdexcept>

// Common components
#include "common/SystemInterface.hpp"
#include "common/SystemConstants.hpp"
#include "common/SystemUtils.hpp"
#include "common/SystemMain.hpp"

// Logging system
#include "log/LogMain.hpp"

// Molecular system
#include "molecular/MolecularComposite.hpp"
#include "molecular/MolecularMain.hpp"

// Monte Carlo system
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

// Re-export all major classes for easy access
// Note: System class is already defined in this namespace in common/SystemMain.hpp
using MolecularSystem = molecular::MolecularMain;
using MonteCarloSystem = montecarlo::MCMain;

// Re-export common types
using SystemKind = common::SystemKind;
using LogLevel = common::LogLevel;
using MovementMolecularInfo = montecarlo::MovementMolecularInfo;

/**
 * @brief Factory function for creating system instances
 * 
 * @param kind Type of system to create
 * @return Unique pointer to the created system
 * @throws std::runtime_error if system type is not supported
 */
inline std::unique_ptr<common::ISystem> createSystem(SystemKind kind) {
    switch (kind) {
        case SystemKind::MOLECULAR:
            // For now, return nullptr as molecular system doesn't implement ISystem yet
            // This can be extended when needed
            throw std::runtime_error("Molecular system factory not yet implemented");
            
        case SystemKind::MONTE_CARLO:
            // For now, return nullptr as Monte Carlo system doesn't implement ISystem yet  
            // This can be extended when needed
            throw std::runtime_error("Monte Carlo system factory not yet implemented");
            
        default:
            throw std::runtime_error("Unknown system kind");
    }
}

/**
 * @brief Initialize logging system with specified settings
 * 
 * @param verbose Enable verbose logging output
 * @param level Minimum logging level to display
 */
inline void initializeLogging(bool verbose = false, LogLevel level = LogLevel::INFO) {
    log::LogMain::set_verbose(verbose);
    log::LogMain::set_log_level(level);
}

/**
 * @brief Get version information for the system module
 * 
 * @return Version string indicating refactored modular design
 */
inline std::string getSystemModuleVersion() {
    return "1.0.0-refactored";
}

/**
 * @brief Quick access to logging functionality
 * 
 * @param verbose Enable verbose mode
 */
inline void setVerbose(bool verbose) {
    log::LogMain::set_verbose(verbose);
}

/**
 * @brief Quick access to log level setting
 * 
 * @param level Log level to set
 */
inline void setLogLevel(LogLevel level) {
    log::LogMain::set_log_level(level);
}

} // namespace system
} // namespace pygcmc 