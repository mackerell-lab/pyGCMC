#pragma once

// Common components
#include "common/SystemInterface.hpp"
#include "common/SystemConstants.hpp"
#include "common/SystemUtils.hpp"

// Logging system
#include "log/LogMain.hpp"

// Molecular system
#include "molecular/MolecularComposite.hpp"

// Monte Carlo system
#include "montecarlo/MCComposite.hpp"

// Backward compatibility layers
#include "common/SystemMain.hpp"
#include "molecular/MolecularMain.hpp"
#include "montecarlo/MCMain.hpp"

/**
 * @file SystemModule.hpp
 * @brief Unified entry point for the system module
 * 
 * This file provides a single include point for all system functionality
 * including common utilities, logging, molecular system management,
 * and Monte Carlo simulations.
 * 
 * The module is organized into the following components:
 * 
 * - common/: Shared utilities, constants, and interfaces
 * - log/: Logging system with configurable levels
 * - molecular/: Molecular system building and management
 * - montecarlo/: Monte Carlo simulation functionality
 * 
 * Backward compatibility is maintained through compatibility layers
 * that preserve the original API while using the new modular architecture.
 */

namespace pygcmc {
namespace system {

/**
 * @brief Factory function for creating system instances
 * 
 * @param kind Type of system to create
 * @return Unique pointer to the created system
 */
std::unique_ptr<common::ISystem> createSystem(common::SystemKind kind);

/**
 * @brief Initialize logging system
 * 
 * @param verbose Enable verbose logging
 * @param level Logging level to use
 */
void initializeLogging(bool verbose = false, common::LogLevel level = common::LogLevel::INFO);

/**
 * @brief Get version information for the system module
 * 
 * @return Version string
 */
std::string getSystemModuleVersion();

} // namespace system
} // namespace pygcmc 