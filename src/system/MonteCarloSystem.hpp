/**
 * @file MonteCarloSystem.hpp
 * @brief Legacy compatibility header for MonteCarloSystem
 * 
 * This file provides backward compatibility by including the new modular
 * MonteCarloSystem implementation. All original functionality is preserved.
 */

#pragma once

// Include the new modular implementation
#include "montecarlo/MCMain.hpp"

// Re-export the MonteCarloSystem class for backward compatibility
// The class is now implemented in montecarlo/MCMain.hpp but the API remains the same
namespace pygcmc {
namespace system {
    using MonteCarloSystem = montecarlo::MCMain;
    
    // Re-export the nested types for compatibility
    using MovementMolecularInfo = montecarlo::MovementMolecularInfo;
} // namespace system
} // namespace pygcmc

