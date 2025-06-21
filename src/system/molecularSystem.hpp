/**
 * @file molecularSystem.hpp
 * @brief Legacy compatibility header for MolecularSystem
 * 
 * This file provides backward compatibility by including the new modular
 * MolecularSystem implementation. All original functionality is preserved.
 */

#pragma once

// Include the new modular implementation
#include "molecular/MolecularMain.hpp"

// Re-export the MolecularSystem class for backward compatibility
// The class is now implemented in molecular/MolecularMain.hpp but the API remains the same
namespace pygcmc {
namespace system {
    using MolecularSystem = molecular::MolecularMain;
} // namespace system
} // namespace pygcmc
