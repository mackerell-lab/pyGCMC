#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_STORAGE_OPTIMIZED_HPP
#define PYGCMC_MODEL_FORCEFIELD_STORAGE_OPTIMIZED_HPP

/**
 * @brief Complete optimized storage implementation with modular design
 * 
 * This file serves as a convenience header that includes all optimized storage
 * components for full backward compatibility.
 * 
 * The ForceFieldStorageOptimized class has been refactored into smaller files:
 * - ForceFieldStorageCore.hpp: Core storage class and main operations
 * - ForceFieldStorageUtils.hpp: Utility methods and helper functions
 * 
 * This maintains 100% backward compatibility while improving maintainability.
 */

// Include core storage implementation
#include "ForceFieldStorageCore.hpp"

// Include utility methods implementation
#include "ForceFieldStorageUtils.hpp"

// The complete ForceFieldStorageOptimized class is now available with all methods
// implemented through the included implementation files

#endif // PYGCMC_MODEL_FORCEFIELD_STORAGE_OPTIMIZED_HPP 