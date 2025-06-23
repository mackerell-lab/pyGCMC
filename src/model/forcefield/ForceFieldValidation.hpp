#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_VALIDATION_HPP
#define PYGCMC_MODEL_FORCEFIELD_VALIDATION_HPP

/**
 * @brief Complete validation implementation with modular design
 * 
 * This file serves as a convenience header that includes all validation
 * components for full backward compatibility.
 * 
 * The ForceFieldValidator class has been refactored into smaller files:
 * - ForceFieldValidationCore.hpp: Core validation class and basic methods
 * - ForceFieldValidationUtils.hpp: Utility methods and advanced validation
 * 
 * This maintains 100% backward compatibility while improving maintainability.
 */

// Include core validation implementation
#include "ForceFieldValidationCore.hpp"

// Include utility methods implementation
#include "ForceFieldValidationUtils.hpp"

// The complete ForceFieldValidator class is now available with all methods
// implemented through the included implementation files

#endif // PYGCMC_MODEL_FORCEFIELD_VALIDATION_HPP 