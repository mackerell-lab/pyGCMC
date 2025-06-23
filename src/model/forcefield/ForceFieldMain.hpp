#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_MAIN_HPP
#define PYGCMC_MODEL_FORCEFIELD_MAIN_HPP

/**
 * @brief Complete ForceField implementation with modular design
 * 
 * This file serves as a convenience header that includes all ForceField
 * implementation components for full backward compatibility.
 * 
 * The ForceField class has been refactored into smaller, specialized files:
 * - ForceFieldCore.hpp: Core class definition and basic methods
 * - ForceFieldParameterOps.hpp: Parameter addition and retrieval implementations
 * - ForceFieldCheckerOps.hpp: Parameter existence and size checking implementations
 * - ForceFieldAnalysisOps.hpp: Analysis and validation implementations
 * 
 * This maintains 100% backward compatibility while improving maintainability.
 */

// Include core class definition
#include "ForceFieldCore.hpp"

// Include all method implementations
#include "ForceFieldParameterOps.hpp"
#include "ForceFieldCheckerOps.hpp"
#include "ForceFieldAnalysisOps.hpp"

// The complete ForceField class is now available with all methods implemented
// through the included implementation files

#endif // PYGCMC_MODEL_FORCEFIELD_MAIN_HPP 