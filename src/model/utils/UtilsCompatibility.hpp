#pragma once

#ifndef PYGCMC_MODEL_UTILS_COMPATIBILITY_HPP
#define PYGCMC_MODEL_UTILS_COMPATIBILITY_HPP

#include "UtilsFactory.hpp"
#include "UtilsValidation.hpp"
#include "UtilsTesting.hpp"
#include "UtilsInfo.hpp"
#include "../atom/AtomMain.hpp"
#include "../residue/ResidueMain.hpp"
#include "../molecule/MolecularMain.hpp"
#include "../structure/StructureMain.hpp"
#include "../topology/TopologyMain.hpp"
#include "../forcefield/ForceFieldModule.hpp"
#include "../montecarlo/MCMain.hpp"
#include "../param/ParamMain.hpp"

namespace pygcmc {
namespace model {

/**
 * @brief Backward Compatibility Type Aliases
 * 
 * These aliases maintain 100% compatibility with existing code that was using
 * individual header files before the modular refactoring.
 */

// === Primary Data Structure Aliases ===
using Atom = atom::Atom;
using Residue = residue::Residue;
using Molecular = molecule::Molecular;
using Structure = structure::Structure;

// === Topology System Aliases ===
using Topology = topology::Topology;

// === Force Field System Aliases ===
using ForceField = forcefield::ForceField;
using NonbondedParams = forcefield::NonbondedParams;
using LJParams = forcefield::LJParams;
using BondParams = forcefield::BondParams;
using AngleParams = forcefield::AngleParams;
using DihedralParams = forcefield::DihedralParams;
using ImproperParams = forcefield::ImproperParams;
using NBFIXParams = forcefield::NBFIXParams;
using ForceFieldStats = forcefield::ForceFieldStats;
using CompletenessResult = forcefield::CompletenessResult;

// === Monte Carlo System Aliases ===
using MCState = montecarlo::MCState;

// === Shared Pointer Aliases for Convenience ===
using AtomPtr = std::shared_ptr<Atom>;
using ResiduePtr = std::shared_ptr<Residue>;
using MolecularPtr = std::shared_ptr<Molecular>;
using StructurePtr = std::shared_ptr<Structure>;
using TopologyPtr = std::shared_ptr<Topology>;
using ForceFieldPtr = std::shared_ptr<ForceField>;
using MCStatePtr = std::shared_ptr<MCState>;

/**
 * @brief Main Parameter Class with Nested Compatibility
 * 
 * This class extends the parameter system with nested type aliases
 * for Python binding compatibility and legacy code support.
 */
class Param : public param::Param {
public:
    // Nested type aliases for Python binding compatibility
    using BasicInfo = param::BasicInfo;
    using SpaceInfo = param::SpaceInfo;
    using MCInfo = param::MCParams;
    using EnergyInfo = param::EnergyInfo;
    using FragmentInfo = param::FragmentInfo;
    using BiasInfo = param::BiasInfo;
    using FileInfo = param::FileInfo;
    
    // Inherit all constructors and functionality
    using param::Param::Param;
};

// === Legacy Function Aliases ===

/**
 * @brief Legacy factory function aliases for backward compatibility
 */
namespace legacy {
    
    // Import factory functions with legacy names
    using namespace factory;
    
    // Legacy validation functions
    using namespace validation;
    
    // Legacy testing functions  
    using namespace testing;
    
} // namespace legacy

/**
 * @brief Compatibility namespace that includes all utilities
 * 
 * This namespace provides a single point of access to all utilities
 * for backward compatibility with older code.
 */
namespace compat {
    
    // Include all factory functions
    using namespace factory;
    
    // Include all validation functions
    using namespace validation;
    
    // Include all testing functions
    using namespace testing;
    
    // Include all information functions
    using namespace info;
    
} // namespace compat

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_UTILS_COMPATIBILITY_HPP 