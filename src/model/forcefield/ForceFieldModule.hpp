#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_MODULE_HPP
#define PYGCMC_MODEL_FORCEFIELD_MODULE_HPP

/**
 * @file ForceFieldModule.hpp
 * @brief Force Field Module - Complete Force Field Parameter Management
 * 
 * This module provides comprehensive force field parameter management for
 * molecular simulations. It includes all CHARMM force field parameter types
 * with full validation, analysis, and backward compatibility.
 * 
 * **Module Components:**
 * 
 * **ForceFieldParams.hpp** - Parameter Structures and Storage
 * - All parameter structures (LJ, Bond, Angle, Dihedral, Improper, NBFIX)
 * - NonbondedParams with complete CHARMM documentation
 * - ForceFieldStorage container class
 * - ParamKeyUtils for consistent parameter keys
 * - Statistics and completeness check structures
 * 
 * **ForceFieldCore.hpp** - Core Management Operations
 * - ForceFieldManager class for parameter operations
 * - Parameter addition with full validation
 * - Parameter retrieval with error checking
 * - Existence checks and size queries
 * - Bidirectional angle parameter lookup
 * 
 * **ForceFieldUtils.hpp** - Analysis and Utility Functions
 * - ForceFieldUtils class for analysis operations
 * - Validation and completeness checking
 * - Statistical analysis and summaries
 * - Consistency validation
 * - Missing parameter detection
 * 
 * **ForceFieldMain.hpp** - Main Interface Class
 * - Complete ForceField class with all functionality
 * - Backward compatibility with original API
 * - Composition-based architecture
 * - IValidatable interface implementation
 * 
 * **Key Features:**
 * - Complete CHARMM force field support
 * - Full backward compatibility with existing code
 * - Modular architecture for easy extension
 * - Comprehensive parameter validation
 * - Detailed analysis and statistics
 * - AI-friendly modular design
 * - Thread-safe read operations
 * 
 * **Usage Examples:**
 * 
 * ```cpp
 * #include "model/forcefield/ForceFieldModule.hpp"
 * using namespace pygcmc::model;
 * 
 * // 1. Basic force field creation
 * ForceField ff;
 * ff.add_atom_mass("CA", 12.011);
 * ff.add_lj_params("CA", 0.070, 1.992);
 * 
 * // 2. Parameter retrieval
 * double mass = ff.get_atom_mass("CA");
 * const auto& lj = ff.get_lj_params("CA");
 * 
 * // 3. Existence checking
 * if (ff.has_atom_mass("CA")) {
 *     std::cout << "CA mass available" << std::endl;
 * }
 * 
 * // 4. Analysis and validation
 * if (ff.is_valid()) {
 *     auto stats = ff.get_statistics();
 *     std::cout << "Atom types: " << stats.num_atom_types << std::endl;
 * }
 * 
 * // 5. Completeness checking
 * std::set<std::string> required_types = {"CA", "CB", "N", "O"};
 * auto result = ff.check_completeness(required_types);
 * if (!result.is_complete) {
 *     std::cout << result.summary << std::endl;
 * }
 * 
 * // 6. Detailed analysis
 * std::string summary = ff.get_detailed_statistics();
 * std::cout << summary << std::endl;
 * ```
 * 
 * **Backward Compatibility:**
 * All original ForceField APIs are preserved. Existing code using the
 * force field classes requires no changes when upgrading to this module.
 * 
 * **Performance Notes:**
 * - Parameter lookups use std::map for O(log n) complexity
 * - Storage is optimized for memory efficiency
 * - Analysis operations cache results where possible
 * - Thread-safe for concurrent read operations
 * 
 * **Integration with Other Modules:**
 * - Used by IO parsers (prmParser, strParser)
 * - Integrated with topology building
 * - Compatible with energy calculation modules
 * - Supports Python bindings via direct map access
 */

// === Core Components ===
#include "ForceFieldParams.hpp"
#include "ForceFieldCore.hpp"
#include "ForceFieldUtils.hpp"
#include "ForceFieldMain.hpp"

namespace pygcmc {
namespace model {

// === Backward Compatibility Type Aliases ===
// These provide seamless compatibility with existing code
using ForceField = forcefield::ForceField;
using NonbondedParams = forcefield::NonbondedParams;
using LJParams = forcefield::LJParams;
using BondParams = forcefield::BondParams;
using AngleParams = forcefield::AngleParams;
using DihedralParams = forcefield::DihedralParams;
using ImproperParams = forcefield::ImproperParams;
using NBFIXParams = forcefield::NBFIXParams;

// Additional utility types
using ForceFieldStats = forcefield::ForceFieldStats;
using CompletenessResult = forcefield::CompletenessResult;
using ForceFieldStorage = forcefield::ForceFieldStorage;
using ParamKeyUtils = forcefield::ParamKeyUtils;

namespace info {
/**
 * @brief Get module information
 */
inline std::string get_forcefield_module_info() {
    return "ForceField Module v1.0 - Complete CHARMM Force Field Parameter Management\n"
           "- Parameter structures: LJ, Bond, Angle, Dihedral, Improper, NBFIX\n"
           "- Core management with full validation\n"
           "- Analysis and utility functions\n"
           "- Backward compatible main interface\n"
           "- Thread-safe read operations\n"
           "- AI-friendly modular architecture";
}

/**
 * @brief Get supported parameter types
 */
inline std::vector<std::string> get_supported_parameter_types() {
    return {
        "Atom Masses",
        "Lennard-Jones (LJ) Parameters", 
        "NBFIX Override Parameters",
        "Bond Parameters",
        "Angle Parameters (with Urey-Bradley)",
        "Dihedral Parameters (multiple per type)",
        "Improper Parameters",
        "Nonbonded Control Parameters"
    };
}
} // namespace info

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_MODULE_HPP 