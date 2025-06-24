#pragma once

#ifndef PYGCMC_MODEL_MODULE_HPP
#define PYGCMC_MODEL_MODULE_HPP

/**
 * @file ModelModule.hpp
 * @brief Model Module - Unified Entry Point for All Data Structures
 * 
 * This is the ONLY header file you need to include to access all molecular modeling
 * data structures and utilities in the GCMC simulation framework. The module provides
 * a complete data layer for molecular systems, force fields, and Monte Carlo states.
 * 
 * **Module Organization:**
 * 
 * **common/ directory** - Core Interfaces and Utilities
 * - ModelInterface.hpp: Abstract interfaces (ICloneable, ISerializable, IValidatable)
 * - ModelConstants.hpp: Physical constants, unit conversions, and molecular data
 * - ModelUtils.hpp: Mathematical utilities, hashing, and validation helpers
 * 
 * **atom/ directory** - Atomic-Level Data Structures
 * - AtomMain.hpp: Complete Atom class with PDB format support and utilities
 * - AtomCore.hpp: Core atom data structure with coordinates, mass, charge
 * - Features: Distance calculations, PDB record generation, atom identification
 * 
 * **residue/ directory** - Residue-Level Composition
 * - ResidueMain.hpp: Main interface with atom management and validation
 * - ResidueComposite.hpp: High-level residue composition and geometry
 * - ResidueValidator.hpp: Comprehensive validation for backbone, charge, mass
 * - Features: Backbone/sidechain selection, geometry validation, statistics
 * 
 * **molecule/ directory** - Molecular System Management
 * - MolecularMain.hpp: Complete molecular system with multi-residue support
 * - MolecularComposite.hpp: System-level composition and coordinate management
 * - MolecularUtils.hpp: Utilities for molecular manipulation and analysis
 * - Features: Center of mass, system statistics, molecular transformations
 * 
 * **topology/ directory** - Force Field and Topology
 * - TopologyMain.hpp: Topology data structures (bonds, angles, dihedrals)
 * - TopologyForceField.hpp: Force field parameters (LJ, bonded, NBFIX)
 * - Features: CHARMM force field support, parameter validation, topology building
 * 
 * **montecarlo/ directory** - Monte Carlo Simulation States
 * - MCStructures.hpp: Pure data structures (TypeMaps, MCInfo, MCAtom, MCResidue)
 * - MCOperations.hpp: Essential operations (add/remove atoms/residues, statistics)
 * - MCMain.hpp: Main interface with complete MCState functionality
 * - Features: Residue insertion/deletion, energy state management, move validation
 * 
 * **param/ directory** - Global Simulation Parameters
 * - ParamMain.hpp: Comprehensive parameter management for all simulation aspects
 * - Features: Basic info, space info, energy parameters, fragment management
 * 
 * **structure/ directory** - Structural Data Management
 * - StructureMain.hpp: High-level structural data organization
 * - Features: Multi-chain structures, structural validation, format conversion
 * 
 * **Usage Examples:**
 * 
 * ```cpp
 * #include "model/ModelModule.hpp"
 * using namespace pygcmc::model;
 * 
 * // Create molecular components
 * auto molecular = std::make_shared<molecule::Molecular>();
 * 
 * // Validate system
 * if (validate_molecular_system(*molecular)) {
 *     std::cout << "System is valid" << std::endl;
 * }
 * 
 * // Version info
 * std::cout << "Version: " << getModelVersion() << std::endl;
 * ```
 * 
 * **Backward Compatibility:**
 * All original APIs (Atom, Residue, Molecular, etc.) are preserved through
 * direct type aliases. Existing code requires no changes when upgrading.
 * 
 * **Design Philosophy:**
 * - Single header inclusion for all data structure functionality
 * - Modular architecture with clear separation of concerns
 * - Each component is independently testable and maintainable
 * - AI-friendly: focused responsibilities with comprehensive documentation
 * - Performance-oriented: efficient data structures with minimal overhead
 * - Follows system module naming patterns for consistency
 */

// === Core Interfaces and Utilities ===
#include <iostream>
#include "common/ModelInterface.hpp"
#include "common/ModelConstants.hpp"
#include "common/ModelUtils.hpp"

// === Data Structure Layers ===
#include "atom/AtomMain.hpp"
#include "residue/ResidueMain.hpp"
#include "molecule/MolecularMain.hpp"
#include "structure/StructureMain.hpp"

// === Topology and Force Fields ===
#include "topology/TopologyMain.hpp"
#include "forcefield/ForceFieldMain.hpp"

// === Monte Carlo Simulation ===
#include "montecarlo/MCMain.hpp"

// === Parameters and Configuration ===
#include "param/ParamMain.hpp"

namespace pygcmc {
namespace model {

// Backward compatibility type aliases
using Atom = atom::Atom;
using Residue = residue::Residue;
using Molecular = molecule::Molecular;
using Structure = structure::Structure;
using Topology = topology::Topology;
using ForceField = forcefield::ForceField;
using MCState = montecarlo::MCState;
using MCAtom = montecarlo::MCAtom;
using MCResidue = montecarlo::MCResidue;
using MCInfo = montecarlo::MCInfo;
using MCForceField = montecarlo::MCForceField;
using Param = param::Param;

// Force field type aliases
using NonbondedParams = forcefield::NonbondedParams;
using LJParams = forcefield::LJParams;
using BondParams = forcefield::BondParams;
using AngleParams = forcefield::AngleParams;
using DihedralParams = forcefield::DihedralParams;
using ImproperParams = forcefield::ImproperParams;
using NBFIXParams = forcefield::NBFIXParams;

/**
 * @brief Get model module version
 */
inline const char* getModelVersion() {
    return "1.0.0";
}

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_MODULE_HPP 