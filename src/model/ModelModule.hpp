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
 * **atom** - Atomic-Level Data Structures
 * - AtomCore.hpp: Core atom data structure with coordinates, mass, charge
 * - AtomMain.hpp: Complete Atom class with PDB format support and utilities
 * 
 * **residue** - Residue-Level Composition
 * - ResidueCore.hpp: Core residue data structure and atom management
 * - ResidueMain.hpp: Main interface with validation and utilities
 * 
 * **molecule** - Molecular System Management
 * - MolecularMain.hpp: Complete molecular system with multi-residue support
 * 
 * **topology** - Topology and Connectivity (follows 4-file pattern)
 * - TopologyStructures.hpp: Data structures (atoms, bonds, angles, dihedrals)
 * - TopologyOperations.hpp: Add/modify operations (static methods)
 * - TopologyQueries.hpp: Find/check/count operations (static methods)
 * - TopologyMain.hpp: Main interface (delegates to operations/queries)
 * 
 * **forcefield** - Force Field Parameters
 * - ForceFieldTypes.hpp: Parameter type definitions
 * - ForceFieldAccessors.hpp: Parameter access utilities
 * - ForceFieldMain.hpp: Main force field interface
 * 
 * **montecarlo** - Monte Carlo Simulation States
 * - MCStructures.hpp: Pure data structures (TypeMaps, MCInfo, MCAtom, MCResidue)
 * - MCMain.hpp: Main interface with complete MCState functionality
 * 
 * **param** - Global Simulation Parameters (follows 4-file pattern)
 * - ParamStructures.hpp: Parameter data structures
 * - ParamOperations.hpp: Parameter modification operations
 * - ParamQueries.hpp: Parameter query operations
 * - ParamMain.hpp: Comprehensive parameter management
 * 
 * **structure** - Structural Data Management
 * - StructureMain.hpp: High-level structural data organization
 * 
 * **Usage Examples:**
 * 
 * ```cpp
 * #include "model/ModelModule.hpp"
 * using namespace pygcmc::model;
 * 
 * // Create molecular components
 * auto molecular = std::make_shared<Molecular>();
 * auto topology = std::make_shared<Topology>();
 * auto mcstate = std::make_shared<MCState>();
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
 * - AI-friendly: focused responsibilities with predictable file organization
 * - Performance-oriented: efficient data structures with minimal overhead
 * - Consistent patterns: topology/ and param/ follow 4-file pattern (Structures/Operations/Queries/Main)
 */

// === Core Utilities ===
#include <iostream>
#include <limits>

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
using StandardCmap = molecule::StandardCmap;
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