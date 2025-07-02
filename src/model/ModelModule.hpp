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
 * **atom** - Atomic-Level Data (Core + Inheritance Pattern)
 * - AtomCore.hpp: Foundation class with 25+ properties (coordinates, force field parameters, PDB fields)
 * - AtomMain.hpp: Extended class adding PDB formatting, validation, and comparison operators
 * 
 * **residue** - Residue-Level Composition (Core + Composition Pattern)
 * - ResidueCore.hpp: Pure data structures including secondary structure and disulfide bond support
 * - ResidueMain.hpp: Container class managing atom collections with validation and center-of-mass calculation
 * 
 * **molecule** - Molecular System Integration (Single File)
 * - MolecularMain.hpp: Data fusion class combining structural (PDB) and topological (PSF/TOP) information
 *   with standardized CMAP handling and comprehensive lookup mappings
 * 
 * **topology** - Molecular Connectivity (4-File Modular Pattern)
 * - TopologyStructures.hpp: Data structures for atoms, bonds, angles, dihedrals, and special features
 * - TopologyOperations.hpp: Static methods for adding topology elements with hierarchical validation
 * - TopologyQueries.hpp: Static methods for finding, checking, and counting topology elements
 * - TopologyMain.hpp: Main interface class delegating to operations and queries with private data storage
 * 
 * **forcefield** - Force Field Parameters (3-File Optimized Pattern)
 * - ForceFieldTypes.hpp: Parameter type definitions for LJ, bonded, and NBFIX interactions
 * - ForceFieldAccessors.hpp: Comprehensive inline implementations with smart key generation and symmetry
 * - ForceFieldMain.hpp: Main class declaration providing parameter access interface
 * 
 * **montecarlo** - Monte Carlo Simulation (2-File Efficient Pattern)
 * - MCStructures.hpp: Performance-optimized data structures for simulation state and type mapping
 * - MCMain.hpp: Complete MCState class with atom/residue management and statistics tracking
 * 
 * **param** - Simulation Parameters (4-File Modular Pattern)
 * - ParamStructures.hpp: Seven parameter structure definitions covering all GCMC simulation aspects
 * - ParamOperations.hpp: Static utility methods for parameter updates and derived value calculations
 * - ParamQueries.hpp: Validation, string conversion, and parameter query methods
 * - ParamMain.hpp: Main interface class with complete delegation to operations and queries
 * 
 * **structure** - Basic Structural Data (Single File)
 * - StructureMain.hpp: Simple container for atoms, residues, and secondary structure elements
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
 * - Optimal pattern selection: each module uses the best file organization for its complexity
 * - Well-designed: predictable organization with clear separation where beneficial
 * - Performance-oriented: efficient data structures with minimal overhead
 * - Backward compatibility: all original APIs preserved through type aliases
 * - Pattern diversity: 4-file (topology/param), 3-file (forcefield), 2-file (montecarlo/atom/residue), single-file (molecule/structure)
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
using AlphaTHoleParams = forcefield::AlphaTHoleParams;
using LonePairParams = forcefield::LonePairParams;
using AnisotropyParams = forcefield::AnisotropyParams;

/**
 * @brief Get model module version
 */
inline const char* getModelVersion() {
    return "1.0.0";
}

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_MODULE_HPP 