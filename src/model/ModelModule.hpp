#pragma once

#ifndef PYGCMC_MODEL_MODULE_HPP
#define PYGCMC_MODEL_MODULE_HPP

/**
 * @file ModelModule.hpp
 * @brief Model Module - Data Layer Unified Entry Point
 * 
 * This is the single entry point for the entire Model module, providing access to
 * all data structures and utilities needed for molecular modeling and GCMC simulation.
 * 
 * The Model module is organized into several functional areas:
 * 
 * 📁 common/     - Common interfaces, constants, and utilities
 * 📁 atom/       - Atom-level data structures and operations  
 * 📁 residue/    - Residue-level composition and validation
 * 📁 molecule/   - Molecular system management and utilities
 * 📁 topology/   - Topology and force field parameters
 * 📁 montecarlo/ - Monte Carlo state and simulation data
 * 📁 param/      - Global simulation parameters
 * 
 * Design Principles:
 * - AI-Friendly: Small, focused files (≤200 lines each)
 * - Maximum Extensibility: Clear separation of concerns
 * - Consistent Naming: Main/Composite/Core/Utils pattern
 * - Backward Compatibility: All original APIs preserved
 */

// === Core Interfaces ===
#include "common/ModelInterface.hpp"
#include "common/ModelConstants.hpp"
#include "common/ModelUtils.hpp"

// === Atom Layer ===
#include "atom/AtomMain.hpp"

// === Residue Layer ===
#include "residue/ResidueMain.hpp"

// === Molecule Layer ===
#include "molecule/MolecularMain.hpp"

// === Topology & Force Field ===
#include "topology/TopologyMain.hpp"
#include "topology/ForceFieldMain.hpp"

// === Monte Carlo Simulation ===
#include "montecarlo/MCMain.hpp"

// === Parameters ===
#include "param/ParamMain.hpp"

// === Structure ===
#include "structure/StructureMain.hpp"

// === Utility & Validation ===
#include "utils/ModelFactory.hpp"

namespace pygcmc {
namespace model {

/**
 * @brief Model Module Information
 */
namespace info {
    constexpr const char* VERSION = "2.0.0";
    constexpr const char* BUILD_DATE = __DATE__;
    
    constexpr int TOTAL_COMPONENTS = 7;
    constexpr const char* COMPONENTS[] = {
        "common", "atom", "residue", "molecule", 
        "topology", "montecarlo", "param"
    };
}

// Detailed usage examples have been moved to docs/ModelExamples.md.

/**
 * @brief Backward Compatibility Type Aliases
 * 
 * These aliases maintain compatibility with code that was using the individual
 * header files (atom.hpp, residue.hpp, etc.) before the refactoring.
 */

// === Atom Module Aliases ===
using Atom = atom::Atom;

// === Residue Module Aliases ===
using Residue = residue::Residue;

// === Molecular Module Aliases ===
using Molecular = molecule::Molecular;

// === Structure Module Aliases ===
using Structure = structure::Structure;

// === Topology Module Aliases ===
using Topology = topology::Topology;
using TopologyAtom = topology::TopologyAtom;
using TopologyResidue = topology::TopologyResidue;
using TopologySegment = topology::TopologySegment;
using TopologyBond = topology::TopologyBond;
using TopologyAngle = topology::TopologyAngle;
using TopologyDihedral = topology::TopologyDihedral;
using TopologyDonor = topology::TopologyDonor;
using TopologyAcceptor = topology::TopologyAcceptor;
using TopologyGroup = topology::TopologyGroup;
using TopologyCmap = topology::TopologyCmap;

// === ForceField Module Aliases ===
using ForceField = topology::ForceField;
using NonbondedParams = topology::NonbondedParams;
using LJParams = topology::LJParams;
using BondParams = topology::BondParams;
using AngleParams = topology::AngleParams;
using DihedralParams = topology::DihedralParams;
using ImproperParams = topology::ImproperParams;
using NBFIXParams = topology::NBFIXParams;

// === Monte Carlo Module Aliases ===
using MCState = montecarlo::MCState;
using MCResidue = montecarlo::MCResidue;
using MCAtom = montecarlo::MCAtom;

// === Parameter Module Aliases ===
using BasicInfo = param::BasicInfo;
using SpaceInfo = param::SpaceInfo;
using EnergyInfo = param::EnergyInfo;
using FragmentInfo = param::FragmentInfo;
using BiasInfo = param::BiasInfo;
using FileInfo = param::FileInfo;

// Main Param class with nested type compatibility
class Param : public param::Param {
public:
    // Re-export types as nested types for Python binding compatibility
    using BasicInfo = param::BasicInfo;
    using SpaceInfo = param::SpaceInfo;
    using MCInfo = param::MCParams;  // Nested alias for backward compatibility
    using EnergyInfo = param::EnergyInfo;
    using FragmentInfo = param::FragmentInfo;
    using BiasInfo = param::BiasInfo;
    using FileInfo = param::FileInfo;
    
    // Inherit all constructors and methods
    using param::Param::Param;
};

} // namespace model
} // namespace pygcmc

/**
 * @brief Module Version and Build Information
 */
#define PYGCMC_MODEL_VERSION_MAJOR 2
#define PYGCMC_MODEL_VERSION_MINOR 0
#define PYGCMC_MODEL_VERSION_PATCH 0
#define PYGCMC_MODEL_VERSION_STRING "2.0.0"

#endif // PYGCMC_MODEL_MODULE_HPP 