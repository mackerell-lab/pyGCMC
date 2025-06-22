#pragma once

#ifndef PYGCMC_MODEL_COMMON_ALIASES_HPP
#define PYGCMC_MODEL_COMMON_ALIASES_HPP

#include "../atom/AtomMain.hpp"
#include "../residue/ResidueMain.hpp"
#include "../molecule/MolecularMain.hpp"
#include "../structure/StructureMain.hpp"
#include "../topology/TopologyMain.hpp"
#include "../topology/ForceFieldMain.hpp"
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
using ForceField = topology::ForceField;

// === Monte Carlo System Aliases ===
using MCState = montecarlo::MCState;

// === Main Parameter Class with Nested Compatibility ===
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

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_COMMON_ALIASES_HPP 