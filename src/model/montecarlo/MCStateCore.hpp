#pragma once

#ifndef PYGCMC_MODEL_MONTECARLO_STATE_CORE_HPP
#define PYGCMC_MODEL_MONTECARLO_STATE_CORE_HPP

#include "MCStructures.hpp"
#include "../common/ModelInterface.hpp"
#include <vector>

namespace pygcmc {
namespace model {
namespace montecarlo {

/**
 * @brief Simplified core state container for Monte Carlo system
 * 
 * Pure data container that holds system state with minimal methods.
 * Operations are handled by MCOperations and queries by MCQueries.
 */
class MCStateCore : public common::IValidatable {
public:
    MCStateCore() : numMovementAtomTypes(0), activeAtomCount(0), activeResidueCount(0) {}

    // IValidatable interface - delegates to MCQueries
    bool is_valid() const override;

    // Core data arrays
    std::vector<MCAtom>    atoms;      ///< Global atom array
    std::vector<MCResidue> residues;   ///< Global residue array

    // Type mappings
    TypeMaps residueTypes;  ///< Residue type mappings
    TypeMaps atomTypes;     ///< Atom type mappings

    // Movement molecule info
    std::vector<MCMovementResidueInfo> movementResidues;  ///< Info for movement residues
    std::vector<int> movementAtomTypes;  ///< Atom types from movement molecules
    int numMovementAtomTypes;  ///< Number of atom types from movement molecules

    // Active counts for swap-and-pop management
    int activeAtomCount;     ///< Current active atom count
    int activeResidueCount;  ///< Current active residue count

    // System info and force field
    MCInfo info;            ///< System information
    MCForceField forcefield;  ///< Force field parameters

    // Energy components
    EwaldEnergy ewald_energy;  ///< Ewald energy breakdown

    // Basic getters - no business logic
    int getActiveAtomCount() const { return activeAtomCount; }
    int getActiveResidueCount() const { return activeResidueCount; }
    int getTotalAtomCount() const { return static_cast<int>(atoms.size()); }
    int getTotalResidueCount() const { return static_cast<int>(residues.size()); }
};

// Implementation of is_valid() - defined here to avoid circular dependencies
inline bool MCStateCore::is_valid() const {
    // Check basic consistency
    if (activeAtomCount < 0 || activeResidueCount < 0) return false;
    if (activeAtomCount > static_cast<int>(atoms.size())) return false;
    if (activeResidueCount > static_cast<int>(residues.size())) return false;
    
    // Check residue-atom consistency
    for (int i = 0; i < activeResidueCount; ++i) {
        const auto& residue = residues[i];
        if (!residue.active) continue;
        if (residue.atomStart < 0 || residue.atomCount <= 0) return false;
        if (residue.atomStart + residue.atomCount > activeAtomCount) return false;
    }
    
    return true;
}

} // namespace montecarlo
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_MONTECARLO_STATE_CORE_HPP 