#pragma once

#include "MCStructures.hpp"
#include "MCOperations.hpp"
#include <memory>
#include <optional>

/**
 * @file   MCMain.hpp
 * @brief  Main interface for Monte Carlo state management
 *
 * Complete Monte Carlo state class that maintains the same API as the original.
 * Provides direct access to all functionality while delegating to operations.
 */

namespace pygcmc {
namespace model {
namespace montecarlo {

/**
 * @brief Current state of the Monte Carlo system
 * 
 * Contains all information about the current state of the system,
 * including atoms, residues, type mappings, and force field parameters.
 */
struct MCState {
    // Core data arrays
    std::vector<MCAtom>    atoms;
    std::vector<MCResidue> residues;
    TypeMaps residueTypes;
    TypeMaps atomTypes;
    std::vector<MCMovementResidueInfo> movementResidues;
    std::vector<int> movementAtomTypes;
    int numMovementAtomTypes{0};
    int activeAtomCount{0};
    int activeResidueCount{0};
    MCInfo info;
    MCForceField forcefield;
    EwaldEnergy ewald_energy;

    // Constructor to initialize counts
    MCState() = default;

    // === Essential Operations ===
    int addAtom(const MCAtom& atom) {
        return MCOperations::addAtom(atoms, activeAtomCount, atom);
    }

    void removeAtom(int index) {
        MCOperations::removeAtom(atoms, activeAtomCount, index);
    }

    int addResidue(const MCResidue& residue) {
        return MCOperations::addResidue(residues, activeResidueCount, residue);
    }

    void removeResidue(int index) {
        MCOperations::removeResidue(residues, activeResidueCount, index);
    }

    // === System Management ===
    void setTemperature(float temperature) { 
        info.setTemperature(temperature); 
    }

    void setBoxDimensions(float x, float y, float z) { 
        MCOperations::setBox(info, x, y, z); 
    }

    void setupForceField(int totalTypes, int movementTypes = 0) { 
        MCOperations::initializeForceField(forcefield, totalTypes, movementTypes); 
        numMovementAtomTypes = movementTypes; 
    }

    void setLJParameters(int type1, int type2, float sigma, float epsilon) { 
        MCOperations::setLJParams(forcefield, type1, type2, sigma, epsilon); 
    }

    // === Type Management ===
    int getOrAddAtomType(const std::string& type) { 
        return atomTypes.getOrAddType(type); 
    }

    int getOrAddResidueType(const std::string& type) { 
        return residueTypes.getOrAddType(type); 
    }

    // === Statistics ===
    void incrementMoveStats(bool accepted) { 
        MCOperations::updateStatistics(info.stats, accepted);
    }

    void incrementInsertionStats(bool accepted) { 
        MCOperations::updateStatistics(info.stats, accepted, true, false);
    }

    void incrementDeletionStats(bool accepted) { 
        MCOperations::updateStatistics(info.stats, accepted, false, true);
    }

    // === Basic Queries ===
    int getActiveAtomCount() const { return activeAtomCount; }
    int getActiveResidueCount() const { return activeResidueCount; }
    int getTotalAtomCount() const { return static_cast<int>(atoms.size()); }
    int getTotalResidueCount() const { return static_cast<int>(residues.size()); }

    std::optional<int> findAtomByType(int type) const {
        for (int i = 0; i < activeAtomCount; ++i) {
            if (atoms[i].type == type) return i;
        }
        return std::nullopt;
    }

    std::optional<int> findResidueByType(int type) const {
        for (int i = 0; i < activeResidueCount; ++i) {
            if (residues[i].type == type && residues[i].active) return i;
        }
        return std::nullopt;
    }

    bool isValid() const {
        if (activeAtomCount < 0 || activeResidueCount < 0) return false;
        if (activeAtomCount > static_cast<int>(atoms.size())) return false;
        if (activeResidueCount > static_cast<int>(residues.size())) return false;
        
        for (int i = 0; i < activeResidueCount; ++i) {
            const auto& residue = residues[i];
            if (!residue.active) continue;
            if (residue.atomStart < 0 || residue.atomCount <= 0) return false;
            if (residue.atomStart + residue.atomCount > activeAtomCount) return false;
        }
        
        return true;
    }

    void clear() {
        MCOperations::clearSystem(atoms, residues, residueTypes, atomTypes,
                                        movementResidues, movementAtomTypes,
                                        activeAtomCount, activeResidueCount, numMovementAtomTypes,
                                        ewald_energy, info.stats);
    }

    std::unique_ptr<MCState> clone() const {
        return std::make_unique<MCState>(*this);
    }
};

} // namespace montecarlo

// Backward compatibility: provide the MC classes in the model namespace
using MCState = montecarlo::MCState;
using MCAtom = montecarlo::MCAtom;
using MCResidue = montecarlo::MCResidue;
using MCInfo = montecarlo::MCInfo;
using MCForceField = montecarlo::MCForceField;
using MCMovementResidueInfo = montecarlo::MCMovementResidueInfo;
using TypeMaps = montecarlo::TypeMaps;

} // namespace model
} // namespace pygcmc