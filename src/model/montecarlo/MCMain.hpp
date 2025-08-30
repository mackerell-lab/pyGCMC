#pragma once

#include "MCStructures.hpp"
#include <memory>
#include <optional>
#include <stdexcept>

/**
 * @file   MCMain.hpp
 * @brief  Main interface for Monte Carlo state management
 */

namespace pygcmc {
namespace model {
namespace montecarlo {

/**
 * @brief Current state of the Monte Carlo system
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
    std::vector<double> periodicBox;  // Box dimensions [x, y, z] in nm

    MCState() = default;

    // === Essential Operations ===
    int addAtom(const MCAtom& atom) {
        if (activeAtomCount >= static_cast<int>(atoms.size())) {
            atoms.push_back(atom);
        } else {
            atoms[activeAtomCount] = atom;
        }
        return activeAtomCount++;
    }

    void removeAtom(int index) {
        if (index >= 0 && index < activeAtomCount) {
            if (index < activeAtomCount - 1) {
                atoms[index] = atoms[activeAtomCount - 1];
            }
            activeAtomCount--;
        }
    }

    int addResidue(const MCResidue& residue) {
        if (activeResidueCount >= static_cast<int>(residues.size())) {
            residues.push_back(residue);
        } else {
            residues[activeResidueCount] = residue;
        }
        residues[activeResidueCount].active = true;
        return activeResidueCount++;
    }

    void removeResidue(int index) {
        if (index >= 0 && index < activeResidueCount) {
            residues[index].active = false;
            if (index < activeResidueCount - 1) {
                residues[index] = residues[activeResidueCount - 1];
                residues[index].active = true;
            }
            activeResidueCount--;
        }
    }

    // === System Management ===
    void setTemperature(float temperature) { 
        info.setTemperature(temperature); 
    }

    void setBoxDimensions(float x, float y, float z) { 
        info.box[0] = x; info.box[1] = y; info.box[2] = z;
        info.volume = x * y * z;
    }

    void setupForceField(int totalTypes, int movementTypes = 0) { 
        forcefield.numTotalTypes = totalTypes;
        forcefield.numMovementTypes = movementTypes;
        size_t matrix_size = totalTypes * totalTypes;
        forcefield.ljSigma.resize(matrix_size, 0.0f);
        forcefield.ljEps.resize(matrix_size, 0.0f);
        numMovementAtomTypes = movementTypes; 
    }

    void setLJParameters(int type1, int type2, float sigma, float epsilon) { 
        if (type1 >= forcefield.numTotalTypes || type2 >= forcefield.numTotalTypes || type1 < 0 || type2 < 0) {
            throw std::out_of_range("Type index out of range");
        }
        int index = type1 * forcefield.numTotalTypes + type2;
        forcefield.ljSigma[index] = sigma;
        forcefield.ljEps[index] = epsilon;
        if (type1 != type2) {
            int sym_index = type2 * forcefield.numTotalTypes + type1;
            forcefield.ljSigma[sym_index] = sigma;
            forcefield.ljEps[sym_index] = epsilon;
        }
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
        info.stats.totalMoves++;
        if (accepted) info.stats.acceptedMoves++;
    }

    void incrementInsertionStats(bool accepted) { 
        info.stats.insertionAttempts++;
        if (accepted) info.stats.acceptedInsertions++;
    }

    void incrementDeletionStats(bool accepted) { 
        info.stats.deletionAttempts++;
        if (accepted) info.stats.acceptedDeletions++;
    }

    // === Basic Queries ===
    int getActiveAtomCount() const { return activeAtomCount; }
    int getActiveResidueCount() const { return activeResidueCount; }

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