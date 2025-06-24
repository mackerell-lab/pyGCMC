#pragma once

#include "MCStructures.hpp"
#include <stdexcept>

/**
 * @file   MCOperations.hpp
 * @brief  Basic operations for Monte Carlo state management
 *
 * Simple operations without complex business logic.
 * Each operation focuses on a single responsibility.
 */

namespace pygcmc {
namespace model {
namespace montecarlo {

/**
 * @brief Essential operations for Monte Carlo state management
 */
class MCOperations {
public:
    // === Force Field Operations ===
    static void initializeForceField(MCForceField& ff, int total_types, int movement_types = 0) {
        ff.numTotalTypes = total_types;
        ff.numMovementTypes = movement_types;
        size_t matrix_size = total_types * total_types;
        ff.ljSigma.resize(matrix_size, 0.0f);
        ff.ljEps.resize(matrix_size, 0.0f);
    }

    static void setLJParams(MCForceField& ff, int type1, int type2, float sigma, float epsilon) {
        if (type1 >= ff.numTotalTypes || type2 >= ff.numTotalTypes || type1 < 0 || type2 < 0) {
            throw std::out_of_range("Type index out of range");
        }
        int index = type1 * ff.numTotalTypes + type2;
        ff.ljSigma[index] = sigma;
        ff.ljEps[index] = epsilon;
        
        if (type1 != type2) {
            int sym_index = type2 * ff.numTotalTypes + type1;
            ff.ljSigma[sym_index] = sigma;
            ff.ljEps[sym_index] = epsilon;
        }
    }

    // === Atom Operations ===
    static int addAtom(std::vector<MCAtom>& atoms, int& activeCount, const MCAtom& atom) {
        if (activeCount >= static_cast<int>(atoms.size())) {
            atoms.push_back(atom);
        } else {
            atoms[activeCount] = atom;
        }
        return activeCount++;
    }

    static void removeAtom(std::vector<MCAtom>& atoms, int& activeCount, int index) {
        if (index >= 0 && index < activeCount) {
            if (index < activeCount - 1) {
                atoms[index] = atoms[activeCount - 1];
            }
            activeCount--;
        }
    }

    // === Residue Operations ===
    static int addResidue(std::vector<MCResidue>& residues, int& activeCount, const MCResidue& residue) {
        if (activeCount >= static_cast<int>(residues.size())) {
            residues.push_back(residue);
        } else {
            residues[activeCount] = residue;
        }
        residues[activeCount].active = true;
        return activeCount++;
    }

    static void removeResidue(std::vector<MCResidue>& residues, int& activeCount, int index) {
        if (index >= 0 && index < activeCount) {
            residues[index].active = false;
            if (index < activeCount - 1) {
                residues[index] = residues[activeCount - 1];
                residues[index].active = true;
            }
            activeCount--;
        }
    }

    // === System Operations ===
    static void setBox(MCInfo& info, float x, float y, float z) {
        info.box[0] = x;
        info.box[1] = y;
        info.box[2] = z;
        info.volume = x * y * z;
    }

    static void updateStatistics(MCInfo::Statistics& stats, bool accepted, bool insertion = false, bool deletion = false) {
        stats.totalMoves++;
        if (accepted) {
            stats.acceptedMoves++;
        }
        
        if (insertion) {
            stats.insertionAttempts++;
            if (accepted) {
                stats.acceptedInsertions++;
            }
        }
        
        if (deletion) {
            stats.deletionAttempts++;
            if (accepted) {
                stats.acceptedDeletions++;
            }
        }
    }

    static void clearSystem(std::vector<MCAtom>& atoms, std::vector<MCResidue>& residues,
                          TypeMaps& residueTypes, TypeMaps& atomTypes,
                          std::vector<MCMovementResidueInfo>& movementResidues,
                          std::vector<int>& movementAtomTypes,
                          int& activeAtomCount, int& activeResidueCount,
                          int& numMovementAtomTypes,
                          EwaldEnergy& ewald_energy, MCInfo::Statistics& stats) {
        atoms.clear();
        residues.clear();
        residueTypes = TypeMaps{};
        atomTypes = TypeMaps{};
        movementResidues.clear();
        movementAtomTypes.clear();
        activeAtomCount = 0;
        activeResidueCount = 0;
        numMovementAtomTypes = 0;
        ewald_energy.reset();
        stats = MCInfo::Statistics{};
    }
};

} // namespace montecarlo
} // namespace model
} // namespace pygcmc