#pragma once

#ifndef PYGCMC_MODEL_MONTECARLO_OPERATIONS_HPP
#define PYGCMC_MODEL_MONTECARLO_OPERATIONS_HPP

#include "MCStructures.hpp"
#include <stdexcept>

namespace pygcmc {
namespace model {
namespace montecarlo {

/**
 * @brief Static operations for Monte Carlo state management
 */
class MCOperations {
public:
    // === MCInfo Operations ===
    static void setTemperature(MCInfo& info, float temperature) {
        if (temperature <= 0.0f) {
            throw std::invalid_argument("Temperature must be positive");
        }
        info.beta = 1.0f / (MCInfo::BOLTZMANN * temperature);
    }

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

    // === MCForceField Operations ===
    static void initializeForceField(MCForceField& ff, int total_types, int movement_types = 0) {
        ff.numTotalTypes = total_types;
        ff.numMovementTypes = movement_types;
        size_t matrix_size = total_types * total_types;
        ff.ljSigma.resize(matrix_size, 0.0f);
        ff.ljEps.resize(matrix_size, 0.0f);
    }

    static void setLJParams(MCForceField& ff, int type1, int type2, float sigma, float epsilon) {
        if (type1 >= ff.numTotalTypes || type2 >= ff.numTotalTypes) {
            throw std::out_of_range("Type index out of range");
        }
        int index = type1 * ff.numTotalTypes + type2;
        ff.ljSigma[index] = sigma;
        ff.ljEps[index] = epsilon;
        
        // Set symmetric entry
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
            // Swap with last active atom and decrease count
            if (index < activeCount - 1) {
                atoms[index] = atoms[activeCount - 1];
            }
            activeCount--;
        }
    }

    static void updateAtomPosition(MCAtom& atom, float x, float y, float z) {
        atom.setPosition(x, y, z);
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
            // Mark as inactive
            residues[index].active = false;
            
            // Swap with last active residue and decrease count
            if (index < activeCount - 1) {
                residues[index] = residues[activeCount - 1];
                residues[index].active = true;
            }
            activeCount--;
        }
    }

    static void updateResidueEnergy(MCResidue& residue, float vdw_energy, float elec_energy) {
        residue.energy_vdw = vdw_energy;
        residue.energy_elec = elec_energy;
    }

    static void updateResidueCenter(MCResidue& residue, float x, float y, float z) {
        residue.setCenter(x, y, z);
    }

    // === Movement Residue Operations ===
    static void addMovementResidue(std::vector<MCMovementResidueInfo>& movementResidues, 
                                 const std::string& resName, int startIndex, int totalCount) {
        MCMovementResidueInfo info(startIndex, 0, totalCount, resName);
        movementResidues.push_back(info);
    }

    static void updateMovementResidueCount(std::vector<MCMovementResidueInfo>& movementResidues, 
                                         size_t index, int activeCount) {
        if (index < movementResidues.size()) {
            movementResidues[index].activeCount = activeCount;
        }
    }

    // === Ewald Energy Operations ===
    static void updateEwaldEnergy(EwaldEnergy& ewald, double real_space, double reciprocal, double self_energy) {
        ewald.real_space = real_space;
        ewald.reciprocal = reciprocal;
        ewald.self = self_energy;
        ewald.updateTotal();
    }

    static void resetEwaldEnergy(EwaldEnergy& ewald) {
        ewald.reset();
    }

    // === Type Operations ===
    static int getOrAddType(TypeMaps& typeMaps, const std::string& type) {
        return typeMaps.getOrAddType(type);
    }

    static void clearTypes(TypeMaps& typeMaps) {
        typeMaps.clear();
    }

    // === System Operations ===
    static void reserveCapacity(std::vector<MCAtom>& atoms, std::vector<MCResidue>& residues, 
                               MCInfo& info, int max_atoms, int max_residues) {
        atoms.reserve(max_atoms);
        residues.reserve(max_residues);
        info.maxAtoms = max_atoms;
        info.maxResidues = max_residues;
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
        residueTypes.clear();
        atomTypes.clear();
        movementResidues.clear();
        movementAtomTypes.clear();
        activeAtomCount = 0;
        activeResidueCount = 0;
        numMovementAtomTypes = 0;
        ewald_energy.reset();
        stats.reset();
    }
};

} // namespace montecarlo
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_MONTECARLO_OPERATIONS_HPP