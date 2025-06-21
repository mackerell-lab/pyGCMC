#pragma once

#ifndef PYGCMC_MODEL_MONTECARLO_MAIN_HPP
#define PYGCMC_MODEL_MONTECARLO_MAIN_HPP

#include "MCStateCore.hpp"
#include "../common/ModelUtils.hpp"
#include <string>
#include <sstream>

namespace pygcmc {
namespace model {
namespace montecarlo {

/**
 * @brief Complete Monte Carlo state class with full functionality and backward compatibility
 * This class maintains the same API as the original montecarlo.hpp MCState
 */
class MCState : public MCStateCore {
public:
    MCState() = default;
    ~MCState() = default;

    // Backward compatibility: provide direct access to structure components
    // These reference the inherited members from MCStateCore

    // Direct access methods for compatibility
    std::vector<MCAtom>& getAtoms() { return atoms; }
    const std::vector<MCAtom>& getAtoms() const { return atoms; }

    std::vector<MCResidue>& getResidues() { return residues; }
    const std::vector<MCResidue>& getResidues() const { return residues; }

    TypeMaps& getResidueTypes() { return residueTypes; }
    const TypeMaps& getResidueTypes() const { return residueTypes; }

    TypeMaps& getAtomTypes() { return atomTypes; }
    const TypeMaps& getAtomTypes() const { return atomTypes; }

    MCInfo& getInfo() { return info; }
    const MCInfo& getInfo() const { return info; }

    MCForceField& getForceField() { return forcefield; }
    const MCForceField& getForceField() const { return forcefield; }

    // Movement residue management
    void addMovementResidue(const std::string& resName, int startIndex, int totalCount) {
        MCMovementResidueInfo info(startIndex, 0, totalCount, resName);
        movementResidues.push_back(info);
    }

    void updateMovementResidueCount(size_t index, int activeCount) {
        if (index < movementResidues.size()) {
            movementResidues[index].activeCount = activeCount;
        }
    }

    const std::vector<MCMovementResidueInfo>& getMovementResidues() const {
        return movementResidues;
    }

    // Atom management with swap-and-pop semantics
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
            // Swap with last active atom and decrease count
            if (index < activeAtomCount - 1) {
                atoms[index] = atoms[activeAtomCount - 1];
            }
            activeAtomCount--;
        }
    }

    // Residue management with swap-and-pop semantics
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
            // Mark as inactive
            residues[index].active = false;
            
            // Swap with last active residue and decrease count
            if (index < activeResidueCount - 1) {
                residues[index] = residues[activeResidueCount - 1];
                residues[index].active = true;
            }
            activeResidueCount--;
        }
    }

    // Energy utilities
    void updateResidueEnergy(int index, float vdw_energy, float elec_energy) {
        if (index >= 0 && index < activeResidueCount) {
            residues[index].energy_vdw = vdw_energy;
            residues[index].energy_elec = elec_energy;
        }
    }

    void updateEwaldEnergy(double real_space, double reciprocal, double self_energy) {
        ewald_energy.real_space = real_space;
        ewald_energy.reciprocal = reciprocal;
        ewald_energy.self = self_energy;
        ewald_energy.updateTotal();
    }

    // Statistics methods
    void incrementMoveStats(bool accepted) {
        info.stats.totalMoves++;
        if (accepted) {
            info.stats.acceptedMoves++;
        }
    }

    void incrementInsertionStats(bool accepted) {
        info.stats.insertionAttempts++;
        if (accepted) {
            info.stats.acceptedInsertions++;
        }
    }

    void incrementDeletionStats(bool accepted) {
        info.stats.deletionAttempts++;
        if (accepted) {
            info.stats.acceptedDeletions++;
        }
    }

    // System properties
    void setBoxDimensions(float x, float y, float z) {
        info.setBox(x, y, z);
    }

    void setTemperature(float temperature) {
        info.setTemperature(temperature);
    }

    void setCutoff(float cutoff) {
        info.cutoff = cutoff;
    }

    void setSwitchingFunction(bool use_switching, float r_on, float r_off) {
        info.use_switching = use_switching;
        info.r_on = r_on;
        info.r_off = r_off;
    }

    // Type management
    int getOrAddAtomType(const std::string& type) {
        return atomTypes.getOrAddType(type);
    }

    int getOrAddResidueType(const std::string& type) {
        return residueTypes.getOrAddType(type);
    }

    // Force field setup
    void setupForceField(int totalTypes, int movementTypes = 0) {
        forcefield.initialize(totalTypes, movementTypes);
        numMovementAtomTypes = movementTypes;
    }

    void setLJParameters(int type1, int type2, float sigma, float epsilon) {
        forcefield.setLJParams(type1, type2, sigma, epsilon);
    }

    // Validation and diagnostics
    bool checkConsistency() const {
        if (!is_valid()) return false;

        // Check active atom consistency
        for (int i = 0; i < activeResidueCount; ++i) {
            const auto& residue = residues[i];
            if (!residue.active) continue;
            
            // Check if residue's atoms are within active range
            if (residue.atomStart + residue.atomCount > activeAtomCount) {
                return false;
            }
        }

        return true;
    }

    std::string getSystemSummary() const {
        std::stringstream ss;
        ss << "MC System: " << activeAtomCount << "/" << atoms.size() << " atoms, "
           << activeResidueCount << "/" << residues.size() << " residues\n";
        ss << "Temperature: " << info.getTemperature() << " K, ";
        ss << "Box: [" << info.box[0] << ", " << info.box[1] << ", " << info.box[2] << "] nm\n";
        ss << "Acceptance rates: Overall=" << info.stats.getAcceptanceRate() 
           << ", Insertion=" << info.stats.getInsertionRate() 
           << ", Deletion=" << info.stats.getDeletionRate();
        return ss.str();
    }

    // Clone method
    std::unique_ptr<MCState> clone() const {
        auto cloned = std::make_unique<MCState>();
        
        // Copy all data
        cloned->atoms = atoms;
        cloned->residues = residues;
        cloned->residueTypes = residueTypes;
        cloned->atomTypes = atomTypes;
        cloned->movementResidues = movementResidues;
        cloned->movementAtomTypes = movementAtomTypes;
        cloned->activeAtomCount = activeAtomCount;
        cloned->activeResidueCount = activeResidueCount;
        cloned->numMovementAtomTypes = numMovementAtomTypes;
        cloned->info = info;
        cloned->forcefield = forcefield;
        cloned->ewald_energy = ewald_energy;
        
        return cloned;
    }

    // String representation
    std::string to_string() const {
        return getSystemSummary();
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

#endif // PYGCMC_MODEL_MONTECARLO_MAIN_HPP 