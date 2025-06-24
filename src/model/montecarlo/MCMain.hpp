#pragma once

#ifndef PYGCMC_MODEL_MONTECARLO_MAIN_HPP
#define PYGCMC_MODEL_MONTECARLO_MAIN_HPP

#include "MCStateCore.hpp"
#include "MCOperations.hpp"
#include "MCQueries.hpp"
#include "../common/ModelUtils.hpp"
#include <string>
#include <sstream>
#include <memory>

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

    // Movement residue management - delegate to MCOperations
    void addMovementResidue(const std::string& resName, int startIndex, int totalCount) {
        MCOperations::addMovementResidue(movementResidues, resName, startIndex, totalCount);
    }

    void updateMovementResidueCount(size_t index, int activeCount) {
        MCOperations::updateMovementResidueCount(movementResidues, index, activeCount);
    }

    const std::vector<MCMovementResidueInfo>& getMovementResidues() const {
        return movementResidues;
    }

    // Atom management - delegate to MCOperations
    int addAtom(const MCAtom& atom) {
        return MCOperations::addAtom(atoms, activeAtomCount, atom);
    }

    void removeAtom(int index) {
        MCOperations::removeAtom(atoms, activeAtomCount, index);
    }

    void updateAtomPosition(int index, float x, float y, float z) {
        if (index >= 0 && index < activeAtomCount) {
            MCOperations::updateAtomPosition(atoms[index], x, y, z);
        }
    }

    // Residue management - delegate to MCOperations
    int addResidue(const MCResidue& residue) {
        return MCOperations::addResidue(residues, activeResidueCount, residue);
    }

    void removeResidue(int index) {
        MCOperations::removeResidue(residues, activeResidueCount, index);
    }

    // Energy utilities - delegate to MCOperations
    void updateResidueEnergy(int index, float vdw_energy, float elec_energy) {
        if (index >= 0 && index < activeResidueCount) {
            MCOperations::updateResidueEnergy(residues[index], vdw_energy, elec_energy);
        }
    }

    void updateEwaldEnergy(double real_space, double reciprocal, double self_energy) {
        MCOperations::updateEwaldEnergy(ewald_energy, real_space, reciprocal, self_energy);
    }

    // Statistics methods - delegate to MCOperations
    void incrementMoveStats(bool accepted) {
        MCOperations::updateStatistics(info.stats, accepted);
    }

    void incrementInsertionStats(bool accepted) {
        MCOperations::updateStatistics(info.stats, accepted, true, false);
    }

    void incrementDeletionStats(bool accepted) {
        MCOperations::updateStatistics(info.stats, accepted, false, true);
    }

    // System properties - delegate to MCOperations
    void setBoxDimensions(float x, float y, float z) {
        MCOperations::setBox(info, x, y, z);
    }

    void setTemperature(float temperature) {
        MCOperations::setTemperature(info, temperature);
    }

    void setCutoff(float cutoff) {
        info.cutoff = cutoff;
    }

    void setSwitchingFunction(bool use_switching, float r_on, float r_off) {
        info.use_switching = use_switching;
        info.r_on = r_on;
        info.r_off = r_off;
    }

    // Type management - delegate to MCOperations
    int getOrAddAtomType(const std::string& type) {
        return MCOperations::getOrAddType(atomTypes, type);
    }

    int getOrAddResidueType(const std::string& type) {
        return MCOperations::getOrAddType(residueTypes, type);
    }

    // Force field setup - delegate to MCOperations
    void setupForceField(int totalTypes, int movementTypes = 0) {
        MCOperations::initializeForceField(forcefield, totalTypes, movementTypes);
        numMovementAtomTypes = movementTypes;
    }

    void setLJParameters(int type1, int type2, float sigma, float epsilon) {
        MCOperations::setLJParams(forcefield, type1, type2, sigma, epsilon);
    }

    // Validation and diagnostics - delegate to MCQueries
    bool checkConsistency() const {
        return MCQueries::isConsistentState(atoms, residues, activeAtomCount, activeResidueCount);
    }

    std::string getSystemSummary() const {
        std::stringstream ss;
        ss << "MC System: " << activeAtomCount << "/" << atoms.size() << " atoms, "
           << activeResidueCount << "/" << residues.size() << " residues\n";
        ss << "Temperature: " << MCQueries::getTemperature(info) << " K, ";
        ss << "Box: [" << info.box[0] << ", " << info.box[1] << ", " << info.box[2] << "] nm\n";
        ss << "Acceptance rates: Overall=" << MCQueries::getAcceptanceRate(info.stats) 
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

    // String representation - delegate to MCQueries
    std::string to_string() const {
        return getSystemSummary();
    }

    // Additional query methods - delegate to MCQueries
    std::optional<int> findAtomByType(int type) const {
        return MCQueries::findAtomByType(atoms, activeAtomCount, type);
    }

    std::optional<int> findResidueByType(int type) const {
        return MCQueries::findResidueByType(residues, activeResidueCount, type);
    }

    int countAtomsByType(int type) const {
        return MCQueries::countAtomsByType(atoms, activeAtomCount, type);
    }

    int countResiduesByType(int type) const {
        return MCQueries::countResiduesByType(residues, activeResidueCount, type);
    }

    double getTotalSystemEnergy() const {
        return ewald_energy.total + MCQueries::getTotalResidueEnergy(residues, activeResidueCount);
    }

    // System management - delegate to MCOperations
    void reserveCapacity(int max_atoms, int max_residues) {
        MCOperations::reserveCapacity(atoms, residues, info, max_atoms, max_residues);
    }

    void clearSystem() {
        MCOperations::clearSystem(atoms, residues, residueTypes, atomTypes,
                                movementResidues, movementAtomTypes,
                                activeAtomCount, activeResidueCount, numMovementAtomTypes,
                                ewald_energy, info.stats);
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