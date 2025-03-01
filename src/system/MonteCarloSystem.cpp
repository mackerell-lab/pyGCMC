#include "MonteCarloSystem.hpp"
#include "system/system.hpp"
#include <cmath>
#include <algorithm>
#include <cctype>
#include <set>
#include <sstream>

namespace pygcmc {
namespace system {

using System = pygcmc::system::System;
using LogLevel = pygcmc::system::LogLevel;

void MonteCarloSystem::addInitialResidues(const model::MCResidue* resVec, int resCount,
                                         const model::MCAtom* atomVec, int atomCount) {
    if (resCount > state.info.maxResidues || atomCount > state.info.maxAtoms) {
        throw std::runtime_error("Initial system exceeds max capacity");
    }

    std::copy(resVec, resVec + resCount, state.residues.begin());
    std::copy(atomVec, atomVec + atomCount, state.atoms.begin());

    state.activeResidueCount = resCount;
    state.activeAtomCount = atomCount;
}

int MonteCarloSystem::insertResidue(const model::MCResidue& res, const model::MCAtom* atoms) {
    if (state.activeResidueCount >= state.info.maxResidues ||
        state.activeAtomCount + res.atomCount > state.info.maxAtoms) {
        return -1;
    }

    int resIdx = state.activeResidueCount++;
    state.residues[resIdx] = res;
    state.residues[resIdx].atomStart = state.activeAtomCount;
    state.residues[resIdx].active = true;

    std::copy(atoms, atoms + res.atomCount, state.atoms.begin() + state.activeAtomCount);
    state.activeAtomCount += res.atomCount;

    return resIdx;
}

bool MonteCarloSystem::removeResidue(int resIdx) {
    if (resIdx < 0 || resIdx >= state.activeResidueCount || !state.residues[resIdx].active) {
        return false;
    }

    // Get residue info and mark it as inactive
    model::MCResidue& res = state.residues[resIdx];
    res.active = false;
    int atomStart = res.atomStart;
    int atomCount = res.atomCount;

    // Move atoms if needed
    if (atomStart != state.activeAtomCount - atomCount) {
        for (int i = 0; i < atomCount; ++i) {
            state.atoms[atomStart + i] = state.atoms[state.activeAtomCount - atomCount + i];
        }
    }
    state.activeAtomCount -= atomCount;

    // If this is not the last residue, move the last active one to this position
    if (resIdx != state.activeResidueCount - 1) {
        state.residues[resIdx] = state.residues[state.activeResidueCount - 1];
        // Update the moved residue's atom start position if needed
        if (atomStart != state.activeAtomCount) {
            state.residues[resIdx].atomStart = atomStart;
        }
    }
    
    state.activeResidueCount--;
    return true;
}

void MonteCarloSystem::translateResidue(int resIdx, float dx, float dy, float dz) {
    if (resIdx >= 0 && resIdx < state.activeResidueCount) {
        model::MCResidue& res = state.residues[resIdx];
        // dx, dy, dz are expected to be in nm, no conversion needed
        for (int i = 0; i < res.atomCount; ++i) {
            model::MCAtom& atom = state.atoms[res.atomStart + i];
            atom.x += dx;  // All coordinates are in nm
            atom.y += dy;
            atom.z += dz;
            applyPBC(atom.x, atom.y, atom.z);
        }
        updateGeometricCenter(res);
    }
}

float MonteCarloSystem::calcNonBondedEnergy(const model::MCResidue& res1, const model::MCResidue& res2) const {
    (void)res1;  // Suppress unused parameter warning
    (void)res2;  // Suppress unused parameter warning
    float energy = 0.0f;
    // TODO: Implement LJ + Coulomb with periodic boundary conditions
    return energy;
}

float MonteCarloSystem::calcTotalEnergy() const {
    float energy = 0.0f;
    // TODO: Implement total energy calculation
    return energy;
}

void MonteCarloSystem::applyPBC(float& x, float& y, float& z) const {
    // Box dimensions and coordinates are in nm, no conversion needed
    x -= state.info.box[0] * std::floor(x / state.info.box[0]);
    y -= state.info.box[1] * std::floor(y / state.info.box[1]);
    z -= state.info.box[2] * std::floor(z / state.info.box[2]);
}

float MonteCarloSystem::getMinImageDistSqr(float dx, float dy, float dz) const {
    // All distances are in nm, no conversion needed
    dx -= state.info.box[0] * std::round(dx / state.info.box[0]);
    dy -= state.info.box[1] * std::round(dy / state.info.box[1]);
    dz -= state.info.box[2] * std::round(dz / state.info.box[2]);
    return dx*dx + dy*dy + dz*dz;  // Returns square of distance in nm²
}

void MonteCarloSystem::updateGeometricCenter(model::MCResidue& res) {
    res.center[0] = res.center[1] = res.center[2] = 0.0f;
    
    // All coordinates are already in nm, no conversion needed
    for (int i = 0; i < res.atomCount; ++i) {
        const model::MCAtom& atom = state.atoms[res.atomStart + i];
        res.center[0] += atom.x;
        res.center[1] += atom.y;
        res.center[2] += atom.z;
    }
    
    if (res.atomCount > 0) {
        float invCount = 1.0f / res.atomCount;
        res.center[0] *= invCount;
        res.center[1] *= invCount;
        res.center[2] *= invCount;
    }
}

void MonteCarloSystem::setSwitchingFunction(bool enable, float r_on, float r_off) {
    // Parameter validation
    if (r_on >= r_off) {
        throw std::runtime_error("Invalid switching function parameters: r_on must be less than r_off");
    }
    if (r_on <= 0.0f || r_off <= 0.0f) {
        throw std::runtime_error("Invalid switching function parameters: radii must be positive");
    }

    // Set parameters in MCState
    state.info.use_switching = enable;
    state.info.r_on = r_on;
    state.info.r_off = r_off;
}

float MonteCarloSystem::calculateSwitchingFunction(float r) const {
    if (!state.info.use_switching || r <= state.info.r_on) {
        return 1.0f;  // No switching function applied when r <= r_on
    }
    if (r >= state.info.r_off) {
        return 0.0f;  // Energy is zero when r >= r_off
    }
    
    // Calculate CHARMM-style switching function
    // S(r) = [(r_off^2 - r^2)^2 * (r_off^2 + 2r^2 - 3r_on^2)] / (r_off^2 - r_on^2)^3
    float r2 = r * r;
    float ron2 = state.info.r_on * state.info.r_on;
    float roff2 = state.info.r_off * state.info.r_off;
    
    float numerator = (roff2 - r2) * (roff2 - r2) * (roff2 + 2.0f*r2 - 3.0f*ron2);
    float denominator = (roff2 - ron2) * (roff2 - ron2) * (roff2 - ron2);
    
    return numerator / denominator;
}

void MonteCarloSystem::initializeFromMolecular(const std::shared_ptr<model::Molecular>& molecular) {
    if (!molecular) {
        throw std::runtime_error("MolecularSystem has no molecular data");
    }
    
    // Store molecular system for later parameter validation
    this->molecular = molecular;
    
    // Unit conversion constant
    const float ANGSTROM_TO_NM = 0.1f;    // 1 Å = 0.1 nm
    [[maybe_unused]] const float KCAL_TO_KJ = 4.184f;    // 1 kcal/mol = 4.184 kJ/mol
    
    // Set box dimensions from molecular system
    state.info.box[0] = molecular->boxDimensions[0] * ANGSTROM_TO_NM;  // Convert Å to nm
    state.info.box[1] = molecular->boxDimensions[1] * ANGSTROM_TO_NM;  // Convert Å to nm
    state.info.box[2] = molecular->boxDimensions[2] * ANGSTROM_TO_NM;  // Convert Å to nm
    state.info.volume = state.info.box[0] * state.info.box[1] * state.info.box[2];

    // Convert residues and atoms
    std::vector<model::MCResidue> tempResidues;
    std::vector<model::MCAtom> tempAtoms;
    
    size_t atomStart = 0;
    const size_t numResidues = molecular->get_num_residues();
    
    for (size_t i = 0; i < numResidues; ++i) {
        const auto& molRes = molecular->residues[i];
        const auto& topRes = molecular->topology_residues[i];
        
        model::MCResidue mcRes;
        mcRes.atomStart = atomStart;
        mcRes.atomCount = molRes->atom_count();
        mcRes.active = true;
        
        // Initialize energy components and GCMC parameters in GROMACS units
        mcRes.energy_vdw = 0.0f;   // kJ/mole
        mcRes.energy_elec = 0.0f;  // kJ/mole
        mcRes.chemPot = 0.0f;      // kJ/mole
        mcRes.concentration = 0.0f; // mol/L
        mcRes.radius = 0.0f;       // nm
        
        // Convert atoms for this residue
        const auto& molAtoms = molRes->get_atoms();
        for (size_t j = 0; j < molAtoms.size(); j++) {
            const auto& molAtom = molAtoms[j];
            const auto& topAtom = molecular->topology_atoms[topRes.atoms[j]];
            
            model::MCAtom mcAtom;
            // Convert coordinates from Å to nm
            mcAtom.x = molAtom->get_x() * ANGSTROM_TO_NM;
            mcAtom.y = molAtom->get_y() * ANGSTROM_TO_NM;
            mcAtom.z = molAtom->get_z() * ANGSTROM_TO_NM;
            mcAtom.charge = topAtom.charge;  // Charge unit (e) remains the same
            mcAtom.type = state.atomTypes.getOrAddType(topAtom.type);
            
            tempAtoms.push_back(mcAtom);
        }

        // Set residue type
        mcRes.type = state.residueTypes.getOrAddType(molRes->get_resname());

        // Calculate center of mass (in nm)
        mcRes.center[0] = mcRes.center[1] = mcRes.center[2] = 0.0f;
        for (const auto& atom : molAtoms) {
            mcRes.center[0] += atom->get_x() * ANGSTROM_TO_NM;
            mcRes.center[1] += atom->get_y() * ANGSTROM_TO_NM;
            mcRes.center[2] += atom->get_z() * ANGSTROM_TO_NM;
        }
        
        if (mcRes.atomCount > 0) {
            float invCount = 1.0f / mcRes.atomCount;
            mcRes.center[0] *= invCount;
            mcRes.center[1] *= invCount;
            mcRes.center[2] *= invCount;
        }
        
        tempResidues.push_back(mcRes);
        atomStart += mcRes.atomCount;
    }

    // Check capacity
    if (tempResidues.size() > static_cast<size_t>(state.info.maxResidues) ||
        tempAtoms.size() > static_cast<size_t>(state.info.maxAtoms)) {
        throw std::runtime_error("Initial system exceeds max capacity");
    }

    // Initialize the system with converted data
    addInitialResidues(tempResidues.data(), tempResidues.size(),
                      tempAtoms.data(), tempAtoms.size());
}

static inline std::string trim(const std::string &s) {
    auto start = s.begin();
    while (start != s.end() && std::isspace(*start)) {
        ++start;
    }
    auto end = s.end();
    while (end != start && std::isspace(*(end - 1))) {
        --end;
    }
    return std::string(start, end);
}

void MonteCarloSystem::addMovementMolecules(const std::vector<MovementMolecularInfo>& molecules)
{
    // --------------------------------------------------------------------
    // [1] First collect residue types and atom types from the "new molecules"
    // --------------------------------------------------------------------
    // Directly use the current state's type maps as a basis
    pygcmc::model::TypeMaps preResidueTypes = state.residueTypes;
    pygcmc::model::TypeMaps preAtomTypes = state.atomTypes;

    // Record atom types of the movement molecules
    std::vector<int> newMovementAtomTypes;
    int numNewMovementTypes = 0;

    // Add types from new molecules
    for (const auto& info : molecules) {
        if (!info.molecular) {
            throw std::runtime_error("Movement molecular data is null");
        }
        const auto& molRes = info.molecular->residues[0];
        std::string rawName = molRes->get_resname();
        std::string nameTrimmed = trim(rawName);
        std::string resName = nameTrimmed;
        std::transform(resName.begin(), resName.end(), resName.begin(),
                       [](unsigned char c){ return std::toupper(c); });

        // Add this new residue name to the map
        preResidueTypes.getOrAddType(resName);

        // Add this new residue's atom types to the map and record their indices
        const auto& topology_atoms = info.molecular->topology_atoms;
        for (const auto& top_atom : topology_atoms) {
            int typeIdx = preAtomTypes.getOrAddType(top_atom.type);
            if (std::find(newMovementAtomTypes.begin(), newMovementAtomTypes.end(), typeIdx) == newMovementAtomTypes.end()) {
                newMovementAtomTypes.push_back(typeIdx);
                numNewMovementTypes++;
            }
        }

        System::log(LogLevel::DEBUG, "Expected movement residue name for molecule: '", 
                   resName, "'");
        System::log(LogLevel::DEBUG, "Movement molecule residues:\n  resname='", 
                   resName, "'");
    }

    // --------------------------------------------------------------------
    // [3] Create a new MCState using the new type maps above
    // --------------------------------------------------------------------
    pygcmc::model::MCState newState;
    newState.info = state.info;
    newState.forcefield = state.forcefield;
    newState.residueTypes = preResidueTypes;
    newState.atomTypes = preAtomTypes;
    
    // Preserve previous movement residues information
    newState.movementResidues = state.movementResidues;
    newState.movementAtomTypes = state.movementAtomTypes;
    newState.numMovementAtomTypes = state.numMovementAtomTypes;

    // Add new movement atom types
    for (int typeIdx : newMovementAtomTypes) {
        if (std::find(newState.movementAtomTypes.begin(), newState.movementAtomTypes.end(), typeIdx) == newState.movementAtomTypes.end()) {
            newState.movementAtomTypes.push_back(typeIdx);
        }
    }
    newState.numMovementAtomTypes = newState.movementAtomTypes.size();

    // --------------------------------------------------------------------
    // [4] Reindex existing residues and atoms
    // --------------------------------------------------------------------
    std::vector<pygcmc::model::MCResidue> oldResidues;
    oldResidues.reserve(state.activeResidueCount);
    std::vector<std::vector<pygcmc::model::MCAtom>> oldResidueAtoms;
    oldResidueAtoms.reserve(state.activeResidueCount);

    // Debug output: first print the old system's residue types and residues
    System::log(LogLevel::DEBUG, "\nOriginal (old) active Residues before reindex:");
    for (int i = 0; i < state.activeResidueCount; i++) {
        const auto& oldRes = state.residues[i];
        if (!oldRes.active) {
            System::log(LogLevel::DEBUG, "  Residue ", i, " inactive, skipping.");
            continue;
        }
        std::string oldTypeName = state.residueTypes.getTypeName(oldRes.type);
        System::log(LogLevel::DEBUG, "  Residue ", i,
                   " old type index ", oldRes.type,
                   " => type name '", oldTypeName, "'");
    }

    for (int i = 0; i < state.activeResidueCount; i++) {
        const auto& oldRes = state.residues[i];
        // If you don't want to copy inactive residues, you can filter them. But your original code copied them all and processed them later
        // So here we don't differentiate active status.
        // Get the old name
        std::string oldResName = state.residueTypes.getTypeName(oldRes.type);
        std::string oldResNameTrimmed = trim(oldResName);
        std::string oldResNameUpper   = oldResNameTrimmed;
        std::transform(
            oldResNameUpper.begin(), oldResNameUpper.end(), oldResNameUpper.begin(),
            [](unsigned char c){ return std::toupper(c); }
        );

        // Construct a new Residue, with type changed to the index in newState.residueTypes
        pygcmc::model::MCResidue newRes = oldRes;
        newRes.type = newState.residueTypes.getOrAddType(oldResNameUpper);

        // Push it to the temporary array
        oldResidues.push_back(newRes);

        // Now collect its atoms and remap the atom types
        std::vector<pygcmc::model::MCAtom> theseAtoms;
        theseAtoms.reserve(newRes.atomCount);
        for (int j = 0; j < oldRes.atomCount; j++) {
            auto oldAtom = state.atoms[oldRes.atomStart + j];
            // Get the atom type string from the old system
            std::string oldAtomTypeName = state.atomTypes.getTypeName(oldAtom.type);
            // Find the corresponding new index in the new typeMaps
            int newAtomTypeIdx = newState.atomTypes.getOrAddType(oldAtomTypeName);
            oldAtom.type = newAtomTypeIdx;
            theseAtoms.push_back(oldAtom);
        }
        oldResidueAtoms.push_back(std::move(theseAtoms));
    }

    // --------------------------------------------------------------------
    // [5] Process new movement molecules
    // --------------------------------------------------------------------
    std::vector<std::string> insertionResNames;
    insertionResNames.reserve(molecules.size());

    for (const auto& info : molecules) {
        const auto& molRes = info.molecular->residues[0];
        std::string rawName = molRes->get_resname();
        std::string nameTrimmed = trim(rawName);
        std::string resName = nameTrimmed;
        std::transform(resName.begin(), resName.end(), resName.begin(),
                       [](unsigned char c){ return std::toupper(c); });
        insertionResNames.push_back(resName);
    }

    // For collection and grouping
    std::vector<std::vector<pygcmc::model::MCResidue>> matchingResidues(molecules.size());
    std::vector<std::vector<std::vector<pygcmc::model::MCAtom>>> matchingAtoms(molecules.size());
    std::vector<pygcmc::model::MCResidue> otherResidues;
    std::vector<std::vector<pygcmc::model::MCAtom>> otherAtoms;

    // First pass: iterate through oldResidues, match all residues (including inactive ones)
    for (int i = 0; i < static_cast<int>(oldResidues.size()); i++) {
        const auto& oldRes = oldResidues[i];
        std::string oldResNameUpper = newState.residueTypes.getTypeName(oldRes.type);

        System::log(LogLevel::DEBUG, "Processing residue ", i, 
                   ": upper='", oldResNameUpper, "'",
                   oldRes.active ? " (active)" : " (inactive)");

        bool matched = false;
        for (size_t m = 0; m < molecules.size(); m++) {
            System::log(LogLevel::DEBUG, "   Comparing with insertionResNames[", m, "]: '", 
                       insertionResNames[m], "'");
            if (oldResNameUpper == insertionResNames[m]) {
                System::log(LogLevel::DEBUG, "   Residue ", i, " matched movement molecule index ", m);
                matchingResidues[m].push_back(oldRes);
                std::vector<pygcmc::model::MCAtom> atoms;
                atoms.reserve(oldRes.atomCount);
                for (int j = 0; j < oldRes.atomCount; j++) {
                    atoms.push_back(oldResidueAtoms[i][j]);
                }
                matchingAtoms[m].push_back(atoms);
                matched = true;
                break;
            }
        }
        if (!matched) {
            System::log(LogLevel::DEBUG, "Residue ", i, " did not match any movement molecule; adding to others.");
            otherResidues.push_back(oldRes);
            std::vector<pygcmc::model::MCAtom> atoms;
            atoms.reserve(oldRes.atomCount);
            for (int j = 0; j < oldRes.atomCount; j++) {
                atoms.push_back(oldResidueAtoms[i][j]);
            }
            otherAtoms.push_back(std::move(atoms));
        }
    }

    // Output the number of residues matched for each movement molecule type
    for (size_t m = 0; m < molecules.size(); m++) {
        System::log(LogLevel::DEBUG, "Movement molecule index ", m, " ('", 
                   insertionResNames[m], "') collected ", 
                   matchingResidues[m].size(), " active residues.");
    }

    // --------------------------------------------------------------------
    // [6] Second pass: construct residues and atoms for newState
    // --------------------------------------------------------------------
    newState.atoms.reserve(state.atoms.size());
    newState.residues.reserve(state.residues.size());
    
    int newAtomStart = 0;
    int newResIdx = 0;

    // First put unmatched residues (otherResidues) into newState
    for (size_t i = 0; i < otherResidues.size(); i++) {
        pygcmc::model::MCResidue newRes = otherResidues[i];
        newRes.atomStart = newAtomStart;
        const auto& atoms = otherAtoms[i];
        for (const auto& atom : atoms) {
            newState.atoms.push_back(atom);
        }
        newAtomStart += static_cast<int>(atoms.size());
        newState.residues.push_back(newRes);
        newResIdx++;
    }

    // Then add old residues matched by each molecule and new inactive copies
    for (size_t m = 0; m < molecules.size(); m++) {
        const auto& molInfo = molecules[m];
        const auto& matches = matchingResidues[m];
        const auto& matchAtoms = matchingAtoms[m];

        // Record the starting position for this molecule type
        int startIndexForThisGroup = newResIdx;
        int activeCountForThisGroup = static_cast<int>(matches.size());

        // Add active residues matched from the "old system"
        for (size_t r = 0; r < matches.size(); r++) {
            pygcmc::model::MCResidue newRes = matches[r];
            newRes.atomStart = newAtomStart;
            newRes.active = true;
            newRes.fixed  = false;
            const auto& atoms = matchAtoms[r];
            for (const auto& atom : atoms) {
                newState.atoms.push_back(atom);
            }
            newAtomStart += static_cast<int>(atoms.size());
            newState.residues.push_back(newRes);
            newResIdx++;
        }

        // Check if this type already exists in movement residues
        bool typeExists = false;
        std::string resName = insertionResNames[m];
        int existingIndex = -1;
        for (size_t i = 0; i < newState.movementResidues.size(); i++) {
            if (newState.movementResidues[i].resName == resName) {
                typeExists = true;
                existingIndex = i;
                break;
            }
        }

        if (typeExists) {
            // If the type already exists, we need to keep the original inactive copies
            const auto& existingInfo = newState.movementResidues[existingIndex];
            // Continue adding from the starting position of the original inactive copies
            startIndexForThisGroup = existingInfo.startIndex;
            activeCountForThisGroup = existingInfo.activeCount;
            // No need to add new movement residue info
        } else {
            // Add inactive copies
            const auto& molRes = molInfo.molecular->residues[0];
            const auto& molAtoms = molRes->get_atoms();
            const auto& topRes = molInfo.molecular->topology_residues[0];
            int atomsPerResidue = static_cast<int>(molAtoms.size());
            
            for (int c = 0; c < molInfo.maxCopies; c++) {
                pygcmc::model::MCResidue newRes;
                newRes.atomStart = newAtomStart;
                newRes.atomCount = atomsPerResidue;
                newRes.active = false;
                newRes.fixed  = false;
                // Residue type index for the new molecule
                std::string rawName = molRes->get_resname();
                std::string nameTrimmed = trim(rawName);
                std::string resName = nameTrimmed;
                std::transform(
                    resName.begin(), resName.end(), resName.begin(),
                    [](unsigned char c){ return std::toupper(c); }
                );
                newRes.type = newState.residueTypes.getOrAddType(resName);

                // New molecule's atom type indices
                for (size_t i = 0; i < molAtoms.size(); i++) {
                    const auto& molAtom = molAtoms[i];
                    const auto& topAtom = molInfo.molecular->topology_atoms[topRes.atoms[i]];
                    
                    pygcmc::model::MCAtom mcAtom;
                    mcAtom.x = molAtom->get_x();
                    mcAtom.y = molAtom->get_y();
                    mcAtom.z = molAtom->get_z();
                    mcAtom.charge = topAtom.charge;
                    mcAtom.type = newState.atomTypes.getOrAddType(topAtom.type);
                    newState.atoms.push_back(mcAtom);
                }
                newAtomStart += atomsPerResidue;
                newState.residues.push_back(newRes);
                newResIdx++;
            }

            // Add new movement residue info
            pygcmc::model::MCMovementResidueInfo moveInfo;
            moveInfo.startIndex = startIndexForThisGroup;
            moveInfo.activeCount = activeCountForThisGroup;
            moveInfo.totalCount = activeCountForThisGroup + molInfo.maxCopies;
            moveInfo.resName = resName;
            newState.movementResidues.push_back(moveInfo);
        }
    }

    // Update newState's active counts
    newState.activeResidueCount = newResIdx;
    newState.activeAtomCount = newAtomStart;

    // Check capacity
    if (newState.activeResidueCount > newState.info.maxResidues ||
        newState.activeAtomCount > newState.info.maxAtoms) {
        throw std::runtime_error("New state exceeds max capacity after movement insertion");
    }

    // --------------------------------------------------------------------
    // [7] Print final mapping and residue information
    // --------------------------------------------------------------------
    System::log(LogLevel::DEBUG, "\nFinal ResidueTypes mapping (newState):");
    for (size_t idx = 0; idx < newState.residueTypes.atomTypes.size(); idx++) {
        System::log(LogLevel::DEBUG, "  index=", idx, 
                   " name='", newState.residueTypes.atomTypes[idx], "'");
    }

    System::log(LogLevel::DEBUG, "\nFinal Residues (newState):");
    for (int i = 0; i < newState.activeResidueCount; i++) {
        const pygcmc::model::MCResidue& newRes = newState.residues[i];
        std::string typeName = newState.residueTypes.getTypeName(newRes.type);
        System::log(LogLevel::DEBUG, "  Residue ", i,
                   " has type index ", newRes.type,
                   " => type name '", typeName, "'",
                   (newRes.active ? " (ACTIVE)" : " (INACTIVE)"));
    }

    // Print movement residues information
    System::log(LogLevel::DEBUG, "\nMovement Residues Info (newState):");
    for (const auto& info : newState.movementResidues) {
        System::log(LogLevel::DEBUG, "  Movement group: name='", info.resName,
                   "' start=", info.startIndex,
                   " active=", info.activeCount,
                   " total=", info.totalCount);
    }

    // Update the class member typeMaps with the new atom types
    state = std::move(newState);
}

void MonteCarloSystem::validateParameters(const model::ForceField& ff, const std::shared_ptr<model::Molecular>& molecular) {
    if (!molecular) {
        throw std::runtime_error("Molecular system is null");
    }

    // Check all atom topology parameters
    std::set<std::string> missingTopoTypes;
    for (const auto& res : molecular->topology_residues) {
        for (int atomIdx : res.atoms) {
            if (atomIdx >= static_cast<int>(molecular->topology_atoms.size())) {
                throw std::runtime_error("Invalid topology atom index: " + std::to_string(atomIdx));
            }
        }
    }

    // Check all atom type LJ parameters
    std::set<std::string> missingLJTypes;
    for (const auto& atom : molecular->topology_atoms) {
        try {
            ff.get_lj_params(atom.type);
        } catch (const std::exception&) {
            missingLJTypes.insert(atom.type);
        }
    }

    // If there are missing parameters, generate detailed error information
    if (!missingTopoTypes.empty() || !missingLJTypes.empty()) {
        std::stringstream error;
        error << "Missing parameters detected:\n";
        
        if (!missingTopoTypes.empty()) {
            error << "Missing topology parameters for atom types:\n";
            for (const auto& type : missingTopoTypes) {
                error << "  - " << type << "\n";
            }
        }
        
        if (!missingLJTypes.empty()) {
            error << "Missing LJ parameters for atom types:\n";
            for (const auto& type : missingLJTypes) {
                error << "  - " << type << "\n";
            }
        }
        
        throw std::runtime_error(error.str());
    }

    System::log(LogLevel::DEBUG, "All topology and force field parameters validated successfully");
}

void MonteCarloSystem::initializeForceField(const model::ForceField& ff) {
    // First verify all parameters
    validateParameters(ff, molecular);  // Need to save molecular as member variable

    const auto& atomTypes = state.atomTypes;
    const int numTypes = atomTypes.atomTypes.size();
    const int numMovementTypes = state.numMovementAtomTypes;
    
    if (numTypes == 0) {
        throw std::runtime_error("No atom types found in the system");
    }
    
    // Initialize force field parameters
    state.forcefield.numTotalTypes = numTypes;
    state.forcefield.numMovementTypes = numMovementTypes;  // Keep for future optimization
    
    // Extend array size to numTotalTypes * numTotalTypes
    state.forcefield.ljSigma.resize(numTypes * numTypes);
    state.forcefield.ljEps.resize(numTypes * numTypes);

    // Unit conversion constants
    const float ANGSTROM_TO_NM = 0.1f;  // 1 Å = 0.1 nm
    [[maybe_unused]] const float KCAL_TO_KJ = 4.184f;    // 1 kcal/mol = 4.184 kJ/mol

    // Loop through all possible type pairs
    for (int i = 0; i < numTypes; ++i) {
        const std::string& type1 = atomTypes.atomTypes[i];
        
        for (int j = 0; j < numTypes; ++j) {
            const std::string& type2 = atomTypes.atomTypes[j];
            const int pairIdx = i * numTypes + j;  // 2D array index
            
            // First try to get NBFIX parameters
            auto [nbfix_params, has_nbfix] = ff.get_nbfix(type1, type2);
            
            try {
                if (has_nbfix) {
                    // Use NBFIX parameters directly
                    // Convert Rmin from Å to nm
                    const float sigma = static_cast<float>(nbfix_params.rmin / std::pow(2.0, 1.0/6.0)) * ANGSTROM_TO_NM;
                    
                    // Store parameters (eps in kJ/mol, sigma in nm)
                    state.forcefield.ljSigma[pairIdx] = sigma;
                    state.forcefield.ljEps[pairIdx] = static_cast<float>(nbfix_params.epsilon) * KCAL_TO_KJ;
                } else {
                    // Get LJ parameters for two types
                    const auto& lj1 = ff.get_lj_params(type1);
                    const auto& lj2 = ff.get_lj_params(type2);
                    
                    // Convert Rmin/2 from Å to nm, then to sigma
                    const float sigma1 = static_cast<float>(2.0 * lj1.rmin_half / std::pow(2.0, 1.0/6.0)) * ANGSTROM_TO_NM;
                    const float sigma2 = static_cast<float>(2.0 * lj2.rmin_half / std::pow(2.0, 1.0/6.0)) * ANGSTROM_TO_NM;
                    
                    // Use Lorentz-Berthelot combination rule
                    const float sigma_avg = 0.5f * (sigma1 + sigma2);
                    const float eps_avg = std::sqrt(lj1.epsilon * lj2.epsilon);
                    
                    // Store parameters (convert eps to kJ/mol)
                    state.forcefield.ljSigma[pairIdx] = sigma_avg;
                    state.forcefield.ljEps[pairIdx] = eps_avg * KCAL_TO_KJ;
                }
            } catch (const std::exception& e) {
                throw std::runtime_error("Missing LJ parameters for atom type pair '" + 
                                      type1 + "'-'" + type2 + "'");
            }
        }
    }
}

} // namespace system
} // namespace pygcmc 