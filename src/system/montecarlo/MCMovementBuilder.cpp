#include "MCMovementBuilder.hpp"
#include "system/log/LogMain.hpp"
#include "system/common/SystemInterface.hpp"
#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace pygcmc {
namespace system {
namespace montecarlo {

void MCMovementBuilder::buildFinalState(
    model::MCState& newState,
    const std::vector<MovementMolecularInfo>& molecules,
    const std::vector<std::string>& processedResidueNames,
    const std::vector<model::MCResidue>& reindexedResidues,
    const std::vector<std::vector<model::MCAtom>>& reindexedAtoms) {

    // Reserve space for efficiency
    newState.atoms.reserve(newState.atoms.capacity());
    newState.residues.reserve(newState.residues.capacity());

    int newAtomStart = 0;
    int newResIdx = 0;

    // Match residues to movement molecules
    auto [matchingResidues, matchingAtoms, otherResidues, otherAtoms] =
        matchResidues(reindexedResidues, reindexedAtoms, processedResidueNames, newState);

    // First put unmatched residues into newState
    addUnmatchedResidues(newState, otherResidues, otherAtoms, newAtomStart, newResIdx);

    // Then add movement molecule groups
    addMovementMoleculeGroups(newState, molecules, processedResidueNames,
                             matchingResidues, matchingAtoms, newAtomStart, newResIdx);

    // Update newState's active counts
    newState.activeResidueCount = newResIdx;
    newState.activeAtomCount = newAtomStart;

    // Check capacity
    if (newState.activeResidueCount > newState.info.maxResidues ||
        newState.activeAtomCount > newState.info.maxAtoms) {
        throw std::runtime_error("New state exceeds max capacity after movement insertion");
    }

    // Print final mapping and residue information
    system::log::LogMain::log(system::common::LogLevel::DEBUG, "\nFinal ResidueTypes mapping (newState):");
    for (size_t idx = 0; idx < newState.residueTypes.atomTypes.size(); idx++) {
        system::log::LogMain::log(system::common::LogLevel::DEBUG, "  index=", idx,
                         " name='", newState.residueTypes.atomTypes[idx], "'");
    }

    system::log::LogMain::log(system::common::LogLevel::DEBUG, "\nFinal Residues (newState):");
    for (int i = 0; i < newState.activeResidueCount; i++) {
        const model::MCResidue& newRes = newState.residues[i];
        std::string typeName = newState.residueTypes.getTypeName(newRes.type);
        system::log::LogMain::log(system::common::LogLevel::DEBUG, "  Residue ", i,
                         " has type index ", newRes.type,
                         " => type name '", typeName, "'",
                         (newRes.active ? " (ACTIVE)" : " (INACTIVE)"));
    }

    // Print movement residues information
    system::log::LogMain::log(system::common::LogLevel::DEBUG, "\nMovement Residues Info (newState):");
    for (const auto& info : newState.movementResidues) {
        system::log::LogMain::log(system::common::LogLevel::DEBUG, "  Movement group: name='", info.resName,
                         "' start=", info.startIndex,
                         " active=", info.activeCount,
                         " total=", info.totalCount);
    }
}

std::tuple<
    std::vector<std::vector<model::MCResidue>>,
    std::vector<std::vector<std::vector<model::MCAtom>>>,
    std::vector<model::MCResidue>,
    std::vector<std::vector<model::MCAtom>>
> MCMovementBuilder::matchResidues(
    const std::vector<model::MCResidue>& reindexedResidues,
    const std::vector<std::vector<model::MCAtom>>& reindexedAtoms,
    const std::vector<std::string>& processedResidueNames,
    const model::MCState& newState) const {

    std::vector<std::vector<model::MCResidue>> matchingResidues(processedResidueNames.size());
    std::vector<std::vector<std::vector<model::MCAtom>>> matchingAtoms(processedResidueNames.size());
    std::vector<model::MCResidue> otherResidues;
    std::vector<std::vector<model::MCAtom>> otherAtoms;

    // Match residues
    for (size_t i = 0; i < reindexedResidues.size(); i++) {
        const auto& oldRes = reindexedResidues[i];
        std::string oldResNameUpper = newState.residueTypes.getTypeName(oldRes.type);

        system::log::LogMain::log(system::common::LogLevel::DEBUG, "Processing residue ", i,
                         ": upper='", oldResNameUpper, "'",
                         oldRes.active ? " (active)" : " (inactive)");

        bool matched = false;
        for (size_t m = 0; m < processedResidueNames.size(); m++) {
            system::log::LogMain::log(system::common::LogLevel::DEBUG, "   Comparing with processedResidueNames[", m, "]: '",
                             processedResidueNames[m], "'");
            if (oldResNameUpper == processedResidueNames[m]) {
                system::log::LogMain::log(system::common::LogLevel::DEBUG, "   Residue ", i, " matched movement molecule index ", m);
                matchingResidues[m].push_back(oldRes);
                matchingAtoms[m].push_back(reindexedAtoms[i]);
                matched = true;
                break;
            }
        }
        if (!matched) {
            system::log::LogMain::log(system::common::LogLevel::DEBUG, "Residue ", i, " did not match any movement molecule; adding to others.");
            otherResidues.push_back(oldRes);
            otherAtoms.push_back(reindexedAtoms[i]);
        }
    }

    // Output the number of residues matched for each movement molecule type
    for (size_t m = 0; m < processedResidueNames.size(); m++) {
        system::log::LogMain::log(system::common::LogLevel::DEBUG, "Movement molecule index ", m, " ('",
                         processedResidueNames[m], "') collected ",
                         matchingResidues[m].size(), " active residues.");
    }

    return std::make_tuple(std::move(matchingResidues), std::move(matchingAtoms),
                          std::move(otherResidues), std::move(otherAtoms));
}

void MCMovementBuilder::addUnmatchedResidues(
    model::MCState& newState,
    const std::vector<model::MCResidue>& otherResidues,
    const std::vector<std::vector<model::MCAtom>>& otherAtoms,
    int& newAtomStart,
    int& newResIdx) const {

    for (size_t i = 0; i < otherResidues.size(); i++) {
        model::MCResidue newRes = otherResidues[i];
        newRes.atomStart = newAtomStart;
        const auto& atoms = otherAtoms[i];
        for (const auto& atom : atoms) {
            newState.atoms.push_back(atom);
        }
        newAtomStart += static_cast<int>(atoms.size());
        newState.residues.push_back(newRes);
        newResIdx++;
    }
}

void MCMovementBuilder::addMovementMoleculeGroups(
    model::MCState& newState,
    const std::vector<MovementMolecularInfo>& molecules,
    const std::vector<std::string>& processedResidueNames,
    const std::vector<std::vector<model::MCResidue>>& matchingResidues,
    const std::vector<std::vector<std::vector<model::MCAtom>>>& matchingAtoms,
    int& newAtomStart,
    int& newResIdx) {

    // Add each movement molecule group
    for (size_t m = 0; m < molecules.size(); m++) {
        const auto& molInfo = molecules[m];
        const auto& matches = matchingResidues[m];
        const auto& matchAtoms = matchingAtoms[m];
        const std::string& resName = processedResidueNames[m];

        // Record the starting position for this molecule type
        int startIndexForThisGroup = newResIdx;
        int activeCountForThisGroup = static_cast<int>(matches.size());

        // Add active residues matched from the "old system"
        for (size_t r = 0; r < matches.size(); r++) {
            model::MCResidue newRes = matches[r];
            newRes.atomStart = newAtomStart;
            newRes.active = true;
            newRes.fixed = false;
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
        int existingIndex = -1;
        for (size_t i = 0; i < newState.movementResidues.size(); i++) {
            if (newState.movementResidues[i].resName == resName) {
                typeExists = true;
                existingIndex = i;
                break;
            }
        }

        if (typeExists) {
            // If the type already exists, update the count
            auto& existingInfo = newState.movementResidues[existingIndex];
            startIndexForThisGroup = existingInfo.startIndex;
            activeCountForThisGroup += existingInfo.activeCount;
            existingInfo.activeCount = activeCountForThisGroup;
            existingInfo.totalCount = activeCountForThisGroup + molInfo.maxCopies;
        } else {
            // Add inactive copies
            const auto& molRes = molInfo.molecular->residues[0];
            const auto& molAtoms = molRes->get_atoms();
            const auto& topRes = molInfo.molecular->topology_residues[0];
            int atomsPerResidue = static_cast<int>(molAtoms.size());

            for (int c = 0; c < molInfo.maxCopies; c++) {
                model::MCResidue newRes;
                newRes.atomStart = newAtomStart;
                newRes.atomCount = atomsPerResidue;
                newRes.active = false;
                newRes.fixed = false;
                newRes.type = newState.residueTypes.getOrAddType(resName);

                // Initialize energy components and GCMC parameters
                newRes.energy_vdw = 0.0f;
                newRes.energy_elec = 0.0f;
                newRes.chemPot = 0.0f;
                newRes.concentration = 0.0f;
                newRes.radius = 0.0f;

                // Add atoms for this inactive copy
                for (size_t i = 0; i < molAtoms.size(); i++) {
                    const auto& molAtom = molAtoms[i];
                    const auto& topAtom = molInfo.molecular->topology_atoms[topRes.atoms[i]];

                    model::MCAtom mcAtom;
                    // Convert coordinates from Å to nm
                    mcAtom.x = molAtom->get_x() * 0.1f;  // ANGSTROM_TO_NM
                    mcAtom.y = molAtom->get_y() * 0.1f;
                    mcAtom.z = molAtom->get_z() * 0.1f;
                    mcAtom.charge = topAtom.charge;
                    mcAtom.type = newState.atomTypes.getOrAddType(topAtom.type);
                    newState.atoms.push_back(mcAtom);
                }

                // Calculate center of mass (in nm)
                newRes.center[0] = newRes.center[1] = newRes.center[2] = 0.0f;
                for (const auto& atom : molAtoms) {
                    newRes.center[0] += atom->get_x() * 0.1f;  // Convert Å to nm
                    newRes.center[1] += atom->get_y() * 0.1f;
                    newRes.center[2] += atom->get_z() * 0.1f;
                }

                if (atomsPerResidue > 0) {
                    float invCount = 1.0f / atomsPerResidue;
                    newRes.center[0] *= invCount;
                    newRes.center[1] *= invCount;
                    newRes.center[2] *= invCount;
                }

                newAtomStart += atomsPerResidue;
                newState.residues.push_back(newRes);
                newResIdx++;
            }

            // Add new movement residue info
            model::MCMovementResidueInfo moveInfo;
            moveInfo.startIndex = startIndexForThisGroup;
            moveInfo.activeCount = activeCountForThisGroup;
            moveInfo.totalCount = activeCountForThisGroup + molInfo.maxCopies;
            moveInfo.resName = resName;
            newState.movementResidues.push_back(moveInfo);

            system::log::LogMain::log(system::common::LogLevel::DEBUG,
                             "Added movement residue info: ", resName,
                             " start=", moveInfo.startIndex,
                             " active=", moveInfo.activeCount,
                             " total=", moveInfo.totalCount);
        }
    }
}

std::string MCMovementBuilder::trim(const std::string& s) const {
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

} // namespace montecarlo
} // namespace system
} // namespace pygcmc
