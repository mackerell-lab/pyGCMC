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

    // Implementation split into smaller parts for maintainability
    // This is a placeholder for the complex movement molecule group addition logic
    // TODO: Implement the full logic from the original addMovementMolecules
    
    (void)newState;
    (void)molecules;
    (void)processedResidueNames;
    (void)matchingResidues;
    (void)matchingAtoms;
    (void)newAtomStart;
    (void)newResIdx;
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