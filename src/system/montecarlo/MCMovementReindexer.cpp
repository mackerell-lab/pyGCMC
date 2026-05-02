#include "MCMovementReindexer.hpp"
#include "system/log/LogMain.hpp"
#include "system/common/SystemInterface.hpp"
#include <algorithm>
#include <cctype>

namespace pygcmc {
namespace system {
namespace montecarlo {

std::pair<std::vector<model::MCResidue>, std::vector<std::vector<model::MCAtom>>>
MCMovementReindexer::reindexExistingResidues(const model::MCState& oldState,
                                           model::TypeMaps& newResidueTypes,
                                           model::TypeMaps& newAtomTypes) {

    std::vector<model::MCResidue> reindexedResidues;
    std::vector<std::vector<model::MCAtom>> reindexedAtoms;

    reindexedResidues.reserve(oldState.activeResidueCount);
    reindexedAtoms.reserve(oldState.activeResidueCount);

    // Debug output: first print the old system's residue types and residues
    system::log::LogMain::log(system::common::LogLevel::DEBUG, "\nOriginal (old) active Residues before reindex:");
    for (int i = 0; i < oldState.activeResidueCount; i++) {
        const auto& oldRes = oldState.residues[i];
        if (!oldRes.active) {
            system::log::LogMain::log(system::common::LogLevel::DEBUG, "  Residue ", i, " inactive, skipping.");
            continue;
        }
        std::string oldTypeName = oldState.residueTypes.getTypeName(oldRes.type);
        system::log::LogMain::log(system::common::LogLevel::DEBUG, "  Residue ", i,
                         " old type index ", oldRes.type,
                         " => type name '", oldTypeName, "'");
    }

    for (int i = 0; i < oldState.activeResidueCount; i++) {
        const auto& oldRes = oldState.residues[i];

        // Get the old residue name and process it
        std::string oldResName = oldState.residueTypes.getTypeName(oldRes.type);
        std::string processedResName = processResidueName(oldResName);

        // Construct a new Residue, with type changed to the index in newResidueTypes
        model::MCResidue newRes = oldRes;
        newRes.type = newResidueTypes.getOrAddType(processedResName);

        // Push it to the result array
        reindexedResidues.push_back(newRes);

        // Now collect its atoms and remap the atom types
        std::vector<model::MCAtom> theseAtoms;
        theseAtoms.reserve(newRes.atomCount);
        for (int j = 0; j < oldRes.atomCount; j++) {
            auto oldAtom = oldState.atoms[oldRes.atomStart + j];
            // Get the atom type string from the old system
            std::string oldAtomTypeName = oldState.atomTypes.getTypeName(oldAtom.type);
            // Find the corresponding new index in the new typeMaps
            int newAtomTypeIdx = newAtomTypes.getOrAddType(oldAtomTypeName);
            oldAtom.type = newAtomTypeIdx;
            theseAtoms.push_back(oldAtom);
        }
        reindexedAtoms.push_back(std::move(theseAtoms));
    }

    return std::make_pair(std::move(reindexedResidues), std::move(reindexedAtoms));
}

std::string MCMovementReindexer::trim(const std::string& s) const {
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

std::string MCMovementReindexer::processResidueName(const std::string& rawName) const {
    std::string nameTrimmed = trim(rawName);
    std::string resName = nameTrimmed;
    std::transform(resName.begin(), resName.end(), resName.begin(),
                   [](unsigned char c){ return std::toupper(c); });
    return resName;
}

} // namespace montecarlo
} // namespace system
} // namespace pygcmc
