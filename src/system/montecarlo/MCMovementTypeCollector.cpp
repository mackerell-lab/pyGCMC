#include "MCMovementTypeCollector.hpp"
#include "system/log/LogMain.hpp"
#include "system/common/SystemInterface.hpp"
#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace pygcmc {
namespace system {
namespace montecarlo {

std::vector<int> MCMovementTypeCollector::collectMovementTypes(
    const std::vector<MovementMolecularInfo>& molecules,
    model::TypeMaps& residueTypes,
    model::TypeMaps& atomTypes) {

    std::vector<int> newMovementAtomTypes;
    processedResidueNames_.clear();
    processedResidueNames_.reserve(molecules.size());

    // Add types from new molecules
    for (const auto& info : molecules) {
        if (!info.molecular) {
            throw std::runtime_error("Movement molecular data is null");
        }

        const auto& molRes = info.molecular->residues[0];
        std::string resName = processResidueName(molRes->get_resname());

        // Add this new residue name to the map
        residueTypes.getOrAddType(resName);
        processedResidueNames_.push_back(resName);

        // Add this new residue's atom types to the map and record their indices
        const auto& topology_atoms = info.molecular->topology_atoms;
        for (const auto& top_atom : topology_atoms) {
            int typeIdx = atomTypes.getOrAddType(top_atom.type);
            if (std::find(newMovementAtomTypes.begin(), newMovementAtomTypes.end(), typeIdx)
                == newMovementAtomTypes.end()) {
                newMovementAtomTypes.push_back(typeIdx);
            }
        }

        system::log::LogMain::log(system::common::LogLevel::DEBUG,
                         "Expected movement residue name for molecule: '", resName, "'");
        system::log::LogMain::log(system::common::LogLevel::DEBUG,
                         "Movement molecule residues:\n  resname='", resName, "'");
    }

    return newMovementAtomTypes;
}

std::string MCMovementTypeCollector::processResidueName(const std::string& rawName) const {
    std::string nameTrimmed = trim(rawName);
    std::string resName = nameTrimmed;
    std::transform(resName.begin(), resName.end(), resName.begin(),
                   [](unsigned char c){ return std::toupper(c); });
    return resName;
}

std::string MCMovementTypeCollector::trim(const std::string& s) const {
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
