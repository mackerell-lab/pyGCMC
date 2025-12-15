#pragma once
#ifndef PYGCMC_IO_TOPOLOGY_FRAGMENTLIBRARY_HPP
#define PYGCMC_IO_TOPOLOGY_FRAGMENTLIBRARY_HPP

#include "../../model/montecarlo/MCStructures.hpp"
#include <map>
#include <string>
#include <vector>

namespace pygcmc {
namespace io {
namespace topology {

class FragmentLibrary {
public:
    struct TemplateData {
        std::string name;
        int typeId = -1;
        std::vector<model::montecarlo::MCAtom> atoms;
        // Per-atom force field type names from the ITP "type" column (same order as atoms).
        // Used later to map atoms onto MCState atom type indices.
        std::vector<std::string> atomTypeNames;
        double radius = 0.0;
        double molecularWeight = 0.0;
    };

    void addTemplate(const TemplateData& t);
    bool loadFromITP(const std::string& path, const std::string& name, int typeId,
                     const std::string& coordinatePdbFile = "");
    bool loadFromDirectory(const std::string& dir); // optional
    const TemplateData* get(const std::string& name) const;
    const TemplateData* getByType(int typeId) const;

private:
    std::map<std::string, TemplateData> byName_;
    std::map<int, std::string> byType_;
};

} // namespace topology
} // namespace io
} // namespace pygcmc

#endif // PYGCMC_IO_TOPOLOGY_FRAGMENTLIBRARY_HPP
