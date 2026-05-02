// src/io/topology/psfParserSectionsConnectivity.hpp

#pragma once

#include "model/ModelModule.hpp"
#include <string>
#include <vector>

namespace pygcmc {
namespace io {

class PSFParserSectionsConnectivity {
public:
    // Parse bonds section
    static bool parse_bonds_from_lines(const std::vector<std::string>& lines,
                                     model::Topology& topology);

    // Parse angles section
    static bool parse_angles_from_lines(const std::vector<std::string>& lines,
                                      model::Topology& topology);

    // Parse dihedrals section
    static bool parse_dihedrals_from_lines(const std::vector<std::string>& dihedral_lines,
                                         model::Topology& topology,
                                         const std::string& section_name);

    // Overloaded version for backward compatibility
    static bool parse_dihedrals_from_lines(const std::vector<std::string>& dihedral_lines,
                                         model::Topology& topology);
};

} // namespace io
} // namespace pygcmc
