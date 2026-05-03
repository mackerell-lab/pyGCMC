// src/io/topology/TopParserSections.hpp

#pragma once

#include "TopParserStructures.hpp"
#include "TopParserUtilities.hpp"

namespace pygcmc {
namespace io {

/**
 * @brief Section parsing functions for TOP file
 */
class TopParserSections {
public:
    // Basic topology sections
    static bool parse_moleculetype_section(const std::vector<LineInfo>& lines, model::Topology& topology,
                                          std::string& current_molecule_type, int& current_molecule_nrexcl);

    static bool parse_atoms_section(const std::vector<LineInfo>& lines, model::Topology& topology,
                                   const std::string& current_molecule_type);

    static bool parse_bonds_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset = 0);

    static bool parse_angles_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset = 0);

    static bool parse_dihedrals_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset = 0);

    static bool parse_impropers_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset = 0);

    static bool parse_molecules_section(const std::vector<LineInfo>& lines, model::Topology& topology,
                                       std::vector<std::pair<std::string, int>>& molecule_order);

    static bool parse_cmaps_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset = 0);
};

} // namespace io
} // namespace pygcmc
