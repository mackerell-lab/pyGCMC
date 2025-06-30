// src/io/topology/psfParserSectionsBasic.hpp

#pragma once

#include "../../model/ModelModule.hpp"
#include <string>
#include <vector>

namespace pygcmc {
namespace io {

class PSFParserSectionsBasic {
public:
    /**
     * Parse title section from PSF lines
     */
    static bool parse_title_from_lines(const std::vector<std::string>& lines, size_t& current_line, model::Topology& topology);
    
    /**
     * Parse atoms section from PSF lines
     */
    static bool parse_atoms_from_lines(const std::vector<std::string>& lines, model::Topology& topology);
    
    /**
     * Parse bonds section from PSF lines
     */
    static bool parse_bonds_from_lines(const std::vector<std::string>& lines, model::Topology& topology);
    
    /**
     * Parse angles section from PSF lines
     */
    static bool parse_angles_from_lines(const std::vector<std::string>& lines, model::Topology& topology);
    
    /**
     * Parse dihedrals section from PSF lines
     */
    static bool parse_dihedrals_from_lines(const std::vector<std::string>& dihedral_lines,
                                         model::Topology& topology, 
                                         const std::string& section_name);
    
    /**
     * Parse dihedrals section from PSF lines (original 2-parameter version)
     */
    static bool parse_dihedrals_from_lines(const std::vector<std::string>& dihedral_lines,
                                         model::Topology& topology);
};

} // namespace io
} // namespace pygcmc