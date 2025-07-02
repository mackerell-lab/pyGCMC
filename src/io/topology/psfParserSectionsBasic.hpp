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
     * Parse atoms section from PSF lines with format flags
     */
    static bool parse_atoms_from_lines(const std::vector<std::string>& lines, model::Topology& topology,
                                     bool is_extended_format, bool is_drude_format);
};

} // namespace io
} // namespace pygcmc