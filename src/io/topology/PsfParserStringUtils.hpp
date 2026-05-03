// src/io/topology/PsfParserStringUtils.hpp

#pragma once

#include "../../model/ModelModule.hpp"
#include <string>

namespace pygcmc {
namespace io {

class PSFParserStringUtils {
public:
    /**
     * Parse a PSF string in memory
     * @param psf_str The PSF string content
     * @return Parsed topology
     */
    static model::Topology parse_string(const std::string& psf_str);

    /**
     * Trim whitespace from both ends of a string
     * @param str Input string
     * @return Trimmed string
     */
    static std::string trim(const std::string& str);
};

} // namespace io
} // namespace pygcmc
