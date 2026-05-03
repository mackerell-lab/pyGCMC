// src/io/topology/PsfParserSections.hpp

#pragma once

#include "PsfParserSectionsBasic.hpp"
#include "PsfParserSectionsConnectivity.hpp"
#include "../../model/ModelModule.hpp"
#include <string>
#include <vector>

namespace pygcmc {
namespace io {

class PSFParserSections : public PSFParserSectionsBasic {
public:
    /**
     * Parse impropers section from PSF lines
     */
    static bool parse_impropers_from_lines(const std::vector<std::string>& lines, model::Topology& topology);

    /**
     * Parse donors section from PSF lines
     */
    static bool parse_donors_from_lines(const std::vector<std::string>& lines, model::Topology& topology);

    /**
     * Parse acceptors section from PSF lines
     */
    static bool parse_acceptors_from_lines(const std::vector<std::string>& lines, model::Topology& topology);

    /**
     * Parse CMAP section from PSF lines
     */
    static bool parse_cmap_from_lines(const std::vector<std::string>& lines, model::Topology& topology);

    /**
     * Parse groups section from PSF lines
     */
    static bool parse_groups_from_lines(const std::vector<std::string>& lines, model::Topology& topology);
};

} // namespace io
} // namespace pygcmc
