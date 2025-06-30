// src/io/topology/topParserMain.hpp

#pragma once

#include "topParserStructures.hpp"
#include "topParserUtilities.hpp"
#include "topParserPreprocessor.hpp"
#include "topParserSections.hpp"
#include "model/ModelModule.hpp"
#include <string>
#include <vector>
#include <map>
#include <set>

namespace pygcmc {
namespace io {

/**
 * @brief Parser for GROMACS topology (.top) file format
 * 
 * This class handles reading topology files and populating a Topology object.
 * The topology file contains structural information about a molecular system,
 * including atoms, bonds, angles, dihedrals, improper dihedrals, and more.
 */
class TOPParser {
public:
    TOPParser() = default;
    ~TOPParser() = default;

    /**
     * @brief Enable or disable debug output
     */
    static void enable_debug(bool enable) {
        TopParserUtilities::getDebugFlag() = enable;
    }

    /**
     * @brief Check if debug output is enabled
     */
    static bool is_debug_enabled() {
        return TopParserUtilities::getDebugFlag();
    }

    /**
     * @brief Parse a topology file and populate a Topology object
     */
    bool parse_to_topology(const std::string& filename, model::Topology& topology);

    /**
     * @brief Static method to parse a topology file and return a new Topology object
     */
    static model::Topology parse_file(const std::string& filename);

    /**
     * @brief Static method to parse a topology string and return a new Topology object
     */
    static model::Topology parse_string(const std::string& top_str);

private:
    // Internal state
    std::set<std::string> processed_files_;
    std::string current_molecule_type_;
    int current_molecule_nrexcl_ = 3;
    std::vector<std::pair<std::string, int>> molecule_order_;
    std::map<std::string, std::string> molecule_to_segment_type_;
    
    // Store molecule definitions
    std::map<std::string, std::vector<LineInfo>> molecule_atoms_;
    std::map<std::string, std::vector<LineInfo>> molecule_bonds_;
    std::map<std::string, std::vector<LineInfo>> molecule_angles_;
    std::map<std::string, std::vector<LineInfo>> molecule_dihedrals_;
    std::map<std::string, std::vector<LineInfo>> molecule_impropers_;
};

} // namespace io
} // namespace pygcmc