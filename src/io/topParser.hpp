// src/io/topParser.hpp

#pragma once

#include "../model/topology.hpp"
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <filesystem>

namespace pygcmc {
namespace io {

/**
 * @brief Parser for GROMACS topology (.top) file format
 * 
 * This class handles reading topology files and populating a Topology object.
 * The topology file contains structural information about a molecular system,
 * including atoms, bonds, angles, dihedrals, improper dihedrals, and more.
 */
class TopParser {
public:
    TopParser() = default;
    ~TopParser() = default;

    /**
     * @brief Parse a topology file and populate a Topology object
     * 
     * @param filename Path to the topology file
     * @param topology Topology object to populate
     * @return true if parsing was successful
     * @return false if there was an error
     */
    bool parse_to_topology(const std::string& filename, model::Topology& topology);

private:
    // Basic topology sections
    bool parse_defaults_section(const std::vector<std::string>& lines, model::Topology& topology);
    bool parse_atomtypes_section(const std::vector<std::string>& lines, model::Topology& topology);
    bool parse_moleculetype_section(const std::vector<std::string>& lines, model::Topology& topology);
    bool parse_atoms_section(const std::vector<std::string>& lines, model::Topology& topology);
    bool parse_bonds_section(const std::vector<std::string>& lines, model::Topology& topology);
    bool parse_angles_section(const std::vector<std::string>& lines, model::Topology& topology);
    bool parse_dihedrals_section(const std::vector<std::string>& lines, model::Topology& topology);
    bool parse_impropers_section(const std::vector<std::string>& lines, model::Topology& topology);
    bool parse_pairs_section(const std::vector<std::string>& lines, model::Topology& topology);
    bool parse_exclusions_section(const std::vector<std::string>& lines, model::Topology& topology);
    bool parse_cmap_section(const std::vector<std::string>& lines, model::Topology& topology);
    bool parse_system_section(const std::vector<std::string>& lines, model::Topology& topology);
    bool parse_molecules_section(const std::vector<std::string>& lines, model::Topology& topology);

    // Include handling
    bool process_includes(const std::string& filename, model::Topology& topology);
    std::string resolve_include_path(const std::string& include_path, const std::string& parent_file);

    // Preprocessor handling
    bool handle_preprocessor_line(const std::string& line);
    bool evaluate_ifdef_condition(const std::string& condition);
    bool should_process_line() const;

    // Helper functions
    std::vector<std::string> read_section(const std::string& filename, const std::string& section_name);
    std::string trim(std::string& str);
    std::vector<std::string> split(const std::string& str);
    void report_error(const std::string& message, bool critical);
    bool validate_topology(const model::Topology& topology);
    bool check_molecule_consistency(const model::Topology& topology);

    // Internal state
    std::set<std::string> processed_files_;
    std::string current_molecule_type_;
    int current_molecule_nrexcl_ = 3;
    std::vector<std::pair<std::string, int>> molecule_order_;
    std::map<std::string, std::string> molecule_to_segment_type_;
    std::set<std::string> preprocessor_defines_;
    bool strict_mode_ = true;

    struct PreprocessorState {
        bool in_ifdef = false;
        bool in_else = false;
        bool ifdef_condition_met = true;
        int ifdef_depth = 0;
        std::string current_ifdef;
    } preproc_state_;
};

} // namespace io
} // namespace pygcmc
