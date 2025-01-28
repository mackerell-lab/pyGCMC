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
    // Structure to track line source information
    struct LineInfo {
        std::string content;      // The actual line content
        std::string source_file;  // Source file path
        int line_number;          // Line number in source file
        
        LineInfo(const std::string& content, const std::string& file, int line) 
            : content(content), source_file(file), line_number(line) {}
    };

    // Preprocessor state
    struct PreprocessorState {
        std::map<std::string, std::string> defines;  // #define macros
        std::vector<bool> ifdef_stack;               // Stack for #ifdef/#ifndef nesting
        std::vector<bool> else_encountered;          // Track if #else was encountered at each nesting level
        bool skip_section = false;                   // Whether to skip current section due to #ifdef
        
        bool should_skip() const {
            // Skip if any level in the stack is false
            for (bool val : ifdef_stack) {
                if (!val) return true;
            }
            return false;
        }
    };

    // New helper functions for flattened include processing
    bool collect_all_lines(const std::string& filename, std::vector<LineInfo>& all_lines, 
                          PreprocessorState& pp_state, bool is_main_file = true);
    void parse_sections(const std::vector<LineInfo>& all_lines, 
                       std::map<std::string, std::vector<LineInfo>>& sections);
    bool process_preprocessor_line(const std::string& line, const std::string& parent_file,
                                 std::vector<LineInfo>& all_lines, PreprocessorState& pp_state,
                                 int line_number);

    // Basic topology sections
    bool parse_defaults_section(const std::vector<LineInfo>& lines, model::Topology& topology);
    bool parse_atomtypes_section(const std::vector<LineInfo>& lines, model::Topology& topology);
    bool parse_moleculetype_section(const std::vector<LineInfo>& lines, model::Topology& topology);
    bool parse_atoms_section(const std::vector<LineInfo>& lines, model::Topology& topology);
    bool parse_bonds_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset = 0);
    bool parse_angles_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset = 0);
    bool parse_dihedrals_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset = 0);
    bool parse_impropers_section(const std::vector<LineInfo>& lines, model::Topology& topology, int atom_offset = 0);
    bool parse_system_section(const std::vector<LineInfo>& lines, model::Topology& topology);
    bool parse_molecules_section(const std::vector<LineInfo>& lines, model::Topology& topology);

    // Include handling
    std::string resolve_include_path(const std::string& include_path, const std::string& parent_file);

    // Helper functions
    std::string trim(std::string& str);
    std::vector<std::string> split(const std::string& str);
    std::string remove_comment(const std::string& line);  // New helper for comment handling

    // Internal state
    std::set<std::string> processed_files_;
    std::string current_molecule_type_;
    int current_molecule_nrexcl_ = 3;
    std::vector<std::pair<std::string, int>> molecule_order_;
    std::map<std::string, std::string> molecule_to_segment_type_;
    
    // Store molecule definitions (now using LineInfo)
    std::map<std::string, std::vector<LineInfo>> molecule_atoms_;
    std::map<std::string, std::vector<LineInfo>> molecule_bonds_;
    std::map<std::string, std::vector<LineInfo>> molecule_angles_;
    std::map<std::string, std::vector<LineInfo>> molecule_dihedrals_;
    std::map<std::string, std::vector<LineInfo>> molecule_impropers_;
};

} // namespace io
} // namespace pygcmc
