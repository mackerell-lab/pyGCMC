// src/io/topParser.hpp

#pragma once

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <filesystem>
#include "../model/topology.hpp"

namespace pygcmc {
namespace io {

/**
 * @brief Parser for GROMACS topology (.top) files
 */
class TopParser {
public:
    TopParser() = default;
    ~TopParser() = default;

    /**
     * @brief Parse a GROMACS topology file and create a Topology object
     * @param filename Path to the .top file
     * @return A populated Topology object
     */
    model::Topology parse(const std::string& filename);

private:
    // Internal data structures for parsing
    struct AtomEntry {
        int id;
        std::string type;
        int residue_number;
        std::string residue_name;
        std::string atom_name;
        int charge_group;
        double charge;
        double mass;
        // Add support for B-state parameters
        std::string typeB;
        double chargeB;
        double massB;
    };

    struct BondEntry {
        int atom1;
        int atom2;
        int function_type;  // Add function type
        double length;
        double force_constant;
    };

    struct AngleEntry {
        int atom1;
        int atom2;
        int atom3;
        int function_type;  // Add function type
        double angle;
        double force_constant;
        double ub_length;
        double ub_constant;
    };

    struct DihedralEntry {
        int atom1;
        int atom2;
        int atom3;
        int atom4;
        int function_type;  // Add function type
        int multiplicity;
        double angle;
        double force_constant;
        bool improper;
    };

    struct PairEntry {
        int atom1;
        int atom2;
        int function_type;  // Add function type
        double c6;  // LJ C6 parameter
        double c12; // LJ C12 parameter
    };

    // Internal parsing methods
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

    // Helper methods
    std::vector<std::string> read_section(const std::string& filename, const std::string& section_name);
    void trim(std::string& str);
    std::vector<std::string> split(const std::string& str);
    
    // Include file handling
    bool process_includes(const std::string& filename, model::Topology& topology);
    std::string resolve_include_path(const std::string& include_path, const std::string& parent_file);
    std::set<std::string> processed_files_; // Keep track of processed files to avoid circular includes

    // Current molecule type being processed
    std::string current_molecule_type_;
    int current_molecule_nrexcl_ = 3;
};

} // namespace io
} // namespace pygcmc
