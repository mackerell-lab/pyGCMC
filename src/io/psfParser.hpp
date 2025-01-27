#pragma once

#include "../model/topology.hpp"
#include <string>
#include <vector>

namespace pygcmc {
namespace io {

/**
 * @brief Parser for CHARMM PSF (Protein Structure File) format
 * 
 * This class handles reading PSF files and populating a Topology object.
 * The PSF file contains structural information about a molecular system,
 * including atoms, bonds, angles, dihedrals, improper dihedrals, and more.
 */
class PSFParser {
public:
    PSFParser() = default;
    ~PSFParser() = default;

    /**
     * @brief Parse a PSF file and populate a Topology object
     * 
     * @param filename Path to the PSF file
     * @param topology Topology object to populate
     * @return true if parsing was successful
     * @return false if there was an error
     */
    bool parse_to_topology(const std::string& filename, model::Topology& topology);

private:
    // Basic topology sections
    bool parse_title_from_lines(const std::vector<std::string>& lines, size_t& current_line, model::Topology& topology);
    bool parse_atoms_from_lines(const std::vector<std::string>& lines, size_t& current_line, model::Topology& topology);
    bool parse_bonds_from_lines(const std::vector<std::string>& lines, size_t& current_line, model::Topology& topology);
    bool parse_angles_from_lines(const std::vector<std::string>& lines, size_t& current_line, model::Topology& topology);
    bool parse_dihedrals_from_lines(const std::vector<std::string>& lines, size_t& current_line, model::Topology& topology);
    bool parse_impropers_from_lines(const std::vector<std::string>& lines, size_t& current_line, model::Topology& topology);
    
    // Nonbonded sections
    bool parse_donors_from_lines(const std::vector<std::string>& lines, size_t& current_line, model::Topology& topology);
    bool parse_acceptors_from_lines(const std::vector<std::string>& lines, size_t& current_line, model::Topology& topology);
    bool parse_nonbonded_exclusions(std::ifstream& file, model::Topology& topology);
    
    // Additional sections
    bool parse_groups_from_lines(const std::vector<std::string>& lines, size_t& current_line, model::Topology& topology);
    bool parse_cmap_from_lines(const std::vector<std::string>& lines, size_t& current_line, model::Topology& topology);

    // Helper functions
    bool read_section_header(std::ifstream& file, const std::string& expected_header, int& count);
    std::vector<int> read_index_block(std::ifstream& file, int expected_count, int indices_per_item);
    static std::string trim(const std::string& str);
};

} // namespace io
} // namespace pygcmc 