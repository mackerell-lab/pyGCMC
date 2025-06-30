// src/io/topology/psfParserMain.hpp

#pragma once

#include "../../model/ModelModule.hpp"
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
     * @brief Parse a PSF file and return a new Topology object
     * 
     * @param filename Path to the PSF file
     * @return model::Topology The parsed topology
     * @throws std::runtime_error if parsing fails
     */
    static model::Topology parse_file(const std::string& filename);

    /**
     * @brief Parse a PSF string and return a new Topology object
     * 
     * @param psf_str String containing PSF data
     * @return model::Topology The parsed topology
     * @throws std::runtime_error if parsing fails
     */
    static model::Topology parse_string(const std::string& psf_str);

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
    // Helper functions
    static std::string trim(const std::string& str);
};

} // namespace io
} // namespace pygcmc 