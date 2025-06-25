#pragma once
#ifndef PYGCMC_IO_STRUCTURE_PDBPARSERMAIN_HPP
#define PYGCMC_IO_STRUCTURE_PDBPARSERMAIN_HPP

#include <string>
#include "model/ModelModule.hpp"
#include "PdbParserStructures.hpp"
#include "PdbParserRecords.hpp"

namespace pygcmc {
namespace io {
namespace structure {

/**
 * @brief Main PDB parser implementation
 */
class PdbParserMain {
public:
    /**
     * @brief Parse PDB file
     * @param filename Path to PDB file
     * @return model::Structure Parsed structure
     * @throw std::runtime_error If parsing fails
     */
    static model::Structure parse_file(const std::string& filename);

    /**
     * @brief Parse PDB string
     * @param pdbStr String containing PDB data
     * @return model::Structure Parsed structure
     * @throw std::runtime_error If parsing fails
     */
    static model::Structure parse_string(const std::string& pdbStr);

    /**
     * @brief Parse PDB file and populate a Structure object
     * @param filename Path to PDB file
     * @param structure Structure object to populate
     * @return bool True if parsing was successful, false otherwise
     */
    static bool parse_to_structure(const std::string& filename, model::Structure& structure);

    /**
     * @brief Parse PDB string and populate a Structure object
     * @param pdbStr String containing PDB data
     * @param structure Structure object to populate
     * @return bool True if parsing was successful, false otherwise
     */
    static bool parse_string_to_structure(const std::string& pdbStr, model::Structure& structure);

private:
    /**
     * @brief Common parsing logic for both file and string inputs
     * @param input Input stream to parse
     * @param structure Structure object to populate
     * @return bool True if parsing was successful, false otherwise
     */
    static bool parse_common(std::istream& input, model::Structure& structure);
};

} // namespace structure
} // namespace io
} // namespace pygcmc

#endif // PYGCMC_IO_STRUCTURE_PDBPARSERMAIN_HPP