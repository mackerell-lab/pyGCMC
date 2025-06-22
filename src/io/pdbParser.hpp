// src/io/pdbParser.hpp

#pragma once

#ifndef PYGCMC_IO_PDBPARSER_HPP
#define PYGCMC_IO_PDBPARSER_HPP

#include <string>
#include <vector>
#include <memory>
#include <map>
#include <optional>
#include "model/ModelModule.hpp"

namespace pygcmc {
namespace io {

/**
 * @brief Parser for PDB format files and strings
 */
class PDBParser {
public:
    // Record types
    enum class RecordType {
        UNKNOWN,
        ATOM,
        HETATM,
        TER,
        HELIX,
        SHEET,
        SSBOND,
        CRYST1
    };

    /**
     * @brief Parse PDB file
     * @param filename Path to PDB file
     * @return Structure containing molecular structure information
     * @throws std::runtime_error if parsing fails
     */
    static model::Structure parse_file(const std::string& filename);

    /**
     * @brief Parse PDB string
     * @param pdbStr String containing PDB data
     * @return Structure containing molecular structure information
     * @throws std::runtime_error if parsing fails
     */
    static model::Structure parse_string(const std::string& pdbStr);

    /**
     * @brief Parse PDB file and populate a Structure object
     * @param filename Path to PDB file
     * @param structure Structure object to populate
     * @return true if parsing was successful, false otherwise
     */
    static bool parse_to_structure(const std::string& filename, model::Structure& structure);

    /**
     * @brief Parse PDB string and populate a Structure object
     * @param pdbStr String containing PDB data
     * @param structure Structure object to populate
     * @return true if parsing was successful, false otherwise
     */
    static bool parse_string_to_structure(const std::string& pdbStr, model::Structure& structure);

private:
    static RecordType getRecordType(const std::string& line);
    
    static bool parseAtomRecord(const std::string& line, RecordType type,
                              model::Structure& structure,
                              std::shared_ptr<model::Residue>& currentResidue);
    
    static bool parseTerRecord(const std::string& line,
                             std::shared_ptr<model::Residue>& currentResidue,
                             model::Structure& structure);
    
    static bool parseHelixRecord(const std::string& line, model::Structure& structure);
    
    static bool parseSheetRecord(const std::string& line, model::Structure& structure);
    
    static bool parseSSBondRecord(const std::string& line, model::Structure& structure);

    static bool parseCryst1Record(const std::string& line, model::Structure& structure);

    // Element mass table
    static const std::map<std::string, double> ELEMENT_MASSES;
};

} // namespace io
} // namespace pygcmc

#endif // PYGCMC_IO_PDBPARSER_HPP


