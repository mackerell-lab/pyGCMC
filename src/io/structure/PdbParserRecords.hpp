#pragma once
#ifndef PYGCMC_IO_STRUCTURE_PDBPARSERRECORDS_HPP
#define PYGCMC_IO_STRUCTURE_PDBPARSERRECORDS_HPP

#include <string>
#include <memory>
#include "model/ModelModule.hpp"
#include "PdbParserStructures.hpp"

namespace pygcmc {
namespace io {
namespace structure {

/**
 * @brief PDB record parsing functions
 */
class PdbParserRecords {
public:
    /**
     * @brief Parse ATOM/HETATM record
     * @param line PDB line to parse
     * @param type Record type (ATOM or HETATM)
     * @param structure Structure to populate
     * @param currentResidue Current residue being processed
     * @return bool True if parsing succeeded
     */
    static bool parseAtomRecord(const std::string& line, 
                               PdbParserStructures::RecordType type,
                               model::Structure& structure,
                               std::shared_ptr<model::Residue>& currentResidue);

    /**
     * @brief Parse TER record
     * @param line PDB line to parse
     * @param currentResidue Current residue being processed
     * @param structure Structure to populate
     * @return bool True if parsing succeeded
     */
    static bool parseTerRecord(const std::string& line,
                              std::shared_ptr<model::Residue>& currentResidue,
                              model::Structure& structure);

    /**
     * @brief Parse HELIX record
     * @param line PDB line to parse
     * @param structure Structure to populate
     * @return bool True if parsing succeeded
     */
    static bool parseHelixRecord(const std::string& line, model::Structure& structure);

    /**
     * @brief Parse SHEET record
     * @param line PDB line to parse
     * @param structure Structure to populate
     * @return bool True if parsing succeeded
     */
    static bool parseSheetRecord(const std::string& line, model::Structure& structure);

    /**
     * @brief Parse SSBOND record
     * @param line PDB line to parse
     * @param structure Structure to populate
     * @return bool True if parsing succeeded
     */
    static bool parseSSBondRecord(const std::string& line, model::Structure& structure);

    /**
     * @brief Parse CRYST1 record
     * @param line PDB line to parse
     * @param structure Structure to populate
     * @return bool True if parsing succeeded
     */
    static bool parseCryst1Record(const std::string& line, model::Structure& structure);
};

} // namespace structure
} // namespace io
} // namespace pygcmc

#endif // PYGCMC_IO_STRUCTURE_PDBPARSERRECORDS_HPP