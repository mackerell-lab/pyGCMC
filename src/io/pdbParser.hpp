// src/io/pdbParser.hpp

#ifndef PYGCMC_IO_PDB_PARSER_HPP
#define PYGCMC_IO_PDB_PARSER_HPP

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>

namespace pygcmc {
namespace data {
    class Atom;     // Forward declaration
    class Residue;  // Forward declaration
}

namespace io {

/**
 * @brief Parser for PDB format files and strings
 */
class PDBParser {
public:
    // Record types
    enum class RecordType {
        ATOM,
        HETATM,
        TER,
        HELIX,
        SHEET,
        SSBOND,
        UNKNOWN
    };

    // Parse results
    struct ParseResult {
        std::vector<std::shared_ptr<data::Atom>> atoms;
        std::vector<std::shared_ptr<data::Residue>> residues;
        std::unordered_map<std::string, std::vector<int>> helices;  // chain -> helix types
        std::unordered_map<std::string, std::vector<std::string>> sheets;  // chain -> sheet info
        std::vector<std::string> ssbonds;  // Disulfide bond information
    };

    /**
     * @brief Parse PDB file
     * @param filename Path to PDB file
     * @return ParseResult containing atoms, residues and structure information
     */
    static ParseResult parseFile(const std::string& filename);

    /**
     * @brief Parse PDB string
     * @param pdbStr String containing PDB data
     * @return ParseResult containing atoms, residues and structure information
     */
    static ParseResult parseString(const std::string& pdbStr);

private:
    static RecordType getRecordType(const std::string& line);
    
    static void parseAtomRecord(const std::string& line, RecordType type,
                              ParseResult& result,
                              std::shared_ptr<data::Residue>& currentResidue);
    
    static void parseTerRecord(const std::string& line,
                             std::shared_ptr<data::Residue>& currentResidue);
    
    static void parseHelixRecord(const std::string& line, ParseResult& result);
    
    static void parseSheetRecord(const std::string& line, ParseResult& result);
    
    static void parseSSBondRecord(const std::string& line, ParseResult& result);
};

} // namespace io
} // namespace pygcmc

#endif // PYGCMC_IO_PDB_PARSER_HPP


