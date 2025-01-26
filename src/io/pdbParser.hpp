// src/io/pdbParser.hpp

#ifndef PYGCMC_IO_PDBPARSER_HPP
#define PYGCMC_IO_PDBPARSER_HPP

#include <string>
#include <vector>
#include <memory>
#include <map>
#include <optional>
#include "model/atom.hpp"
#include "model/residue.hpp"

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

    // Terminal information for TER records
    struct TerminalInfo {
        char chainId;
        int resSeq;
        char iCode;
        std::string resName;
    };

    // Parse results
    struct ParseResult {
        std::vector<std::shared_ptr<model::Atom>> atoms;
        std::vector<std::shared_ptr<model::Residue>> residues;
        std::map<std::string, std::vector<int>> helices;  // Chain -> helix classes
        std::map<std::string, std::vector<std::string>> sheets;  // Chain -> sheet info
        std::vector<std::string> ssbonds;  // Disulfide bond info
        std::optional<std::vector<double>> boxDimensions;  // Unit cell parameters
        std::vector<TerminalInfo> terminals;  // TER record information
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
                       std::shared_ptr<model::Residue>& currentResidue);
    
    static void parseTerRecord(const std::string& line,
                       std::shared_ptr<model::Residue>& currentResidue,
                       ParseResult& result);
    
    static void parseHelixRecord(const std::string& line, ParseResult& result);
    
    static void parseSheetRecord(const std::string& line, ParseResult& result);
    
    static void parseSSBondRecord(const std::string& line, ParseResult& result);

    static void parseCryst1Record(const std::string& line, ParseResult& result);
};

} // namespace io
} // namespace pygcmc

#endif // PYGCMC_IO_PDBPARSER_HPP


