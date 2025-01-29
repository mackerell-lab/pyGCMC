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

    // Structure to hold helix information
    struct HelixInfo {
        std::string helixId;
        std::string initResName;
        char initChainId;
        int initSeqNum;
        char initICode;
        std::string endResName;
        char endChainId;
        int endSeqNum;
        char endICode;
        int helixClass;
    };

    // Structure to hold terminal record information
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
        std::vector<TerminalInfo> terminals;
        std::map<std::string, std::vector<HelixInfo>> helices;  // Map of chain ID to helix information
        std::map<std::string, std::vector<std::string>> sheets;
        std::vector<std::string> ssbonds;
        std::vector<double> boxDimensions;
    };

    /**
     * @brief Parse PDB file
     * @param filename Path to PDB file
     * @return ParseResult containing atoms, residues and structure information
     * @throws std::runtime_error if parsing fails
     */
    static ParseResult parse_file(const std::string& filename);

    /**
     * @brief Parse PDB string
     * @param pdbStr String containing PDB data
     * @return ParseResult containing atoms, residues and structure information
     * @throws std::runtime_error if parsing fails
     */
    static ParseResult parse_string(const std::string& pdbStr);

    /**
     * @brief Parse PDB file and populate a ParseResult object
     * @param filename Path to PDB file
     * @param result ParseResult object to populate
     * @return true if parsing was successful, false otherwise
     */
    static bool parse_to_result(const std::string& filename, ParseResult& result);

    /**
     * @brief Parse PDB string and populate a ParseResult object
     * @param pdbStr String containing PDB data
     * @param result ParseResult object to populate
     * @return true if parsing was successful, false otherwise
     */
    static bool parse_string_to_result(const std::string& pdbStr, ParseResult& result);

private:
    static RecordType getRecordType(const std::string& line);
    
    static bool parseAtomRecord(const std::string& line, RecordType type,
                              ParseResult& result,
                              std::shared_ptr<model::Residue>& currentResidue);
    
    static bool parseTerRecord(const std::string& line,
                             std::shared_ptr<model::Residue>& currentResidue,
                             ParseResult& result);
    
    static bool parseHelixRecord(const std::string& line, ParseResult& result);
    
    static bool parseSheetRecord(const std::string& line, ParseResult& result);
    
    static bool parseSSBondRecord(const std::string& line, ParseResult& result);

    static bool parseCryst1Record(const std::string& line, ParseResult& result);

    // Element mass table
    static const std::map<std::string, double> ELEMENT_MASSES;
};

} // namespace io
} // namespace pygcmc

#endif // PYGCMC_IO_PDBPARSER_HPP


