#pragma once
#ifndef PYGCMC_IO_STRUCTURE_PDBPARSERSTRUCTURES_HPP
#define PYGCMC_IO_STRUCTURE_PDBPARSERSTRUCTURES_HPP

#include <string>
#include <map>

namespace pygcmc {
namespace io {
namespace structure {

/**
 * @brief Constants and data structures for PDB parsing
 */
class PdbParserStructures {
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
     * @brief Get record type from PDB line
     * @param line PDB line to analyze
     * @return RecordType Identified record type
     */
    static RecordType getRecordType(const std::string& line);

    /**
     * @brief Get element mass table
     * @return const std::map<std::string, double>& Reference to element masses
     */
    static const std::map<std::string, double>& getElementMasses();

private:
    // Element mass table
    static const std::map<std::string, double> ELEMENT_MASSES;
};

} // namespace structure
} // namespace io
} // namespace pygcmc

#endif // PYGCMC_IO_STRUCTURE_PDBPARSERSTRUCTURES_HPP
