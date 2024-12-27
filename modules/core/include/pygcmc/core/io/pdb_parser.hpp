// modules/core/include/pygcmc/core/io/pdb_parser.hpp

#ifndef PYGCMC_CORE_IO_PDB_PARSER_HPP
#define PYGCMC_CORE_IO_PDB_PARSER_HPP

#include "parser_common.hpp"

namespace pygcmc {
namespace core {
namespace io {

/**
 * @brief PDB Parser
 * 
 * Parses PDB files to extract crystal parameters and atom information.
 */
class PDBParser {
public:
    /**
     * @brief Parse a PDB file
     * 
     * @param filename Path to the PDB file
     * @return std::pair<std::vector<double>, std::vector<PDBAtom>> 
     *         A pair containing crystal parameters and a list of atoms
     */
    static std::pair<std::vector<double>, std::vector<PDBAtom>> parse(const std::string& filename);

private:
    static bool parse_cryst1_line(const std::string& line, std::vector<double>& cell_params);
    static bool parse_atom_line(const std::string& line, PDBAtom& atom);
    static bool validate_pdb_structure(const std::vector<PDBAtom>& atoms);

    // Private helper functions
    static std::string derive_element_from_name(const std::string& name);
    static bool validate_atom(const PDBAtom& atom);
    static bool validate_chain_structure(const std::unordered_map<char, 
        std::map<std::string, std::set<std::pair<int, char>>>>& chain_residues);
};

} // namespace io
} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_IO_PDB_PARSER_HPP
