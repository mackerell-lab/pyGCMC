// modules/core/include/pygcmc/core/io/pdb_parser.hpp

#ifndef PYGCMC_CORE_IO_PDB_PARSER_HPP
#define PYGCMC_CORE_IO_PDB_PARSER_HPP

#include <string>
#include <vector>
#include <utility>
#include <unordered_map>
#include <map>
#include <set>
#include <optional>
#include "pygcmc/core/io/parser_common.hpp"

namespace pygcmc {
namespace core {
namespace io {

/**
 * @brief PDB Parser
 * 
 * Parses PDB files to extract crystal parameters and atom information.
 * Crystal parameters are optional - if CRYST1 record is not found or invalid,
 * returns std::nullopt for the crystal parameters.
 */
class PDBParser {
public:
    /**
     * @brief Parse a PDB file
     * 
     * @param filename Path to the PDB file
     * @return std::pair<std::optional<std::vector<double>>, std::vector<IOResidue>> 
     *         A pair containing optional crystal parameters (a, b, c, alpha, beta, gamma)
     *         and a list of IOResidues. Crystal parameters will be std::nullopt if
     *         CRYST1 record is not found or invalid.
     */
    static std::pair<std::optional<std::vector<double>>, std::vector<IOResidue>> parse(const std::string& filename);

private:
    static bool parse_atom_line(const std::string& line, PDBAtom& atom);
    static bool parse_cryst1_line(const std::string& line, std::vector<double>& cell_params);
    static bool is_atom_line(const std::string& line);
    static bool is_cryst1_line(const std::string& line);
    static std::string derive_element_from_name(const std::string& name);
    static bool validate_atom(const PDBAtom& atom);
    static bool validate_chain_structure(const std::unordered_map<char, 
        std::map<std::string, std::set<std::pair<int, char>>>>& chain_residues);
    static bool validate_pdb_structure(const std::vector<IOResidue>& residues);
};

} // namespace io
} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_IO_PDB_PARSER_HPP
