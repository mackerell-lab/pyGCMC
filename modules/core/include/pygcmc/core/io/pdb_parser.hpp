// modules/core/include/pygcmc/core/io/pdb_parser.hpp

#ifndef PYGCMC_CORE_IO_PDB_PARSER_HPP
#define PYGCMC_CORE_IO_PDB_PARSER_HPP

#include "parser_common.hpp"

namespace pygcmc {
namespace core {
namespace io {

class PDBParser {
public:
    static std::pair<std::vector<double>, std::vector<PDBAtom>> parse(const std::string& filename);

private:
    static bool parse_cryst1_line(const std::string& line, std::vector<double>& cell_params);
    static bool parse_atom_line(const std::string& line, PDBAtom& atom);
    static bool validate_pdb_structure(const std::vector<PDBAtom>& atoms);
};

} // namespace io
} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_IO_PDB_PARSER_HPP
