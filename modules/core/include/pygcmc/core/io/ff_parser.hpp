// modules/core/include/pygcmc/core/io/ff_parser.hpp

#ifndef PYGCMC_CORE_IO_FF_PARSER_HPP
#define PYGCMC_CORE_IO_FF_PARSER_HPP

#include <string>
#include <unordered_map>
#include <utility>
#include <stdexcept>
#include <vector>
#include "pygcmc/core/io/parser_common.hpp"

namespace pygcmc {
namespace core {
namespace io {

class FFParser {
public:
    // Parse CHARMM parameter file
    bool parse(const std::string& filename);

    // Get nonbonded parameters
    const std::unordered_map<std::string, ForceFieldPair>& get_nonbonded_params() const {
        return nonbonded_params_;
    }

    // Get NBFIX parameters
    const std::unordered_map<
        std::pair<std::string, std::string>, 
        ForceFieldPair, 
        PairStringHash
    >& get_nbfix_params() const {
        return nbfix_params_;
    }

    /**
     * @brief Update PDB atoms with force field parameters
     * @param atoms Vector of PDB atoms to update
     * @return Number of atoms successfully updated
     */
    int update_pdb_atoms(std::vector<PDBAtom>& atoms) const;

    /**
     * @brief Update PDB atoms with force field parameters
     * @param atoms Vector of PDB atoms to update
     * @return Number of atoms successfully updated
     */
    int update_pdb_atoms(std::vector<PDBAtom*>& atoms) const;   

private:
    // Storage for parameters
    std::unordered_map<std::string, ForceFieldPair> nonbonded_params_;
    std::unordered_map<
        std::pair<std::string, std::string>, 
        ForceFieldPair, 
        PairStringHash
    > nbfix_params_;

    // Helper functions for parsing sections
    void parse_nonbonded_line(const std::string& line);
    void parse_nbfix_line(const std::string& line);

    // Nonbonded parameters
    double cutnb_ = 14.0;
    double ctofnb_ = 12.0;
    double ctonnb_ = 10.0;
    double eps_ = 1.0;
    double e14fac_ = 1.0;
    double wmin_ = 1.5;
};

} // namespace io
} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_IO_FF_PARSER_HPP

