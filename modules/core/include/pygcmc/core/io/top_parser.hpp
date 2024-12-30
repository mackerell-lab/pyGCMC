// modules/core/include/pygcmc/core/io/top_parser.hpp

#ifndef PYGCMC_CORE_IO_TOP_PARSER_HPP
#define PYGCMC_CORE_IO_TOP_PARSER_HPP

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include "pygcmc/core/io/parser_common.hpp"

namespace pygcmc {
namespace core {
namespace io {

/**
 * @brief Parser for GROMACS topology files
 */
class TopParser {
public:
    /**
     * @brief Parse a GROMACS topology file
     * @param filename Path to the topology file
     * @return True if parsing was successful
     */
    bool parse(const std::string& filename);

    /**
     * @brief Get the charge and mass for a specific atom
     * @param residue_name Residue name
     * @param atom_name Atom name
     * @param charge Output parameter for charge
     * @param mass Output parameter for mass
     * @return True if the atom was found
     */
    bool get_atom_properties(const std::string& residue_name, 
                           const std::string& atom_name,
                           double& charge,
                           double& mass) const;

    /**
     * @brief Update PDB atoms with charge and mass from topology
     * @param pdb_atoms Vector of PDB atoms to update
     * @return Number of atoms successfully updated
     */
    int update_pdb_atoms(std::vector<PDBAtom>& pdb_atoms) const;

private:
    struct TopAtom {
        std::string residue;
        std::string name;
        std::string type;    ///< Atom type from topology
        int residue_number;  ///< Residue number from topology
        double charge;
        double mass;
    };

    std::vector<TopAtom> atoms_;
    // Index structure: residue -> residue_number -> atom_name -> index
    std::unordered_map<std::string, 
        std::map<int, 
            std::unordered_map<std::string, size_t>>> atom_index_;

    /**
     * @brief Parse the atoms section of the topology file
     * @param lines Vector of lines from the atoms section
     * @return True if parsing was successful
     */
    bool parse_atoms_section(const std::vector<std::string>& lines);
};

} // namespace io
} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_IO_TOP_PARSER_HPP

