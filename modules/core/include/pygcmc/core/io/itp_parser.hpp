// modules/core/include/pygcmc/core/io/itp_parser.hpp

#ifndef PYGCMC_CORE_IO_ITP_PARSER_HPP
#define PYGCMC_CORE_IO_ITP_PARSER_HPP

#include <string>
#include <vector>
#include <map>
#include <set>
#include "pygcmc/core/io/parser_common.hpp"

namespace pygcmc {
namespace core {
namespace io {

/**
 * @brief Parser for ITP (GROMACS include topology) files
 */
class ITPParser {
public:
    ITPParser() = default;
    ~ITPParser() = default;

    /**
     * @brief Parse an ITP file
     * @param filename Path to the ITP file
     * @return True if parsing was successful
     */
    bool parse(const std::string& filename);

    /**
     * @brief Get atom properties from ITP file
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
     * @brief Update PDB atoms with topology information from ITP
     * @param pdb_atoms Vector of PDB atoms to update
     * @return Number of atoms successfully updated
     */
    int update_pdb_atoms(std::vector<PDBAtom>& pdb_atoms) const;

    /**
     * @brief Update PDB atoms with topology information from ITP
     * @param pdb_atoms Vector of pointers to PDB atoms to update
     * @return Number of atoms successfully updated
     */
    int update_pdb_atoms(std::vector<PDBAtom*>& pdb_atoms) const;

    /**
     * @brief Get missing topology information for atoms
     * @param atoms Vector of PDB atoms to check
     * @return Map of residue names to sets of atom names missing topology info
     */
    std::map<std::string, std::set<std::string>> get_missing_topology_info(
        const std::vector<PDBAtom>& atoms) const;

private:
    struct ITPAtom {
        std::string name;
        std::string type;
        std::string resname;
        int resid;
        double charge;
        double mass;
    };

    struct ResidueKey {
        std::string resname;
        int resid;

        bool operator<(const ResidueKey& other) const {
            if (resname != other.resname) {
                return resname < other.resname;
            }
            return resid < other.resid;
        }
    };

    std::vector<ITPAtom> itp_atoms_;
    std::map<ResidueKey, std::map<std::string, size_t>> atom_index_;

    /**
     * @brief Parse the atoms section of the ITP file
     * @param lines Vector of lines from the atoms section
     * @return True if parsing was successful
     */
    bool parse_atoms_section(const std::vector<std::string>& lines);
};

} // namespace io
} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_IO_ITP_PARSER_HPP

