// modules/core/include/pygcmc/core/io/psf_parser.hpp

#ifndef PYGCMC_CORE_IO_PSF_PARSER_HPP
#define PYGCMC_CORE_IO_PSF_PARSER_HPP

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include "pygcmc/core/io/parser_common.hpp"

namespace pygcmc {
namespace core {
namespace io {

/**
 * @brief Parser for PSF files
 */
class PSFParser {
public:
    PSFParser() = default;
    ~PSFParser() = default;

    /**
     * @brief Parse a PSF file
     * @param filename Path to the PSF file
     * @return True if parsing was successful
     */
    bool parse(const std::string& filename);

    /**
     * @brief Get atom properties from PSF file
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
     * @brief Update PDB atoms with charge and mass from PSF
     * @param pdb_atoms Vector of PDB atoms to update
     * @return Number of atoms successfully updated
     */
    int update_pdb_atoms(std::vector<PDBAtom>& pdb_atoms) const;

    /**
     * @brief Get residues and their atoms that are missing topology information
     * @param atoms Vector of PDB atoms to check
     * @return Map of residue names to sets of atom names that are missing topology info
     */
    std::map<std::string, std::set<std::string>> get_missing_topology_info(
        const std::vector<PDBAtom>& atoms) const;

private:
    struct PSFAtom {
        std::string segment;
        std::string residue;
        std::string name;
        std::string type;
        int residue_number;
        double charge;
        double mass;
    };

    std::vector<PSFAtom> atoms_;
    // Index structure: residue -> residue_number -> atom_name -> index
    std::unordered_map<std::string,
        std::map<int,
            std::unordered_map<std::string, size_t>>> atom_index_;

    /**
     * @brief Parse the atoms section of the PSF file
     * @param lines Vector of lines from the atoms section
     * @return True if parsing was successful
     */
    bool parse_atoms_section(const std::vector<std::string>& lines);
};

} // namespace io
} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_IO_PSF_PARSER_HPP

