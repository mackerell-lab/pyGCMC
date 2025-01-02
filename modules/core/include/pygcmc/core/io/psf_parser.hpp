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
     * @brief Get atom properties from PSF file
     * @param residue_name Residue name
     * @param residue_number Residue number
     * @param atom_name Atom name
     * @param charge Output parameter for charge
     * @param mass Output parameter for mass
     * @return True if the atom was found
     */
    bool get_atom_properties(const std::string& residue_name,
                           int residue_number,
                           const std::string& atom_name,
                           double& charge,
                           double& mass) const;

    /**
     * @brief Update PDB atoms with topology information from PSF file
     * @param pdb_atoms Vector of PDB atoms to update
     * @return Number of atoms updated
     */
    int update_pdb_atoms(std::vector<PDBAtom>& pdb_atoms) const;

    /**
     * @brief Update PDB atoms with topology information from PSF file
     * @param pdb_atoms Vector of pointers to PDB atoms to update
     * @return Number of atoms updated
     */
    int update_pdb_atoms(std::vector<PDBAtom*>& pdb_atoms) const;

    /**
     * @brief Get missing topology information for a set of atoms
     * @param atoms Vector of atoms to check
     * @return Map of residue names to sets of atom names that are missing topology info
     */
    std::map<std::string, std::set<std::string>> get_missing_topology_info(
        const std::vector<PDBAtom>& atoms) const;

private:
    struct PSFAtom {
        int id;                 ///< PSF 中的全局原子编号
        std::string segment;    ///< 对应 SEGID
        std::string residue;    ///< 残基名 (RESNAME)
        std::string name;       ///< 原子名 (ATOMNAME)
        std::string type;       ///< 原子类型 (ATOMTYPE)
        int residue_number;     ///< 残基序号 (RESID)
        double charge;          ///< 原子电荷
        double mass;            ///< 原子质量
    };

    std::vector<PSFAtom> atoms_;
    // Index structure: residue -> residue_number -> atom_name -> index
    std::unordered_map<std::string,
        std::map<int,
            std::unordered_map<std::string, size_t>>> atom_index_;
    // Index structure: atom_id -> index
    std::unordered_map<int, size_t> id_index_;

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

