// modules/core/include/pygcmc/core/io/psf_parser.hpp

#pragma once

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <optional>
#include "pygcmc/core/io/parser_common.hpp"

namespace pygcmc {
namespace core {
namespace io {

enum class PSFParsingMode {
    Exact,   // Parse with exact residue numbers and all fields
    Rough    // Parse only essential fields like ITP parser
};

/**
 * @brief Parser for PSF files
 */
class PSFParser {
public:
    PSFParser() : is_first_file_(true) {}
    ~PSFParser() = default;

    /**
     * @brief Reset the parser state
     * This clears all stored atoms and resets the first file flag
     */
    void reset() {
        atoms_.clear();
        atom_index_.clear();
        id_index_.clear();
        is_first_file_ = true;
    }

    /**
     * @brief Parse a PSF file
     * @param filename Path to the PSF file
     * @param mode Parsing mode
     * @return True if parsing was successful
     */
    bool parse(const std::string& filename, PSFParsingMode mode = PSFParsingMode::Exact);

    /**
     * @brief Parse PSF file in rough mode (less strict parsing)
     * @param filename Path to PSF file
     * @return True if parsing was successful
     */
    bool parse_rough(const std::string& filename);

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

    /**
     * @brief Update PDB atoms with PSF topology information based on residue and atom order
     * @param atoms Vector of pointers to PDB atoms to update
     * @return Number of atoms updated
     */
    int update_pdb_atoms_by_order(std::vector<PDBAtom*>& atoms);

    /**
     * @brief Update PDB atoms with topology information from multiple PSF files
     * @param pdb_atoms Vector of pointers to PDB atoms to update
     * @param psf_files Vector of PSF file paths
     * @return Total number of atoms updated across all PSF files
     */
    static int update_pdb_atoms_from_multiple_psf(std::vector<PDBAtom*>& pdb_atoms, 
                                                const std::vector<std::string>& psf_files);

    // New functions for handling different PSF types
    static int update_pdb_atoms_multi_residue(std::vector<PDBAtom*>& pdb_atoms,
                                            const std::string& psf_file);
    
    static int update_pdb_atoms_single_residue(std::vector<PDBAtom*>& pdb_atoms,
                                             const std::string& psf_file,
                                             const std::string& target_residue);

private:
    bool is_first_file_;  // Track if we're parsing the first file
    struct PSFAtom {
        int id;                 ///< PSF atom ID
        std::string segment;    ///< Segment ID
        std::string residue;    ///< Residue name
        std::string name;       ///< Atom name
        std::string type;       ///< Atom type
        int residue_number;     ///< Residue number
        double charge;          ///< Atom charge
        double mass;            ///< Atom mass
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
     * @param mode Parsing mode
     * @return True if parsing was successful
     */
    bool parse_atoms_section(const std::vector<std::string>& lines, PSFParsingMode mode);

    /**
     * @brief Parse the atoms section of the PSF file
     * @param lines Vector of lines from the atoms section
     * @return True if parsing was successful
     */
    bool parse_atoms_section_rough(const std::vector<std::string>& lines);
};

} // namespace io
} // namespace core
} // namespace pygcmc

