#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_OPERATIONS_HPP
#define PYGCMC_MODEL_FORCEFIELD_OPERATIONS_HPP

#include "ForceFieldParams.hpp"
#include "ForceFieldInterface.hpp"

namespace pygcmc {
namespace model {
namespace forcefield {

/**
 * @brief Force field operations for existence checks and size queries
 * 
 * This class provides operations to check parameter existence and
 * get size information without retrieving the actual parameters.
 * All operations are read-only and optimized for performance.
 */
class ForceFieldOperations {
public:
    /**
     * @brief Constructor
     * @param storage Reference to the force field storage container
     */
    explicit ForceFieldOperations(const ForceFieldStorage& storage) : storage_(storage) {}
    ~ForceFieldOperations() = default;

    // === Existence Check Methods ===

    /**
     * @brief Check if atom mass exists for given type
     * @param type Atom type name
     * @return true if mass is defined for this type
     */
    bool has_atom_mass(const std::string& type) const {
        return storage_.atom_masses.find(type) != storage_.atom_masses.end();
    }

    /**
     * @brief Check if LJ parameters exist for given type
     * @param type Atom type name
     * @return true if LJ parameters are defined for this type
     */
    bool has_lj_params(const std::string& type) const {
        return storage_.lj_params.find(type) != storage_.lj_params.end();
    }

    /**
     * @brief Check if NBFIX parameters exist for given pair
     * @param type1 First atom type
     * @param type2 Second atom type
     * @return true if NBFIX parameters are defined for this pair
     */
    bool has_nbfix(const std::string& type1, const std::string& type2) const {
        auto key = ParamKeyUtils::make_pair(type1, type2);
        return storage_.nbfix.find(key) != storage_.nbfix.end();
    }

    /**
     * @brief Check if bond parameters exist for given pair
     * @param type1 First atom type
     * @param type2 Second atom type
     * @return true if bond parameters are defined for this pair
     */
    bool has_bond_params(const std::string& type1, const std::string& type2) const {
        auto key = ParamKeyUtils::make_pair(type1, type2);
        return storage_.bond_params.find(key) != storage_.bond_params.end();
    }

    /**
     * @brief Check if angle parameters exist for given triple
     * @param type1 First atom type
     * @param type2 Second atom type (center)
     * @param type3 Third atom type
     * @return true if angle parameters are defined for this triple (either orientation)
     */
    bool has_angle_params(const std::string& type1, const std::string& type2,
                         const std::string& type3) const {
        auto key1 = std::make_tuple(type1, type2, type3);
        auto key2 = std::make_tuple(type3, type2, type1);
        return storage_.angle_params.find(key1) != storage_.angle_params.end() ||
               storage_.angle_params.find(key2) != storage_.angle_params.end();
    }

    /**
     * @brief Check if dihedral parameters exist for given quartet
     * @param type1 First atom type
     * @param type2 Second atom type
     * @param type3 Third atom type
     * @param type4 Fourth atom type
     * @return true if dihedral parameters are defined for this quartet
     */
    bool has_dihedral_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4) const {
        auto key = ParamKeyUtils::make_quad(type1, type2, type3, type4);
        return storage_.dihedral_params.find(key) != storage_.dihedral_params.end();
    }

    /**
     * @brief Check if improper parameters exist for given quartet
     * @param type1 First atom type
     * @param type2 Second atom type
     * @param type3 Third atom type
     * @param type4 Fourth atom type
     * @return true if improper parameters are defined for this quartet
     */
    bool has_improper_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4) const {
        auto key = ParamKeyUtils::make_quad(type1, type2, type3, type4);
        return storage_.improper_params.find(key) != storage_.improper_params.end();
    }

    // === Size Query Methods ===

    /**
     * @brief Get number of atom types with masses defined
     * @return Number of atom types
     */
    size_t get_num_atom_types() const { 
        return storage_.atom_masses.size(); 
    }

    /**
     * @brief Get number of atom types with LJ parameters defined
     * @return Number of LJ parameter sets
     */
    size_t get_num_lj_params() const { 
        return storage_.lj_params.size(); 
    }

    /**
     * @brief Get number of NBFIX parameter sets defined
     * @return Number of NBFIX pairs
     */
    size_t get_num_nbfix() const { 
        return storage_.nbfix.size(); 
    }

    /**
     * @brief Get number of bond parameter sets defined
     * @return Number of bond types
     */
    size_t get_num_bond_types() const { 
        return storage_.bond_params.size(); 
    }

    /**
     * @brief Get number of angle parameter sets defined
     * @return Number of angle types
     */
    size_t get_num_angle_types() const { 
        return storage_.angle_params.size(); 
    }

    /**
     * @brief Get number of dihedral parameter sets defined
     * @return Number of dihedral types
     */
    size_t get_num_dihedral_types() const { 
        return storage_.dihedral_params.size(); 
    }

    /**
     * @brief Get number of improper parameter sets defined
     * @return Number of improper types
     */
    size_t get_num_improper_types() const { 
        return storage_.improper_params.size(); 
    }

    // === Summary Operations ===

    /**
     * @brief Check if force field has basic required parameters
     * @return true if force field has at least some atom masses and LJ parameters
     */
    bool has_basic_parameters() const {
        return !storage_.atom_masses.empty() && !storage_.lj_params.empty();
    }

    /**
     * @brief Get all atom types that have any parameters defined
     * @return Set of all referenced atom types
     */
    std::set<std::string> get_all_atom_types() const {
        std::set<std::string> types;

        // From atom masses
        for (const auto& pair : storage_.atom_masses) {
            types.insert(pair.first);
        }

        // From LJ parameters
        for (const auto& pair : storage_.lj_params) {
            types.insert(pair.first);
        }

        return types;
    }

    /**
     * @brief Get all atom types referenced in bonded parameters
     * @return Set of atom types used in bonds, angles, dihedrals, impropers
     */
    std::set<std::string> get_bonded_atom_types() const {
        std::set<std::string> types;

        // From bond parameters
        for (const auto& pair : storage_.bond_params) {
            types.insert(pair.first.first);
            types.insert(pair.first.second);
        }

        // From angle parameters
        for (const auto& pair : storage_.angle_params) {
            types.insert(std::get<0>(pair.first));
            types.insert(std::get<1>(pair.first));
            types.insert(std::get<2>(pair.first));
        }

        // From dihedral parameters
        for (const auto& pair : storage_.dihedral_params) {
            types.insert(std::get<0>(pair.first));
            types.insert(std::get<1>(pair.first));
            types.insert(std::get<2>(pair.first));
            types.insert(std::get<3>(pair.first));
        }

        // From improper parameters
        for (const auto& pair : storage_.improper_params) {
            types.insert(std::get<0>(pair.first));
            types.insert(std::get<1>(pair.first));
            types.insert(std::get<2>(pair.first));
            types.insert(std::get<3>(pair.first));
        }

        return types;
    }

    /**
     * @brief Get all atom types referenced in nonbonded parameters (LJ, NBFIX)
     * @return Set of atom types used in nonbonded interactions
     */
    std::set<std::string> get_nonbonded_atom_types() const {
        std::set<std::string> types;

        // From LJ parameters
        for (const auto& pair : storage_.lj_params) {
            types.insert(pair.first);
        }

        // From NBFIX parameters
        for (const auto& pair : storage_.nbfix) {
            types.insert(pair.first.first);
            types.insert(pair.first.second);
        }

        return types;
    }

private:
    const ForceFieldStorage& storage_;  ///< Reference to parameter storage
};

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_OPERATIONS_HPP 