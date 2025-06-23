#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_OPERATIONS_HPP
#define PYGCMC_MODEL_FORCEFIELD_OPERATIONS_HPP

#include "ForceFieldParams.hpp"
#include "ForceFieldInterface.hpp"
#include "ForceFieldOperationUtils.hpp"

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

    // === Existence Check Methods (using utility templates) ===

    bool has_atom_mass(const std::string& type) const {
        return has_parameter(storage_.atom_masses, type);
    }

    bool has_lj_params(const std::string& type) const {
        return has_parameter(storage_.lj_params, type);
    }

    bool has_nbfix(const std::string& type1, const std::string& type2) const {
        return has_parameter(storage_.nbfix, ParamKeyUtils::make_pair(type1, type2));
    }

    bool has_bond_params(const std::string& type1, const std::string& type2) const {
        return has_parameter(storage_.bond_params, ParamKeyUtils::make_pair(type1, type2));
    }

    bool has_angle_params(const std::string& type1, const std::string& type2,
                         const std::string& type3) const {
        return has_angle_bidirectional(storage_.angle_params, type1, type2, type3);
    }

    bool has_dihedral_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4) const {
        return has_parameter(storage_.dihedral_params, ParamKeyUtils::make_quad(type1, type2, type3, type4));
    }

    bool has_improper_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4) const {
        return has_parameter(storage_.improper_params, ParamKeyUtils::make_quad(type1, type2, type3, type4));
    }

    // === Size Query Methods (using utility templates) ===

    size_t get_num_atom_types() const { return get_container_size(storage_.atom_masses); }
    size_t get_num_lj_params() const { return get_container_size(storage_.lj_params); }
    size_t get_num_nbfix() const { return get_container_size(storage_.nbfix); }
    size_t get_num_bond_types() const { return get_container_size(storage_.bond_params); }
    size_t get_num_angle_types() const { return get_container_size(storage_.angle_params); }
    size_t get_num_dihedral_types() const { return get_container_size(storage_.dihedral_params); }
    size_t get_num_improper_types() const { return get_container_size(storage_.improper_params); }

    // === Summary Operations (using utility templates) ===

    bool has_basic_parameters() const {
        return !is_container_empty(storage_.atom_masses) && !is_container_empty(storage_.lj_params);
    }

    std::set<std::string> get_all_atom_types() const {
        return merge_type_sets(
            extract_single_types_all(storage_.atom_masses),
            extract_single_types_all(storage_.lj_params)
        );
    }

    std::set<std::string> get_bonded_atom_types() const {
        return merge_type_sets(
            extract_pair_types_all(storage_.bond_params),
            extract_triple_types_all(storage_.angle_params),
            extract_quad_types_all(storage_.dihedral_params),
            extract_quad_types_all(storage_.improper_params)
        );
    }

    std::set<std::string> get_nonbonded_atom_types() const {
        return merge_type_sets(
            extract_single_types_all(storage_.lj_params),
            extract_pair_types_all(storage_.nbfix)
        );
    }

private:
    const ForceFieldStorage& storage_;  ///< Reference to parameter storage
};

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_OPERATIONS_HPP 