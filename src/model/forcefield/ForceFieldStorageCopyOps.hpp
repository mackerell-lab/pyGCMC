#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_STORAGE_COPY_OPS_HPP
#define PYGCMC_MODEL_FORCEFIELD_STORAGE_COPY_OPS_HPP

#include "ForceFieldStorageTypes.hpp"
#include "ForceFieldParams.hpp"
#include <unordered_map>
#include <vector>
#include <memory>

namespace pygcmc {
namespace model {
namespace forcefield {

/**
 * @brief Copy operations template functions for ForceFieldStorageOptimized
 * 
 * These template functions provide efficient deep copying operations
 * for the complex data structures used in optimized storage.
 */

/**
 * @brief Deep copy a unique_ptr to hash map
 */
template<typename Key, typename Value, typename Hash = std::hash<Key>>
inline void deep_copy_hash_map(
    std::unique_ptr<std::unordered_map<Key, Value, Hash>>& dest,
    const std::unique_ptr<std::unordered_map<Key, Value, Hash>>& src) {
    
    if (src) {
        dest = std::make_unique<std::unordered_map<Key, Value, Hash>>(*src);
    } else {
        dest.reset();
    }
}

/**
 * @brief Copy all basic vector and map data
 */
template<typename Storage>
inline void copy_basic_storage_data(Storage& dest, const Storage& src) {
    dest.atom_types_ = src.atom_types_;
    dest.atom_type_to_index_ = src.atom_type_to_index_;
    dest.atom_masses_ = src.atom_masses_;
    dest.lj_params_ = src.lj_params_;
    dest.nonbonded_params_ = src.nonbonded_params_;
}

/**
 * @brief Deep copy all hash map unique_ptrs
 */
template<typename Storage>
inline void deep_copy_hash_maps(Storage& dest, const Storage& src) {
    deep_copy_hash_map(dest.nbfix_params_, src.nbfix_params_);
    deep_copy_hash_map(dest.bond_params_, src.bond_params_);
    deep_copy_hash_map(dest.angle_params_, src.angle_params_);
    deep_copy_hash_map(dest.dihedral_params_, src.dihedral_params_);
    deep_copy_hash_map(dest.improper_params_, src.improper_params_);
}

/**
 * @brief Reset all hash map unique_ptrs
 */
template<typename Storage>
inline void reset_all_hash_maps(Storage& storage) {
    storage.nbfix_params_.reset();
    storage.bond_params_.reset();
    storage.angle_params_.reset();
    storage.dihedral_params_.reset();
    storage.improper_params_.reset();
}

/**
 * @brief Complete copy operation combining basic data and hash maps
 */
template<typename Storage>
inline void complete_copy(Storage& dest, const Storage& src) {
    copy_basic_storage_data(dest, src);
    deep_copy_hash_maps(dest, src);
}

/**
 * @brief Complete assignment operation with reset and copy
 */
template<typename Storage>
inline void complete_assignment(Storage& dest, const Storage& src) {
    if (&dest != &src) {
        copy_basic_storage_data(dest, src);
        reset_all_hash_maps(dest);
        deep_copy_hash_maps(dest, src);
    }
}

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_STORAGE_COPY_OPS_HPP 