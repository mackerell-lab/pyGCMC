#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_OPERATION_UTILS_HPP
#define PYGCMC_MODEL_FORCEFIELD_OPERATION_UTILS_HPP

#include "ForceFieldParams.hpp"
#include <set>
#include <map>
#include <tuple>

namespace pygcmc {
namespace model {
namespace forcefield {

/**
 * @brief Template utilities for common force field operations
 * 
 * These template functions provide reusable patterns for existence checks,
 * type extraction, and container operations across different parameter types.
 */

/**
 * @brief Generic existence check for map-like containers
 */
template<typename Container, typename Key>
inline bool has_parameter(const Container& container, const Key& key) {
    return container.find(key) != container.end();
}

/**
 * @brief Extract first element from pair keys
 */
template<typename Container>
inline std::set<std::string> extract_pair_types_first(const Container& container) {
    std::set<std::string> types;
    for (const auto& pair : container) {
        types.insert(pair.first.first);
    }
    return types;
}

/**
 * @brief Extract second element from pair keys
 */
template<typename Container>
inline std::set<std::string> extract_pair_types_second(const Container& container) {
    std::set<std::string> types;
    for (const auto& pair : container) {
        types.insert(pair.first.second);
    }
    return types;
}

/**
 * @brief Extract all types from pair keys
 */
template<typename Container>
inline std::set<std::string> extract_pair_types_all(const Container& container) {
    std::set<std::string> types;
    for (const auto& pair : container) {
        types.insert(pair.first.first);
        types.insert(pair.first.second);
    }
    return types;
}

/**
 * @brief Extract all types from triple keys (3-tuples)
 */
template<typename Container>
inline std::set<std::string> extract_triple_types_all(const Container& container) {
    std::set<std::string> types;
    for (const auto& pair : container) {
        types.insert(std::get<0>(pair.first));
        types.insert(std::get<1>(pair.first));
        types.insert(std::get<2>(pair.first));
    }
    return types;
}

/**
 * @brief Extract all types from quad keys (4-tuples)
 */
template<typename Container>
inline std::set<std::string> extract_quad_types_all(const Container& container) {
    std::set<std::string> types;
    for (const auto& pair : container) {
        types.insert(std::get<0>(pair.first));
        types.insert(std::get<1>(pair.first));
        types.insert(std::get<2>(pair.first));
        types.insert(std::get<3>(pair.first));
    }
    return types;
}

/**
 * @brief Extract all first elements from map (for single-string keys)
 */
template<typename Container>
inline std::set<std::string> extract_single_types_all(const Container& container) {
    std::set<std::string> types;
    for (const auto& pair : container) {
        types.insert(pair.first);
    }
    return types;
}

/**
 * @brief Merge multiple sets of strings
 */
template<typename... Sets>
inline std::set<std::string> merge_type_sets(const Sets&... sets) {
    std::set<std::string> merged;
    auto merge_one = [&merged](const std::set<std::string>& s) {
        merged.insert(s.begin(), s.end());
    };
    (merge_one(sets), ...);
    return merged;
}

/**
 * @brief Check bidirectional existence for angle parameters
 * 
 * Angles can be defined in either direction (A-B-C or C-B-A)
 */
template<typename Container>
inline bool has_angle_bidirectional(const Container& container,
                                   const std::string& type1,
                                   const std::string& type2,
                                   const std::string& type3) {
    auto key1 = std::make_tuple(type1, type2, type3);
    auto key2 = std::make_tuple(type3, type2, type1);
    return has_parameter(container, key1) || has_parameter(container, key2);
}

/**
 * @brief Container size helper
 */
template<typename Container>
inline size_t get_container_size(const Container& container) {
    return container.size();
}

/**
 * @brief Check if container is empty
 */
template<typename Container>
inline bool is_container_empty(const Container& container) {
    return container.empty();
}

/**
 * @brief Get common statistics from storage
 */
struct ParameterCounts {
    size_t atom_masses = 0;
    size_t lj_params = 0;
    size_t nbfix = 0;
    size_t bonds = 0;
    size_t angles = 0;
    size_t dihedrals = 0;
    size_t impropers = 0;
};

inline ParameterCounts get_parameter_counts(const ForceFieldStorage& storage) {
    ParameterCounts counts;
    counts.atom_masses = get_container_size(storage.atom_masses);
    counts.lj_params = get_container_size(storage.lj_params);
    counts.nbfix = get_container_size(storage.nbfix);
    counts.bonds = get_container_size(storage.bond_params);
    counts.angles = get_container_size(storage.angle_params);
    counts.dihedrals = get_container_size(storage.dihedral_params);
    counts.impropers = get_container_size(storage.improper_params);
    return counts;
}

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_OPERATION_UTILS_HPP 