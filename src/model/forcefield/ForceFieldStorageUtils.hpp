#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_STORAGE_UTILS_HPP
#define PYGCMC_MODEL_FORCEFIELD_STORAGE_UTILS_HPP

#include "ForceFieldStorageCore.hpp"
#include <sstream>
#include <iomanip>
#include <limits>

namespace pygcmc {
namespace model {
namespace forcefield {

// === Implementation of ForceFieldStorageOptimized utility methods ===

inline ForceFieldStorageOptimized::MemoryInfo ForceFieldStorageOptimized::get_memory_info() const {
    MemoryInfo info{};
    
    // Atom types and indices
    info.atom_types_memory = atom_types_.size() * sizeof(std::string) + 
                            atom_type_to_index_.size() * (sizeof(std::string) + sizeof(size_t));
    
    // Parameter arrays
    info.parameter_arrays_memory = atom_masses_.size() * sizeof(double) + 
                                 lj_params_.size() * sizeof(LJParams);
    
    // Hash maps (estimated)
    info.hash_maps_memory = 0;
    if (nbfix_params_) info.hash_maps_memory += nbfix_params_->size() * (sizeof(PairKey) + sizeof(NBFIXParams));
    if (bond_params_) info.hash_maps_memory += bond_params_->size() * (sizeof(PairKey) + sizeof(BondParams));
    
    info.total_memory = info.atom_types_memory + info.parameter_arrays_memory + info.hash_maps_memory;
    
    return info;
}

inline ForceFieldStorage ForceFieldStorageOptimized::to_original_format() const {
    ForceFieldStorage storage;
    
    // Convert atom masses
    for (size_t i = 0; i < atom_types_.size(); ++i) {
        if (atom_masses_[i] != 0.0) {
            storage.atom_masses[atom_types_[i]] = atom_masses_[i];
        }
    }
    
    // Convert LJ parameters
    for (size_t i = 0; i < atom_types_.size(); ++i) {
        if (lj_params_[i].epsilon != 0.0 || lj_params_[i].rmin_half != 0.0) {
            storage.lj_params[atom_types_[i]] = lj_params_[i];
        }
    }
    
    // Convert NBFIX parameters
    if (nbfix_params_) {
        for (const auto& [key, params] : *nbfix_params_) {
            const std::string& type1 = atom_types_[key.type1_idx];
            const std::string& type2 = atom_types_[key.type2_idx];
            auto storage_key = ParamKeyUtils::make_pair(type1, type2);
            storage.nbfix[storage_key] = params;
        }
    }
    
    storage.nonbonded_params = nonbonded_params_;
    
    return storage;
}

/**
 * @brief Utility class for advanced ForceFieldStorageOptimized operations
 * 
 * This class provides additional utility functions for the optimized storage
 * that are separate from the core storage operations.
 */
class ForceFieldStorageOptimizedUtils {
public:
    /**
     * @brief Create an optimized storage from original storage format
     */
    static ForceFieldStorageOptimized from_original_format(const ForceFieldStorage& original) {
        ForceFieldStorageOptimized optimized;
        
        // Copy atom masses
        for (const auto& [type, mass] : original.atom_masses) {
            optimized.set_atom_mass(type, mass);
        }
        
        // Copy LJ parameters
        for (const auto& [type, params] : original.lj_params) {
            optimized.set_lj_params(type, params);
        }
        
        // Copy NBFIX parameters
        for (const auto& [key, params] : original.nbfix) {
            optimized.set_nbfix(key.first, key.second, params);
        }
        
        // Copy nonbonded parameters
        optimized.get_nonbonded_params() = original.nonbonded_params;
        
        return optimized;
    }

    /**
     * @brief Get detailed memory breakdown
     */
    static std::string get_memory_report(const ForceFieldStorageOptimized& storage) {
        auto info = storage.get_memory_info();
        std::ostringstream oss;
        
        oss << "=== Optimized Storage Memory Report ===\n";
        oss << "Atom types memory: " << info.atom_types_memory << " bytes\n";
        oss << "Parameter arrays memory: " << info.parameter_arrays_memory << " bytes\n";
        oss << "Hash maps memory: " << info.hash_maps_memory << " bytes\n";
        oss << "Total memory: " << info.total_memory << " bytes\n";
        
        // Calculate memory efficiency
        size_t original_estimate = estimate_original_memory(storage);
        if (original_estimate > 0) {
            double efficiency = (double)info.total_memory / original_estimate;
            oss << "Memory efficiency: " << std::fixed << std::setprecision(2) 
                << (efficiency * 100) << "% of original\n";
        }
        
        return oss.str();
    }

    /**
     * @brief Compare performance characteristics
     */
    static std::string get_performance_comparison(const ForceFieldStorageOptimized& storage) {
        std::ostringstream oss;
        
        oss << "=== Performance Characteristics ===\n";
        oss << "Cache-friendly arrays: " << storage.get_atom_types().size() << " types\n";
        oss << "Index-based lookups: O(1) for masses and LJ params\n";
        oss << "Hash-based lookups: O(1) average for pair parameters\n";
        oss << "Memory locality: Excellent for sequential access\n";
        
        return oss.str();
    }

    /**
     * @brief Validate storage integrity
     */
    static std::vector<std::string> validate_storage_integrity(const ForceFieldStorageOptimized& storage) {
        std::vector<std::string> issues;
        
        const auto& types = storage.get_atom_types();
        const auto& masses = storage.get_atom_masses_array();
        const auto& lj_params = storage.get_lj_params_array();
        
        // Check array size consistency
        if (types.size() != masses.size()) {
            issues.push_back("Atom types and masses arrays have different sizes");
        }
        if (types.size() != lj_params.size()) {
            issues.push_back("Atom types and LJ params arrays have different sizes");
        }
        
        // Check for duplicate types
        std::set<std::string> unique_types(types.begin(), types.end());
        if (unique_types.size() != types.size()) {
            issues.push_back("Duplicate atom types found in storage");
        }
        
        return issues;
    }

private:
    static size_t estimate_original_memory(const ForceFieldStorageOptimized& storage) {
        // Rough estimate of original storage memory usage
        const auto& types = storage.get_atom_types();
        size_t estimate = 0;
        
        // std::map overhead is roughly 3-4x the data size
        estimate += types.size() * (sizeof(std::string) + sizeof(double)) * 4;  // atom_masses map
        estimate += types.size() * (sizeof(std::string) + sizeof(LJParams)) * 4;  // lj_params map
        
        return estimate;
    }
};

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_STORAGE_UTILS_HPP 