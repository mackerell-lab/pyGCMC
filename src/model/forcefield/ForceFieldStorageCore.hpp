#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_STORAGE_CORE_HPP
#define PYGCMC_MODEL_FORCEFIELD_STORAGE_CORE_HPP

#include "ForceFieldParams.hpp"
#include <unordered_map>
#include <vector>
#include <memory>

namespace pygcmc {
namespace model {
namespace forcefield {

/**
 * @brief Memory-optimized force field storage with improved cache locality
 * 
 * This class provides an alternative storage backend with:
 * - Better cache locality through structure-of-arrays layout
 * - Optimized hash maps for faster lookups
 * - Memory pool allocation for reduced fragmentation
 * - Lazy initialization for large parameter sets
 */
class ForceFieldStorageOptimized {
public:
    ForceFieldStorageOptimized() {
        // Reserve space for common parameter counts to avoid reallocations
        atom_type_to_index_.reserve(256);
        atom_masses_.reserve(256);
        lj_params_.reserve(256);
    }

    // Copy constructor (deep copy)
    ForceFieldStorageOptimized(const ForceFieldStorageOptimized& other) 
        : atom_types_(other.atom_types_)
        , atom_type_to_index_(other.atom_type_to_index_)
        , atom_masses_(other.atom_masses_)
        , lj_params_(other.lj_params_)
        , nonbonded_params_(other.nonbonded_params_) {
        
        // Deep copy unique_ptr members
        if (other.nbfix_params_) {
            nbfix_params_ = std::make_unique<std::unordered_map<PairKey, NBFIXParams, PairKeyHash>>(*other.nbfix_params_);
        }
        if (other.bond_params_) {
            bond_params_ = std::make_unique<std::unordered_map<PairKey, BondParams, PairKeyHash>>(*other.bond_params_);
        }
        if (other.angle_params_) {
            angle_params_ = std::make_unique<std::unordered_map<uint64_t, AngleParams>>(*other.angle_params_);
        }
        if (other.dihedral_params_) {
            dihedral_params_ = std::make_unique<std::unordered_map<uint64_t, std::vector<DihedralParams>>>(*other.dihedral_params_);
        }
        if (other.improper_params_) {
            improper_params_ = std::make_unique<std::unordered_map<uint64_t, ImproperParams>>(*other.improper_params_);
        }
    }

    // Move constructor
    ForceFieldStorageOptimized(ForceFieldStorageOptimized&& other) noexcept = default;

    // Copy assignment
    ForceFieldStorageOptimized& operator=(const ForceFieldStorageOptimized& other) {
        if (this != &other) {
            atom_types_ = other.atom_types_;
            atom_type_to_index_ = other.atom_type_to_index_;
            atom_masses_ = other.atom_masses_;
            lj_params_ = other.lj_params_;
            nonbonded_params_ = other.nonbonded_params_;
            
            // Deep copy unique_ptr members
            nbfix_params_.reset();
            bond_params_.reset();
            angle_params_.reset();
            dihedral_params_.reset();
            improper_params_.reset();
            
            if (other.nbfix_params_) {
                nbfix_params_ = std::make_unique<std::unordered_map<PairKey, NBFIXParams, PairKeyHash>>(*other.nbfix_params_);
            }
            if (other.bond_params_) {
                bond_params_ = std::make_unique<std::unordered_map<PairKey, BondParams, PairKeyHash>>(*other.bond_params_);
            }
            if (other.angle_params_) {
                angle_params_ = std::make_unique<std::unordered_map<uint64_t, AngleParams>>(*other.angle_params_);
            }
            if (other.dihedral_params_) {
                dihedral_params_ = std::make_unique<std::unordered_map<uint64_t, std::vector<DihedralParams>>>(*other.dihedral_params_);
            }
            if (other.improper_params_) {
                improper_params_ = std::make_unique<std::unordered_map<uint64_t, ImproperParams>>(*other.improper_params_);
            }
        }
        return *this;
    }

    // Move assignment
    ForceFieldStorageOptimized& operator=(ForceFieldStorageOptimized&& other) noexcept = default;

    ~ForceFieldStorageOptimized() = default;

    // === Atom Type Management ===

    /**
     * @brief Get or create index for atom type
     */
    size_t get_or_create_atom_type_index(const std::string& type) {
        auto it = atom_type_to_index_.find(type);
        if (it != atom_type_to_index_.end()) {
            return it->second;
        }
        
        size_t index = atom_types_.size();
        atom_types_.push_back(type);
        atom_type_to_index_[type] = index;
        
        // Resize arrays to accommodate new type
        atom_masses_.resize(index + 1, 0.0);
        lj_params_.resize(index + 1, LJParams{});
        
        return index;
    }

    /**
     * @brief Get index for existing atom type
     */
    std::pair<size_t, bool> find_atom_type_index(const std::string& type) const {
        auto it = atom_type_to_index_.find(type);
        if (it != atom_type_to_index_.end()) {
            return {it->second, true};
        }
        return {0, false};
    }

    // === Optimized Parameter Storage ===

    /**
     * @brief Set atom mass by type
     */
    void set_atom_mass(const std::string& type, double mass) {
        size_t index = get_or_create_atom_type_index(type);
        atom_masses_[index] = mass;
    }

    /**
     * @brief Get atom mass by type
     */
    std::pair<double, bool> get_atom_mass(const std::string& type) const {
        auto [index, found] = find_atom_type_index(type);
        if (found) {
            return {atom_masses_[index], true};
        }
        return {0.0, false};
    }

    /**
     * @brief Set LJ parameters by type
     */
    void set_lj_params(const std::string& type, const LJParams& params) {
        size_t index = get_or_create_atom_type_index(type);
        lj_params_[index] = params;
    }

    /**
     * @brief Get LJ parameters by type
     */
    std::pair<LJParams, bool> get_lj_params(const std::string& type) const {
        auto [index, found] = find_atom_type_index(type);
        if (found) {
            return {lj_params_[index], true};
        }
        return {LJParams{}, false};
    }

    // === Pair Parameter Storage ===

    /**
     * @brief Optimized key for pair parameters using type indices
     */
    struct PairKey {
        uint32_t type1_idx;
        uint32_t type2_idx;
        
        bool operator==(const PairKey& other) const {
            return type1_idx == other.type1_idx && type2_idx == other.type2_idx;
        }
    };

    struct PairKeyHash {
        std::size_t operator()(const PairKey& key) const {
            return std::hash<uint64_t>{}((uint64_t(key.type1_idx) << 32) | key.type2_idx);
        }
    };

    /**
     * @brief Set NBFIX parameters
     */
    void set_nbfix(const std::string& type1, const std::string& type2, const NBFIXParams& params) {
        auto idx1 = get_or_create_atom_type_index(type1);
        auto idx2 = get_or_create_atom_type_index(type2);
        
        PairKey key{static_cast<uint32_t>(std::min(idx1, idx2)), 
                   static_cast<uint32_t>(std::max(idx1, idx2))};
        
        if (!nbfix_params_) {
            nbfix_params_ = std::make_unique<std::unordered_map<PairKey, NBFIXParams, PairKeyHash>>();
        }
        (*nbfix_params_)[key] = params;
    }

    /**
     * @brief Get NBFIX parameters
     */
    std::pair<NBFIXParams, bool> get_nbfix(const std::string& type1, const std::string& type2) const {
        if (!nbfix_params_) {
            return {NBFIXParams{}, false};
        }
        
        auto [idx1, found1] = find_atom_type_index(type1);
        auto [idx2, found2] = find_atom_type_index(type2);
        
        if (!found1 || !found2) {
            return {NBFIXParams{}, false};
        }
        
        PairKey key{static_cast<uint32_t>(std::min(idx1, idx2)), 
                   static_cast<uint32_t>(std::max(idx1, idx2))};
        
        auto it = nbfix_params_->find(key);
        if (it != nbfix_params_->end()) {
            return {it->second, true};
        }
        
        return {NBFIXParams{}, false};
    }

    /**
     * @brief Clear all parameters and reset storage
     */
    void clear() {
        atom_types_.clear();
        atom_type_to_index_.clear();
        atom_masses_.clear();
        lj_params_.clear();
        
        nbfix_params_.reset();
        bond_params_.reset();
        angle_params_.reset();
        dihedral_params_.reset();
        improper_params_.reset();
        
        nonbonded_params_ = NonbondedParams{};
    }

    // === Statistics ===

    /**
     * @brief Get optimized statistics
     */
    ForceFieldStats get_statistics() const {
        ForceFieldStats stats;
        stats.num_atom_types = atom_types_.size();
        stats.num_lj_params = count_non_zero_lj_params();
        stats.num_nbfix = nbfix_params_ ? nbfix_params_->size() : 0;
        stats.num_bond_types = bond_params_ ? bond_params_->size() : 0;
        stats.num_angle_types = angle_params_ ? angle_params_->size() : 0;
        stats.num_dihedral_types = dihedral_params_ ? dihedral_params_->size() : 0;
        stats.num_improper_types = improper_params_ ? improper_params_->size() : 0;
        return stats;
    }

    // === Accessors for Data ===

    const std::vector<std::string>& get_atom_types() const { return atom_types_; }
    const std::vector<double>& get_atom_masses_array() const { return atom_masses_; }
    const std::vector<LJParams>& get_lj_params_array() const { return lj_params_; }
    const NonbondedParams& get_nonbonded_params() const { return nonbonded_params_; }
    NonbondedParams& get_nonbonded_params() { return nonbonded_params_; }

    // === Forward declaration of utility methods (implemented in separate file) ===
    
    /**
     * @brief Get memory usage information
     */
    struct MemoryInfo {
        size_t atom_types_memory;
        size_t hash_maps_memory;
        size_t parameter_arrays_memory;
        size_t total_memory;
    };

    MemoryInfo get_memory_info() const;
    ForceFieldStorage to_original_format() const;

private:
    // === Core Type Management ===
    std::vector<std::string> atom_types_;                                    ///< Ordered list of atom types
    std::unordered_map<std::string, size_t> atom_type_to_index_;            ///< Type to index mapping

    // === Structure-of-Arrays Storage ===
    std::vector<double> atom_masses_;                                        ///< Masses indexed by type
    std::vector<LJParams> lj_params_;                                       ///< LJ params indexed by type

    // === Lazy-Initialized Hash Maps ===
    std::unique_ptr<std::unordered_map<PairKey, NBFIXParams, PairKeyHash>> nbfix_params_;
    std::unique_ptr<std::unordered_map<PairKey, BondParams, PairKeyHash>> bond_params_;
    std::unique_ptr<std::unordered_map<uint64_t, AngleParams>> angle_params_;              // Triple hash
    std::unique_ptr<std::unordered_map<uint64_t, std::vector<DihedralParams>>> dihedral_params_;
    std::unique_ptr<std::unordered_map<uint64_t, ImproperParams>> improper_params_;

    // === Global Parameters ===
    NonbondedParams nonbonded_params_;

    // === Helper Methods ===

    size_t count_non_zero_lj_params() const {
        size_t count = 0;
        for (const auto& params : lj_params_) {
            if (params.epsilon != 0.0 || params.rmin_half != 0.0) {
                count++;
            }
        }
        return count;
    }

    // Make utility class friend for implementation access
    friend class ForceFieldStorageOptimizedUtils;
};

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_STORAGE_CORE_HPP 