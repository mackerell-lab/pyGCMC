#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_STORAGE_TYPES_HPP
#define PYGCMC_MODEL_FORCEFIELD_STORAGE_TYPES_HPP

#include <cstdint>
#include <functional>

namespace pygcmc {
namespace model {
namespace forcefield {

/**
 * @brief Optimized key for pair parameters using type indices
 * 
 * This structure provides efficient storage and lookup for parameter pairs
 * by using integer indices instead of string keys, reducing memory usage
 * and improving cache performance.
 */
struct PairKey {
    uint32_t type1_idx;  ///< Index of first atom type
    uint32_t type2_idx;  ///< Index of second atom type
    
    /**
     * @brief Equality comparison for PairKey
     */
    bool operator==(const PairKey& other) const {
        return type1_idx == other.type1_idx && type2_idx == other.type2_idx;
    }
    
    /**
     * @brief Inequality comparison for PairKey
     */
    bool operator!=(const PairKey& other) const {
        return !(*this == other);
    }
    
    /**
     * @brief Less-than comparison for ordered containers
     */
    bool operator<(const PairKey& other) const {
        if (type1_idx != other.type1_idx) {
            return type1_idx < other.type1_idx;
        }
        return type2_idx < other.type2_idx;
    }
};

/**
 * @brief Hash function for PairKey
 * 
 * Combines the two 32-bit indices into a single 64-bit hash value
 * for efficient hash table lookups.
 */
struct PairKeyHash {
    std::size_t operator()(const PairKey& key) const {
        // Combine two 32-bit values into one 64-bit hash
        return std::hash<uint64_t>{}((uint64_t(key.type1_idx) << 32) | key.type2_idx);
    }
};

/**
 * @brief Triple key for three-parameter interactions (angles, etc.)
 */
struct TripleKey {
    uint32_t type1_idx;
    uint32_t type2_idx;
    uint32_t type3_idx;
    
    bool operator==(const TripleKey& other) const {
        return type1_idx == other.type1_idx && 
               type2_idx == other.type2_idx && 
               type3_idx == other.type3_idx;
    }
    
    bool operator<(const TripleKey& other) const {
        if (type1_idx != other.type1_idx) return type1_idx < other.type1_idx;
        if (type2_idx != other.type2_idx) return type2_idx < other.type2_idx;
        return type3_idx < other.type3_idx;
    }
};

/**
 * @brief Hash function for TripleKey
 */
struct TripleKeyHash {
    std::size_t operator()(const TripleKey& key) const {
        // Simple hash combination
        std::size_t h1 = std::hash<uint32_t>{}(key.type1_idx);
        std::size_t h2 = std::hash<uint32_t>{}(key.type2_idx);
        std::size_t h3 = std::hash<uint32_t>{}(key.type3_idx);
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

/**
 * @brief Quad key for four-parameter interactions (dihedrals, impropers)
 */
struct QuadKey {
    uint32_t type1_idx;
    uint32_t type2_idx;
    uint32_t type3_idx;
    uint32_t type4_idx;
    
    bool operator==(const QuadKey& other) const {
        return type1_idx == other.type1_idx && 
               type2_idx == other.type2_idx && 
               type3_idx == other.type3_idx && 
               type4_idx == other.type4_idx;
    }
    
    bool operator<(const QuadKey& other) const {
        if (type1_idx != other.type1_idx) return type1_idx < other.type1_idx;
        if (type2_idx != other.type2_idx) return type2_idx < other.type2_idx;
        if (type3_idx != other.type3_idx) return type3_idx < other.type3_idx;
        return type4_idx < other.type4_idx;
    }
};

/**
 * @brief Hash function for QuadKey
 */
struct QuadKeyHash {
    std::size_t operator()(const QuadKey& key) const {
        // Simple hash combination
        std::size_t h1 = std::hash<uint32_t>{}(key.type1_idx);
        std::size_t h2 = std::hash<uint32_t>{}(key.type2_idx);
        std::size_t h3 = std::hash<uint32_t>{}(key.type3_idx);
        std::size_t h4 = std::hash<uint32_t>{}(key.type4_idx);
        return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3);
    }
};

/**
 * @brief Utility functions for creating optimized keys
 */
class StorageKeyUtils {
public:
    /**
     * @brief Create ordered pair key (smaller index first)
     */
    static PairKey make_ordered_pair(uint32_t idx1, uint32_t idx2) {
        return {std::min(idx1, idx2), std::max(idx1, idx2)};
    }
    
    /**
     * @brief Create triple key
     */
    static TripleKey make_triple(uint32_t idx1, uint32_t idx2, uint32_t idx3) {
        return {idx1, idx2, idx3};
    }
    
    /**
     * @brief Create quad key
     */
    static QuadKey make_quad(uint32_t idx1, uint32_t idx2, uint32_t idx3, uint32_t idx4) {
        return {idx1, idx2, idx3, idx4};
    }
    
    /**
     * @brief Convert uint64_t to pair key (for backward compatibility)
     */
    static PairKey from_uint64(uint64_t combined) {
        uint32_t idx1 = static_cast<uint32_t>(combined >> 32);
        uint32_t idx2 = static_cast<uint32_t>(combined & 0xFFFFFFFF);
        return {idx1, idx2};
    }
    
    /**
     * @brief Convert pair key to uint64_t (for backward compatibility)
     */
    static uint64_t to_uint64(const PairKey& key) {
        return (uint64_t(key.type1_idx) << 32) | key.type2_idx;
    }
};

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_STORAGE_TYPES_HPP 