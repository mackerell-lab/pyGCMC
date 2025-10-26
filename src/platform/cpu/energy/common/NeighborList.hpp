#pragma once

#include <vector>
#include <cstddef>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Simple Verlet neighbor list for nonbonded calculations
 *
 * Stores a list of neighboring atoms within a cutoff distance for each atom.
 * Reduces O(N²) all-pairs calculations to O(N*neighbors).
 */
class NeighborList {
public:
    /**
     * @brief Construct empty neighbor list
     */
    NeighborList() : cutoff_(0.0f), cutoff_sq_(0.0f), num_atoms_(0), rebuild_counter_(0) {}

    /**
     * @brief Initialize neighbor list with cutoff
     * @param cutoff Cutoff distance for neighbor list (nm)
     * @param num_atoms Number of atoms
     */
    void initialize(float cutoff, size_t num_atoms);

    /**
     * @brief Build neighbor list from atom coordinates
     * @param coords Atom coordinates [x1,y1,z1, x2,y2,z2, ...]
     * @param box Box dimensions [x, y, z] (nm)
     * @param use_pbc Whether to use periodic boundary conditions
     */
    void build(const std::vector<float>& coords, const std::vector<float>& box, bool use_pbc);

    /**
     * @brief Get neighbors for a specific atom
     * @param atom_idx Atom index
     * @return Vector of neighbor atom indices
     */
    const std::vector<int>& getNeighbors(int atom_idx) const {
        if (atom_idx < 0 || static_cast<size_t>(atom_idx) >= neighbors_.size()) {
            static const std::vector<int> empty;
            return empty;
        }
        return neighbors_[atom_idx];
    }

    /**
     * @brief Check if neighbor list is empty/uninitialized
     */
    bool empty() const { return neighbors_.empty(); }

    /**
     * @brief Get number of atoms
     */
    size_t numAtoms() const { return num_atoms_; }

    /**
     * @brief Get cutoff distance
     */
    float getCutoff() const { return cutoff_; }

    /**
     * @brief Get total number of neighbor pairs
     */
    size_t getTotalNeighborPairs() const;

    /**
     * @brief Increment rebuild counter
     */
    void incrementRebuildCounter() { rebuild_counter_++; }

    /**
     * @brief Get rebuild counter value
     */
    unsigned int getRebuildCounter() const { return rebuild_counter_; }

    /**
     * @brief Reset rebuild counter
     */
    void resetRebuildCounter() { rebuild_counter_ = 0; }

private:
    float cutoff_;                              ///< Neighbor list cutoff (nm)
    float cutoff_sq_;                           ///< Squared cutoff for faster distance checks
    size_t num_atoms_;                          ///< Number of atoms
    std::vector<std::vector<int>> neighbors_;   ///< neighbors_[i] = list of neighbor indices for atom i
    unsigned int rebuild_counter_;              ///< Counter for tracking rebuilds
};

/**
 * @brief Calculate minimum image distance squared with PBC
 * @param dx Distance vector component
 * @param box_size Box size in that dimension
 * @return Minimum image distance component
 */
inline float applyPBC(float dx, float box_size) {
    if (dx > box_size * 0.5f) {
        dx -= box_size;
    } else if (dx < -box_size * 0.5f) {
        dx += box_size;
    }
    return dx;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
