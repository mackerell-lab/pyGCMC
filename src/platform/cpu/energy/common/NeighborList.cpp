#include "NeighborList.hpp"
#include <cmath>
#include <algorithm>

namespace pygcmc {
namespace platform {
namespace cpu {

void NeighborList::initialize(float cutoff, size_t num_atoms) {
    cutoff_ = cutoff;
    cutoff_sq_ = cutoff * cutoff;
    num_atoms_ = num_atoms;

    // Allocate neighbor vectors
    neighbors_.clear();
    neighbors_.resize(num_atoms);

    rebuild_counter_ = 0;
}

void NeighborList::build(const std::vector<float>& coords, const std::vector<float>& box, bool use_pbc) {
    if (coords.size() < num_atoms_ * 3) {
        return; // Invalid input
    }

    // Clear existing neighbor lists
    for (auto& nlist : neighbors_) {
        nlist.clear();
    }

    // Build neighbor list: for each atom, find all atoms within cutoff
    for (size_t i = 0; i < num_atoms_; ++i) {
        float xi = coords[i * 3];
        float yi = coords[i * 3 + 1];
        float zi = coords[i * 3 + 2];

        for (size_t j = i + 1; j < num_atoms_; ++j) {
            float xj = coords[j * 3];
            float yj = coords[j * 3 + 1];
            float zj = coords[j * 3 + 2];

            // Calculate distance with optional PBC
            float dx = xj - xi;
            float dy = yj - yi;
            float dz = zj - zi;

            if (use_pbc && !box.empty() && box.size() >= 3) {
                dx = applyPBC(dx, box[0]);
                dy = applyPBC(dy, box[1]);
                dz = applyPBC(dz, box[2]);
            }

            float dist_sq = dx * dx + dy * dy + dz * dz;

            // If within cutoff, add to both neighbor lists
            if (dist_sq < cutoff_sq_) {
                neighbors_[i].push_back(static_cast<int>(j));
                neighbors_[j].push_back(static_cast<int>(i));
            }
        }
    }

    rebuild_counter_++;
}

size_t NeighborList::getTotalNeighborPairs() const {
    size_t total = 0;
    for (const auto& nlist : neighbors_) {
        total += nlist.size();
    }
    // Each pair is counted twice (in both atoms' lists), so divide by 2
    return total / 2;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
