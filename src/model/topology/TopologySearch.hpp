#pragma once

#ifndef PYGCMC_MODEL_TOPOLOGY_SEARCH_HPP
#define PYGCMC_MODEL_TOPOLOGY_SEARCH_HPP

#include "TopologyCore.hpp"
#include <algorithm>
#include <optional>

namespace pygcmc {
namespace model {
namespace topology {

/**
 * @brief Topology search and finder utilities
 */
class TopologySearcher {
public:
    TopologySearcher(const TopologyStorage& storage) : storage_(storage) {}

    /**
     * @brief Get atoms in a specific residue
     */
    std::vector<int> get_residue_atoms(int residue_id) const {
        if (residue_id < 0 || residue_id >= static_cast<int>(storage_.residues.size())) {
            return {};
        }
        return storage_.residues[residue_id].atoms;
    }

    /**
     * @brief Get residues in a specific segment
     */
    std::vector<int> get_segment_residues(int segment_id) const {
        if (segment_id < 0 || segment_id >= static_cast<int>(storage_.segments.size())) {
            return {};
        }
        return storage_.segments[segment_id].residues;
    }

    /**
     * @brief Get bonded neighbors of an atom
     */
    std::vector<int> get_bonded_atoms(int atom_id) const {
        std::vector<int> neighbors;
        neighbors.reserve(8); // Typical coordination number
        for (const auto& bond : storage_.bonds) {
            if (bond.atom1 == atom_id) {
                neighbors.push_back(bond.atom2);
            } else if (bond.atom2 == atom_id) {
                neighbors.push_back(bond.atom1);
            }
        }
        return neighbors;
    }

    /**
     * @brief Find atom by residue and atom name
     */
    std::optional<int> find_atom(const std::string& residue_name, int residue_number,
                                const std::string& atom_name) const {
        for (const auto& segment : storage_.segments) {
            auto it = storage_.atom_map.find(std::make_tuple(residue_name, residue_number, segment.name, atom_name));
            if (it != storage_.atom_map.end()) {
                return it->second;
            }
        }
        return std::nullopt;
    }

    /**
     * @brief Find residue by name and number
     */
    std::optional<int> find_residue(const std::string& name, int number) const {
        for (const auto& segment : storage_.segments) {
            auto it = storage_.residue_map.find(std::make_tuple(name, number, segment.name));
            if (it != storage_.residue_map.end()) {
                return it->second;
            }
        }
        return std::nullopt;
    }

    /**
     * @brief Find segment by name
     */
    std::optional<int> find_segment(const std::string& name) const {
        auto it = storage_.segment_map.find(name);
        if (it != storage_.segment_map.end()) {
            return it->second;
        }
        return std::nullopt;
    }

    /**
     * @brief Check if two atoms are bonded
     */
    bool are_bonded(int atom1, int atom2) const {
        return std::any_of(storage_.bonds.begin(), storage_.bonds.end(),
            [atom1, atom2](const TopologyBond& bond) {
                return (bond.atom1 == atom1 && bond.atom2 == atom2) ||
                       (bond.atom1 == atom2 && bond.atom2 == atom1);
            });
    }

    /**
     * @brief Get titles
     */
    const std::vector<std::string>& get_titles() const { 
        return storage_.titles; 
    }

private:
    const TopologyStorage& storage_;
};

} // namespace topology
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_TOPOLOGY_SEARCH_HPP 