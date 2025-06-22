#pragma once

#ifndef PYGCMC_MODEL_TOPOLOGY_UTILS_HPP
#define PYGCMC_MODEL_TOPOLOGY_UTILS_HPP

#include "TopologyCore.hpp"
#include <algorithm>
#include <sstream>

namespace pygcmc {
namespace model {
namespace topology {

/**
 * @brief Topology utilities and helper functions
 */
class TopologyUtils {
public:
    TopologyUtils(const TopologyStorage& storage) : storage_(storage) {}

    /**
     * @brief Add title to topology
     */
    static void add_title(TopologyStorage& storage, const std::string& title) {
        storage.titles.push_back(title);
    }

    /**
     * @brief Get titles
     */
    const std::vector<std::string>& get_titles() const { 
        return storage_.titles; 
    }

    /**
     * @brief Check for specific connectivity patterns
     */
    bool has_cmap(const std::vector<int>& atoms) const {
        if (atoms.size() < 5) return false;
        
        return std::any_of(storage_.cmaps.begin(), storage_.cmaps.end(),
            [&atoms](const TopologyCmap& cmap) {
                // Check if first 5 atoms match
                for (size_t i = 0; i < 5 && i < atoms.size(); ++i) {
                    if (cmap.atoms[i] != atoms[i]) return false;
                }
                return true;
            });
    }

    /**
     * @brief Validate topology consistency
     */
    bool validate_topology() const {
        // Check basic consistency
        for (const auto& atom : storage_.atoms) {
            if (atom.residue_id < 0 || atom.residue_id >= static_cast<int>(storage_.residues.size()) ||
                atom.segment_id < 0 || atom.segment_id >= static_cast<int>(storage_.segments.size())) {
                return false;
            }
        }

        // Check bonds reference valid atoms
        for (const auto& bond : storage_.bonds) {
            if (bond.atom1 < 0 || bond.atom1 >= static_cast<int>(storage_.atoms.size()) ||
                bond.atom2 < 0 || bond.atom2 >= static_cast<int>(storage_.atoms.size())) {
                return false;
            }
        }

        // Check angles reference valid atoms
        for (const auto& angle : storage_.angles) {
            if (angle.atom1 < 0 || angle.atom1 >= static_cast<int>(storage_.atoms.size()) ||
                angle.atom2 < 0 || angle.atom2 >= static_cast<int>(storage_.atoms.size()) ||
                angle.atom3 < 0 || angle.atom3 >= static_cast<int>(storage_.atoms.size())) {
                return false;
            }
        }

        // Check dihedrals reference valid atoms
        for (const auto& dihedral : storage_.dihedrals) {
            if (dihedral.atom1 < 0 || dihedral.atom1 >= static_cast<int>(storage_.atoms.size()) ||
                dihedral.atom2 < 0 || dihedral.atom2 >= static_cast<int>(storage_.atoms.size()) ||
                dihedral.atom3 < 0 || dihedral.atom3 >= static_cast<int>(storage_.atoms.size()) ||
                dihedral.atom4 < 0 || dihedral.atom4 >= static_cast<int>(storage_.atoms.size())) {
                return false;
            }
        }

        return true;
    }

    /**
     * @brief Get topology statistics
     */
    struct TopologyStats {
        size_t num_atoms;
        size_t num_residues;
        size_t num_segments;
        size_t num_bonds;
        size_t num_angles;
        size_t num_dihedrals;
        size_t num_impropers;
        size_t num_donors;
        size_t num_acceptors;
        size_t num_groups;
        size_t num_cmaps;
        size_t num_exclusions;
    };

    TopologyStats get_statistics() const {
        TopologyStats stats;
        stats.num_atoms = storage_.atoms.size();
        stats.num_residues = storage_.residues.size();
        stats.num_segments = storage_.segments.size();
        stats.num_bonds = storage_.bonds.size();
        stats.num_angles = storage_.angles.size();
        
        stats.num_dihedrals = 0;
        stats.num_impropers = 0;
        for (const auto& dih : storage_.dihedrals) {
            if (dih.improper) {
                stats.num_impropers++;
            } else {
                stats.num_dihedrals++;
            }
        }
        
        stats.num_donors = storage_.donors.size();
        stats.num_acceptors = storage_.acceptors.size();
        stats.num_groups = storage_.groups.size();
        stats.num_cmaps = storage_.cmaps.size();
        
        stats.num_exclusions = 0;
        for (const auto& pair : storage_.exclusions) {
            stats.num_exclusions += pair.second.size();
        }
        stats.num_exclusions /= 2; // Each exclusion is counted twice
        
        return stats;
    }

    /**
     * @brief Get topology summary string
     */
    std::string get_summary() const {
        auto stats = get_statistics();
        std::stringstream ss;
        ss << "Topology Summary:\n";
        ss << "  Atoms: " << stats.num_atoms << "\n";
        ss << "  Residues: " << stats.num_residues << "\n";
        ss << "  Segments: " << stats.num_segments << "\n";
        ss << "  Bonds: " << stats.num_bonds << "\n";
        ss << "  Angles: " << stats.num_angles << "\n";
        ss << "  Dihedrals: " << stats.num_dihedrals << "\n";
        ss << "  Impropers: " << stats.num_impropers << "\n";
        if (stats.num_donors > 0) ss << "  H-bond donors: " << stats.num_donors << "\n";
        if (stats.num_acceptors > 0) ss << "  H-bond acceptors: " << stats.num_acceptors << "\n";
        if (stats.num_groups > 0) ss << "  Groups: " << stats.num_groups << "\n";
        if (stats.num_cmaps > 0) ss << "  CMAP terms: " << stats.num_cmaps << "\n";
        if (stats.num_exclusions > 0) ss << "  Exclusions: " << stats.num_exclusions << "\n";
        return ss.str();
    }

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
     * @brief Get all atoms in a specific segment
     */
    std::vector<int> get_segment_atoms(int segment_id) const {
        std::vector<int> atoms;
        auto residue_ids = get_segment_residues(segment_id);
        for (int res_id : residue_ids) {
            auto res_atoms = get_residue_atoms(res_id);
            atoms.insert(atoms.end(), res_atoms.begin(), res_atoms.end());
        }
        return atoms;
    }

    /**
     * @brief Get bonds involving a specific atom
     */
    std::vector<int> get_atom_bonds(int atom_id) const {
        std::vector<int> bond_indices;
        for (size_t i = 0; i < storage_.bonds.size(); ++i) {
            if (storage_.bonds[i].atom1 == atom_id || storage_.bonds[i].atom2 == atom_id) {
                bond_indices.push_back(static_cast<int>(i));
            }
        }
        return bond_indices;
    }

    /**
     * @brief Get bonded neighbors of an atom
     */
    std::vector<int> get_bonded_atoms(int atom_id) const {
        std::vector<int> neighbors;
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
     * @brief Find shortest path between two atoms (number of bonds)
     */
    int get_bond_distance(int atom1, int atom2) const {
        if (atom1 == atom2) return 0;
        
        std::vector<bool> visited(storage_.atoms.size(), false);
        std::vector<int> queue;
        std::vector<int> distance(storage_.atoms.size(), -1);
        
        queue.push_back(atom1);
        visited[atom1] = true;
        distance[atom1] = 0;
        
        size_t front = 0;
        while (front < queue.size()) {
            int current = queue[front++];
            
            auto neighbors = get_bonded_atoms(current);
            for (int neighbor : neighbors) {
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    distance[neighbor] = distance[current] + 1;
                    queue.push_back(neighbor);
                    
                    if (neighbor == atom2) {
                        return distance[neighbor];
                    }
                }
            }
        }
        
        return -1; // Not connected
    }

    /**
     * @brief Count atoms by type
     */
    std::map<std::string, int> count_atom_types() const {
        std::map<std::string, int> counts;
        for (const auto& atom : storage_.atoms) {
            counts[atom.type]++;
        }
        return counts;
    }

    /**
     * @brief Count residues by name
     */
    std::map<std::string, int> count_residue_types() const {
        std::map<std::string, int> counts;
        for (const auto& residue : storage_.residues) {
            counts[residue.name]++;
        }
        return counts;
    }

private:
    const TopologyStorage& storage_;
};

} // namespace topology
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_TOPOLOGY_UTILS_HPP 