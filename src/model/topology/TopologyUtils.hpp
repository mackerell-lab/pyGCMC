#pragma once

#ifndef PYGCMC_MODEL_TOPOLOGY_UTILS_HPP
#define PYGCMC_MODEL_TOPOLOGY_UTILS_HPP

#include "TopologyCore.hpp"
#include <algorithm>
#include <sstream>
#include <map>

namespace pygcmc {
namespace model {
namespace topology {

/**
 * @brief Topology utilities and analysis functions
 */
class TopologyUtils {
public:
    /**
     * @brief Topology statistics structure
     */
    struct TopologyStats {
        size_t num_atoms = 0;
        size_t num_residues = 0;
        size_t num_segments = 0;
        size_t num_bonds = 0;
        size_t num_angles = 0;
        size_t num_dihedrals = 0;
        size_t num_impropers = 0;
        size_t num_donors = 0;
        size_t num_acceptors = 0;
        size_t num_groups = 0;
        size_t num_cmaps = 0;
        size_t num_exclusions = 0;
    };

    TopologyUtils(const TopologyStorage& storage) : storage_(storage) {}

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
        const int num_atoms = static_cast<int>(storage_.atoms.size());
        for (const auto& bond : storage_.bonds) {
            if (bond.atom1 < 0 || bond.atom1 >= num_atoms ||
                bond.atom2 < 0 || bond.atom2 >= num_atoms) {
                return false;
            }
        }

        // Check angles reference valid atoms
        for (const auto& angle : storage_.angles) {
            if (angle.atom1 < 0 || angle.atom1 >= num_atoms ||
                angle.atom2 < 0 || angle.atom2 >= num_atoms ||
                angle.atom3 < 0 || angle.atom3 >= num_atoms) {
                return false;
            }
        }

        // Check dihedrals reference valid atoms
        for (const auto& dihedral : storage_.dihedrals) {
            if (dihedral.atom1 < 0 || dihedral.atom1 >= num_atoms ||
                dihedral.atom2 < 0 || dihedral.atom2 >= num_atoms ||
                dihedral.atom3 < 0 || dihedral.atom3 >= num_atoms ||
                dihedral.atom4 < 0 || dihedral.atom4 >= num_atoms) {
                return false;
            }
        }

        return true;
    }

    /**
     * @brief Get comprehensive topology statistics
     */
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
     * @brief Get human-readable topology summary
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

    /**
     * @brief Find shortest path between atoms (bond distance)
     */
    int get_bond_distance(int atom1, int atom2) const {
        if (atom1 == atom2) return 0;
        
        const size_t num_atoms = storage_.atoms.size();
        std::vector<bool> visited(num_atoms, false);
        std::vector<int> queue, distance(num_atoms, -1);
        
        queue.push_back(atom1);
        visited[atom1] = true;
        distance[atom1] = 0;
        
        for (size_t front = 0; front < queue.size(); ++front) {
            int current = queue[front];
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
     * @brief Check for CMAP with specific atoms
     */
    bool has_cmap(const std::vector<int>& atoms) const {
        if (atoms.size() < 5) return false;
        
        return std::any_of(storage_.cmaps.begin(), storage_.cmaps.end(),
            [&atoms](const TopologyCmap& cmap) {
                for (size_t i = 0; i < 5 && i < atoms.size(); ++i) {
                    if (cmap.atoms[i] != atoms[i]) return false;
                }
                return true;
            });
    }

    /**
     * @brief Get titles
     */
    const std::vector<std::string>& get_titles() const { 
        return storage_.titles; 
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
     * @brief Check if topology element exists
     */
    bool has_atom(int index) const {
        return index >= 0 && index < static_cast<int>(storage_.atoms.size());
    }

    bool has_residue(int index) const {
        return index >= 0 && index < static_cast<int>(storage_.residues.size());
    }

    bool has_segment(int index) const {
        return index >= 0 && index < static_cast<int>(storage_.segments.size());
    }

    /**
     * @brief Check existence of bonds, angles, dihedrals
     */
    bool has_bond(int atom1, int atom2) const {
        return std::any_of(storage_.bonds.begin(), storage_.bonds.end(),
            [atom1, atom2](const TopologyBond& bond) {
                return (bond.atom1 == atom1 && bond.atom2 == atom2) ||
                       (bond.atom1 == atom2 && bond.atom2 == atom1);
            });
    }

    bool has_angle(int atom1, int atom2, int atom3) const {
        return std::any_of(storage_.angles.begin(), storage_.angles.end(),
            [atom1, atom2, atom3](const TopologyAngle& angle) {
                return (angle.atom1 == atom1 && angle.atom2 == atom2 && angle.atom3 == atom3) ||
                       (angle.atom1 == atom3 && angle.atom2 == atom2 && angle.atom3 == atom1);
            });
    }

    bool has_dihedral(int atom1, int atom2, int atom3, int atom4) const {
        return std::any_of(storage_.dihedrals.begin(), storage_.dihedrals.end(),
            [atom1, atom2, atom3, atom4](const TopologyDihedral& dihedral) {
                return !dihedral.improper &&
                       ((dihedral.atom1 == atom1 && dihedral.atom2 == atom2 && 
                         dihedral.atom3 == atom3 && dihedral.atom4 == atom4) ||
                        (dihedral.atom1 == atom4 && dihedral.atom2 == atom3 && 
                         dihedral.atom3 == atom2 && dihedral.atom4 == atom1));
            });
    }

    bool has_improper(int atom1, int atom2, int atom3, int atom4) const {
        return std::any_of(storage_.dihedrals.begin(), storage_.dihedrals.end(),
            [atom1, atom2, atom3, atom4](const TopologyDihedral& dihedral) {
                return dihedral.improper &&
                       dihedral.atom1 == atom1 && dihedral.atom2 == atom2 && 
                       dihedral.atom3 == atom3 && dihedral.atom4 == atom4;
            });
    }

    /**
     * @brief Check advanced topology features
     */
    bool has_donor(int donor_atom) const {
        return std::any_of(storage_.donors.begin(), storage_.donors.end(),
            [donor_atom](const TopologyDonor& d) { return d.donor_atom == donor_atom; });
    }
    
    bool has_donor(int donor_atom, int hydrogen_atom) const {
        return std::any_of(storage_.donors.begin(), storage_.donors.end(),
            [donor_atom, hydrogen_atom](const TopologyDonor& d) { 
                return d.donor_atom == donor_atom && d.hydrogen_atom == hydrogen_atom; 
            });
    }
    
    bool has_acceptor(int acceptor_atom) const {
        return std::any_of(storage_.acceptors.begin(), storage_.acceptors.end(),
            [acceptor_atom](const TopologyAcceptor& a) { return a.acceptor_atom == acceptor_atom; });
    }
    
    bool has_group(int group_id) const {
        return std::any_of(storage_.groups.begin(), storage_.groups.end(),
            [group_id](const TopologyGroup& g) { return g.id == group_id; });
    }

private:
    const TopologyStorage& storage_;
};

} // namespace topology
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_TOPOLOGY_UTILS_HPP 