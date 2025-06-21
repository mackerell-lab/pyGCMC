#pragma once

#ifndef PYGCMC_MODEL_MOLECULE_UTILS_HPP
#define PYGCMC_MODEL_MOLECULE_UTILS_HPP

#include "MolecularComposite.hpp"
#include "../common/ModelUtils.hpp"
#include <vector>
#include <string>
#include <functional>
#include <algorithm>
#include <map>
#include <set>
#include <array>

namespace pygcmc {
namespace model {
namespace molecule {
namespace utils {

/**
 * @brief Selection utilities for molecular systems
 */
namespace selection {

/**
 * @brief Find atoms by predicate function
 */
template<typename Predicate>
std::vector<std::shared_ptr<atom::Atom>> find_atoms_if(
    const MolecularComposite& molecular, 
    Predicate&& predicate) {
    std::vector<std::shared_ptr<atom::Atom>> result;
    const auto& atoms = molecular.get_atoms();
    
    std::copy_if(atoms.begin(), atoms.end(), std::back_inserter(result),
                [&](const std::shared_ptr<atom::Atom>& atom) {
                    return atom && predicate(*atom);
                });
    return result;
}

/**
 * @brief Find residues by predicate function
 */
template<typename Predicate>
std::vector<std::shared_ptr<residue::Residue>> find_residues_if(
    const MolecularComposite& molecular, 
    Predicate&& predicate) {
    std::vector<std::shared_ptr<residue::Residue>> result;
    const auto& residues = molecular.get_residues();
    
    std::copy_if(residues.begin(), residues.end(), std::back_inserter(result),
                [&](const std::shared_ptr<residue::Residue>& residue) {
                    return residue && predicate(*residue);
                });
    return result;
}

/**
 * @brief Find atoms by segment ID
 */
inline std::vector<std::shared_ptr<atom::Atom>> find_atoms_by_segment(
    const MolecularComposite& molecular, 
    const std::string& segment_id) {
    return find_atoms_if(molecular, 
        [&segment_id](const atom::Atom& atom) {
            return atom.get_segid() == segment_id;
        });
}

/**
 * @brief Find atoms by residue name
 */
inline std::vector<std::shared_ptr<atom::Atom>> find_atoms_by_resname(
    const MolecularComposite& molecular, 
    const std::string& resname) {
    return find_atoms_if(molecular, 
        [&resname](const atom::Atom& atom) {
            return atom.get_resname() == resname;
        });
}

/**
 * @brief Find atoms by atom type
 */
inline std::vector<std::shared_ptr<atom::Atom>> find_atoms_by_type(
    const MolecularComposite& molecular, 
    const std::string& atom_type) {
    return find_atoms_if(molecular, 
        [&atom_type](const atom::Atom& atom) {
            return atom.get_type() == atom_type;
        });
}

/**
 * @brief Find atoms in residue range
 */
inline std::vector<std::shared_ptr<atom::Atom>> find_atoms_in_residue_range(
    const MolecularComposite& molecular, 
    int start_residue, int end_residue) {
    return find_atoms_if(molecular, 
        [start_residue, end_residue](const atom::Atom& atom) {
            int ires = atom.get_ires();
            return ires >= start_residue && ires <= end_residue;
        });
}

/**
 * @brief Find residues by segment ID
 */
inline std::vector<std::shared_ptr<residue::Residue>> find_residues_by_segment(
    const MolecularComposite& molecular, 
    const std::string& segment_id) {
    return find_residues_if(molecular, 
        [&segment_id](const residue::Residue& residue) {
            return residue.get_segid() == segment_id;
        });
}

/**
 * @brief Find residues by chain ID
 */
inline std::vector<std::shared_ptr<residue::Residue>> find_residues_by_chain(
    const MolecularComposite& molecular, 
    char chain_id) {
    return find_residues_if(molecular, 
        [chain_id](const residue::Residue& residue) {
            return residue.get_chain() == chain_id;
        });
}

/**
 * @brief Find residues by name
 */
inline std::vector<std::shared_ptr<residue::Residue>> find_residues_by_name(
    const MolecularComposite& molecular, 
    const std::string& resname) {
    return find_residues_if(molecular, 
        [&resname](const residue::Residue& residue) {
            return residue.get_resname() == resname;
        });
}

/**
 * @brief Find protein residues
 */
inline std::vector<std::shared_ptr<residue::Residue>> find_protein_residues(
    const MolecularComposite& molecular) {
    return find_residues_if(molecular, 
        [](const residue::Residue& residue) {
            return residue.is_protein_residue();
        });
}

/**
 * @brief Find nucleic acid residues
 */
inline std::vector<std::shared_ptr<residue::Residue>> find_nucleic_residues(
    const MolecularComposite& molecular) {
    return find_residues_if(molecular, 
        [](const residue::Residue& residue) {
            return residue.is_nucleic_acid_residue();
        });
}

/**
 * @brief Find heavy atoms (non-hydrogen)
 */
inline std::vector<std::shared_ptr<atom::Atom>> find_heavy_atoms(
    const MolecularComposite& molecular) {
    return find_atoms_if(molecular, 
        [](const atom::Atom& atom) {
            const std::string& type = atom.get_type();
            return !type.empty() && type[0] != 'H';
        });
}

/**
 * @brief Find hydrogen atoms
 */
inline std::vector<std::shared_ptr<atom::Atom>> find_hydrogen_atoms(
    const MolecularComposite& molecular) {
    return find_atoms_if(molecular, 
        [](const atom::Atom& atom) {
            const std::string& type = atom.get_type();
            return !type.empty() && type[0] == 'H';
        });
}

} // namespace selection

/**
 * @brief Grouping utilities for molecular systems
 */
namespace grouping {

/**
 * @brief Group atoms by segment ID
 */
inline std::map<std::string, std::vector<std::shared_ptr<atom::Atom>>> 
group_atoms_by_segment(const MolecularComposite& molecular) {
    std::map<std::string, std::vector<std::shared_ptr<atom::Atom>>> groups;
    
    for (const auto& atom : molecular.get_atoms()) {
        if (atom) {
            groups[atom->get_segid()].push_back(atom);
        }
    }
    return groups;
}

/**
 * @brief Group atoms by residue name
 */
inline std::map<std::string, std::vector<std::shared_ptr<atom::Atom>>> 
group_atoms_by_resname(const MolecularComposite& molecular) {
    std::map<std::string, std::vector<std::shared_ptr<atom::Atom>>> groups;
    
    for (const auto& atom : molecular.get_atoms()) {
        if (atom) {
            groups[atom->get_resname()].push_back(atom);
        }
    }
    return groups;
}

/**
 * @brief Group atoms by atom type
 */
inline std::map<std::string, std::vector<std::shared_ptr<atom::Atom>>> 
group_atoms_by_type(const MolecularComposite& molecular) {
    std::map<std::string, std::vector<std::shared_ptr<atom::Atom>>> groups;
    
    for (const auto& atom : molecular.get_atoms()) {
        if (atom) {
            groups[atom->get_type()].push_back(atom);
        }
    }
    return groups;
}

/**
 * @brief Group residues by segment ID
 */
inline std::map<std::string, std::vector<std::shared_ptr<residue::Residue>>> 
group_residues_by_segment(const MolecularComposite& molecular) {
    std::map<std::string, std::vector<std::shared_ptr<residue::Residue>>> groups;
    
    for (const auto& residue : molecular.get_residues()) {
        if (residue) {
            groups[residue->get_segid()].push_back(residue);
        }
    }
    return groups;
}

/**
 * @brief Group residues by chain ID
 */
inline std::map<char, std::vector<std::shared_ptr<residue::Residue>>> 
group_residues_by_chain(const MolecularComposite& molecular) {
    std::map<char, std::vector<std::shared_ptr<residue::Residue>>> groups;
    
    for (const auto& residue : molecular.get_residues()) {
        if (residue) {
            groups[residue->get_chain()].push_back(residue);
        }
    }
    return groups;
}

/**
 * @brief Group residues by name
 */
inline std::map<std::string, std::vector<std::shared_ptr<residue::Residue>>> 
group_residues_by_name(const MolecularComposite& molecular) {
    std::map<std::string, std::vector<std::shared_ptr<residue::Residue>>> groups;
    
    for (const auto& residue : molecular.get_residues()) {
        if (residue) {
            groups[residue->get_resname()].push_back(residue);
        }
    }
    return groups;
}

} // namespace grouping

/**
 * @brief Analysis utilities for molecular systems
 */
namespace analysis {

/**
 * @brief Calculate total molecular mass
 */
inline double calculate_total_mass(const MolecularComposite& molecular) {
    double total = 0.0;
    for (const auto& atom : molecular.get_atoms()) {
        if (atom) {
            total += atom->get_mass();
        }
    }
    return total;
}

/**
 * @brief Calculate total molecular charge
 */
inline double calculate_total_charge(const MolecularComposite& molecular) {
    double total = 0.0;
    for (const auto& atom : molecular.get_atoms()) {
        if (atom) {
            total += atom->get_charge();
        }
    }
    return total;
}

/**
 * @brief Calculate center of mass
 */
inline std::array<double, 3> calculate_center_of_mass(const MolecularComposite& molecular) {
    std::array<double, 3> com = {0.0, 0.0, 0.0};
    double total_mass = 0.0;

    for (const auto& atom : molecular.get_atoms()) {
        if (!atom) continue;
        double mass = atom->get_mass();
        const auto& coor = atom->get_coor();
        for (int i = 0; i < 3; ++i) {
            com[i] += mass * coor[i];
        }
        total_mass += mass;
    }

    if (total_mass > 0.0) {
        for (double& x : com) x /= total_mass;
    }

    return com;
}

/**
 * @brief Calculate geometric center
 */
inline std::array<double, 3> calculate_geometric_center(const MolecularComposite& molecular) {
    std::array<double, 3> center = {0.0, 0.0, 0.0};
    size_t count = 0;

    for (const auto& atom : molecular.get_atoms()) {
        if (!atom) continue;
        const auto& coor = atom->get_coor();
        for (int i = 0; i < 3; ++i) {
            center[i] += coor[i];
        }
        count++;
    }

    if (count > 0) {
        for (double& x : center) x /= static_cast<double>(count);
    }

    return center;
}

/**
 * @brief Calculate bounding box
 */
inline std::pair<std::array<double, 3>, std::array<double, 3>> 
calculate_bounding_box(const MolecularComposite& molecular) {
    std::array<double, 3> min_coords = {
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max()
    };
    std::array<double, 3> max_coords = {
        std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::lowest()
    };

    for (const auto& atom : molecular.get_atoms()) {
        if (!atom) continue;
        const auto& coor = atom->get_coor();
        for (int i = 0; i < 3; ++i) {
            min_coords[i] = std::min(min_coords[i], coor[i]);
            max_coords[i] = std::max(max_coords[i], coor[i]);
        }
    }

    return std::make_pair(min_coords, max_coords);
}

/**
 * @brief Get chain identifiers
 */
inline std::set<char> get_chain_ids(const MolecularComposite& molecular) {
    std::set<char> chains;
    for (const auto& residue : molecular.get_residues()) {
        if (residue) {
            chains.insert(residue->get_chain());
        }
    }
    return chains;
}

/**
 * @brief Get segment identifiers
 */
inline std::set<std::string> get_segment_ids(const MolecularComposite& molecular) {
    std::set<std::string> segments;
    for (const auto& residue : molecular.get_residues()) {
        if (residue) {
            segments.insert(residue->get_segid());
        }
    }
    return segments;
}

/**
 * @brief Get residue names
 */
inline std::set<std::string> get_residue_names(const MolecularComposite& molecular) {
    std::set<std::string> names;
    for (const auto& residue : molecular.get_residues()) {
        if (residue) {
            names.insert(residue->get_resname());
        }
    }
    return names;
}

/**
 * @brief Get atom types
 */
inline std::set<std::string> get_atom_types(const MolecularComposite& molecular) {
    std::set<std::string> types;
    for (const auto& atom : molecular.get_atoms()) {
        if (atom) {
            types.insert(atom->get_type());
        }
    }
    return types;
}

/**
 * @brief Calculate system statistics
 */
struct SystemStatistics {
    size_t num_atoms;
    size_t num_residues;
    size_t num_chains;
    size_t num_segments;
    size_t num_heavy_atoms;
    size_t num_hydrogen_atoms;
    double total_mass;
    double total_charge;
    std::array<double, 3> center_of_mass;
    std::array<double, 3> geometric_center;
    std::pair<std::array<double, 3>, std::array<double, 3>> bounding_box;
};

inline SystemStatistics calculate_system_statistics(const MolecularComposite& molecular) {
    SystemStatistics stats;
    
    stats.num_atoms = molecular.get_num_atoms();
    stats.num_residues = molecular.get_num_residues();
    stats.num_chains = get_chain_ids(molecular).size();
    stats.num_segments = get_segment_ids(molecular).size();
    stats.num_heavy_atoms = selection::find_heavy_atoms(molecular).size();
    stats.num_hydrogen_atoms = selection::find_hydrogen_atoms(molecular).size();
    stats.total_mass = calculate_total_mass(molecular);
    stats.total_charge = calculate_total_charge(molecular);
    stats.center_of_mass = calculate_center_of_mass(molecular);
    stats.geometric_center = calculate_geometric_center(molecular);
    stats.bounding_box = calculate_bounding_box(molecular);
    
    return stats;
}

} // namespace analysis

/**
 * @brief Distance calculation utilities
 */
namespace distance {

/**
 * @brief Calculate distance between two atoms
 */
inline double atom_distance(const atom::Atom& atom1, const atom::Atom& atom2) {
    return atom1.distance_to(atom2);
}

/**
 * @brief Calculate distance between two residues (center of mass)
 */
inline double residue_distance(const residue::Residue& res1, const residue::Residue& res2) {
    return res1.distance_to(res2);
}

/**
 * @brief Find atoms within distance of a point
 */
inline std::vector<std::shared_ptr<atom::Atom>> find_atoms_within_distance(
    const MolecularComposite& molecular,
    const std::array<double, 3>& point,
    double max_distance) {
    
    std::vector<std::shared_ptr<atom::Atom>> result;
    double max_dist_sq = max_distance * max_distance;
    
    for (const auto& atom : molecular.get_atoms()) {
        if (!atom) continue;
        
        const auto& coor = atom->get_coor();
        double dist_sq = 0.0;
        for (int i = 0; i < 3; ++i) {
            double diff = coor[i] - point[i];
            dist_sq += diff * diff;
        }
        
        if (dist_sq <= max_dist_sq) {
            result.push_back(atom);
        }
    }
    
    return result;
}

/**
 * @brief Find atoms within distance of another atom
 */
inline std::vector<std::shared_ptr<atom::Atom>> find_atoms_within_distance(
    const MolecularComposite& molecular,
    const atom::Atom& reference_atom,
    double max_distance) {
    
    return find_atoms_within_distance(molecular, reference_atom.get_coor(), max_distance);
}

} // namespace distance

} // namespace utils
} // namespace molecule
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_MOLECULE_UTILS_HPP 