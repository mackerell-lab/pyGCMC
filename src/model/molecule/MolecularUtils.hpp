#pragma once

#ifndef PYGCMC_MODEL_MOLECULAR_UTILS_HPP
#define PYGCMC_MODEL_MOLECULAR_UTILS_HPP

#include <vector>
#include <memory>
#include <functional>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include "../atom/AtomMain.hpp"
#include "../residue/ResidueMain.hpp"
#include "../common/ModelUtils.hpp"

namespace pygcmc {
namespace model {

// Forward declaration
class MolecularComposite;

/**
 * @brief Selection criteria for atoms and residues
 */
namespace selection {
    
    // Atom selection predicates
    inline std::function<bool(const Atom&)> by_element(const std::string& element) {
        return [element](const Atom& atom) {
            return atom.get_element() == element;
        };
    }
    
    inline std::function<bool(const Atom&)> by_atom_type(const std::string& type) {
        return [type](const Atom& atom) {
            return atom.get_type() == type;
        };
    }
    
    inline std::function<bool(const Atom&)> by_residue_name(const std::string& resname) {
        return [resname](const Atom& atom) {
            return atom.get_resname() == resname;
        };
    }
    
    inline std::function<bool(const Atom&)> by_chain(char chain) {
        return [chain](const Atom& atom) {
            return atom.get_chain() == chain;
        };
    }
    
    inline std::function<bool(const Atom&)> by_segment(const std::string& segid) {
        return [segid](const Atom& atom) {
            return atom.get_segid() == segid;
        };
    }
    
    inline std::function<bool(const Atom&)> is_hetatm() {
        return [](const Atom& atom) {
            return atom.is_hetatm();
        };
    }
    
    inline std::function<bool(const Atom&)> within_sphere(
        const std::array<double, 3>& center, double radius) {
        return [center, radius](const Atom& atom) {
            return utils::math::distance(atom.get_coor(), center) <= radius;
        };
    }
    
    inline std::function<bool(const Atom&)> charge_range(double min_charge, double max_charge) {
        return [min_charge, max_charge](const Atom& atom) {
            double charge = atom.get_charge();
            return charge >= min_charge && charge <= max_charge;
        };
    }
    
    // Residue selection predicates
    inline std::function<bool(const Residue&)> residue_by_name(const std::string& resname) {
        return [resname](const Residue& residue) {
            return residue.get_resname() == resname;
        };
    }
    
    inline std::function<bool(const Residue&)> residue_by_chain(char chain) {
        return [chain](const Residue& residue) {
            return residue.get_chain() == chain;
        };
    }
    
    inline std::function<bool(const Residue&)> residue_by_number_range(int min_num, int max_num) {
        return [min_num, max_num](const Residue& residue) {
            int resnum = residue.get_ires();
            return resnum >= min_num && resnum <= max_num;
        };
    }
    
    inline std::function<bool(const Residue&)> residue_is_hetatm() {
        return [](const Residue& residue) {
            return residue.is_hetatm();
        };
    }
    
    // Compound predicates
    template<typename... Predicates>
    inline std::function<bool(const Atom&)> atom_and(Predicates... preds) {
        return [preds...](const Atom& atom) {
            return (preds(atom) && ...);
        };
    }
    
    template<typename... Predicates>
    inline std::function<bool(const Atom&)> atom_or(Predicates... preds) {
        return [preds...](const Atom& atom) {
            return (preds(atom) || ...);
        };
    }

} // namespace selection

/**
 * @brief Molecular analysis and manipulation utilities
 */
class MolecularUtils {
public:
    // Atom selection methods
    static std::vector<std::shared_ptr<Atom>> select_atoms(
        const MolecularComposite& mol,
        const std::function<bool(const Atom&)>& predicate);

    static std::vector<std::shared_ptr<Atom>> select_atoms_by_type(
        const MolecularComposite& mol,
        const std::vector<std::string>& types);

    static std::vector<std::shared_ptr<Atom>> select_atoms_by_element(
        const MolecularComposite& mol,
        const std::vector<std::string>& elements);

    static std::vector<std::shared_ptr<Atom>> select_atoms_within_distance(
        const MolecularComposite& mol,
        const std::array<double, 3>& center,
        double distance);

    // Residue selection methods
    static std::vector<std::shared_ptr<Residue>> select_residues(
        const MolecularComposite& mol,
        const std::function<bool(const Residue&)>& predicate);

    static std::vector<std::shared_ptr<Residue>> select_residues_by_chain(
        const MolecularComposite& mol,
        char chain_id);

    static std::vector<std::shared_ptr<Residue>> select_residues_by_names(
        const MolecularComposite& mol,
        const std::vector<std::string>& names);

    // Grouping methods
    static std::unordered_map<std::string, std::vector<std::shared_ptr<Atom>>> 
    group_atoms_by_residue_name(const MolecularComposite& mol);

    static std::unordered_map<char, std::vector<std::shared_ptr<Atom>>> 
    group_atoms_by_chain(const MolecularComposite& mol);

    static std::unordered_map<std::string, std::vector<std::shared_ptr<Atom>>> 
    group_atoms_by_element(const MolecularComposite& mol);

    static std::unordered_map<char, std::vector<std::shared_ptr<Residue>>> 
    group_residues_by_chain(const MolecularComposite& mol);

    // Distance and geometry analysis
    static std::vector<std::pair<std::shared_ptr<Atom>, double>> 
    find_atoms_within_distance(
        const MolecularComposite& mol,
        const std::array<double, 3>& center,
        double max_distance);

    static std::vector<std::pair<std::shared_ptr<Atom>, std::shared_ptr<Atom>>> 
    find_atom_pairs_within_distance(
        const MolecularComposite& mol,
        double max_distance);

    static double calculate_minimum_distance(
        const std::vector<std::shared_ptr<Atom>>& atoms1,
        const std::vector<std::shared_ptr<Atom>>& atoms2);

    static double calculate_center_distance(
        const std::vector<std::shared_ptr<Atom>>& atoms1,
        const std::vector<std::shared_ptr<Atom>>& atoms2);

    // Statistical analysis
    struct MolecularStatistics {
        size_t total_atoms = 0;
        size_t total_residues = 0;
        size_t hetatm_count = 0;
        double total_mass = 0.0;
        double total_charge = 0.0;
        std::array<double, 3> center_of_mass = {0.0, 0.0, 0.0};
        std::array<double, 3> geometric_center = {0.0, 0.0, 0.0};
        std::unordered_map<std::string, size_t> element_counts;
        std::unordered_map<std::string, size_t> residue_counts;
        std::unordered_map<char, size_t> chain_counts;
    };

    static MolecularStatistics calculate_statistics(const MolecularComposite& mol);

    // Connectivity analysis
    static std::vector<std::vector<size_t>> find_connected_components(
        const MolecularComposite& mol,
        double bond_distance_threshold = 2.0);

    static bool are_atoms_bonded(
        const Atom& atom1, 
        const Atom& atom2,
        double bond_distance_threshold = 2.0);

    // Transformation utilities
    static void translate_atoms(
        const std::vector<std::shared_ptr<Atom>>& atoms,
        const std::array<double, 3>& translation);

    static void center_molecule(MolecularComposite& mol);

    static void align_molecule_to_axis(
        MolecularComposite& mol,
        const std::array<double, 3>& axis);

    // Validation utilities
    static std::vector<std::string> check_molecular_consistency(const MolecularComposite& mol);
    
    static bool has_duplicate_atoms(const MolecularComposite& mol);
    
    static std::vector<std::pair<size_t, size_t>> find_overlapping_atoms(
        const MolecularComposite& mol,
        double overlap_threshold = 0.5);

    // Export utilities
    static std::string to_xyz_format(const MolecularComposite& mol);
    
    static std::string to_pdb_format(const MolecularComposite& mol);

    // Molecular property calculations
    static double calculate_radius_of_gyration(const MolecularComposite& mol);
    
    static std::array<double, 6> calculate_bounding_box(const MolecularComposite& mol);
    
    static double calculate_molecular_volume(
        const MolecularComposite& mol,
        double probe_radius = 1.4);

private:
    // Helper methods
    static std::array<double, 3> calculate_geometric_center(
        const std::vector<std::shared_ptr<Atom>>& atoms);
    
    static double estimate_atomic_radius(const std::string& element);
};

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_MOLECULAR_UTILS_HPP 