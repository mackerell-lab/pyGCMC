#pragma once

#ifndef PYGCMC_MODEL_RESIDUE_MAIN_HPP
#define PYGCMC_MODEL_RESIDUE_MAIN_HPP

#include "ResidueComposite.hpp"
#include "ResidueValidator.hpp"
#include "../common/ModelUtils.hpp"
#include <string>
#include <sstream>
#include <iomanip>

namespace pygcmc {
namespace model {
namespace residue {

/**
 * @brief Complete Residue class with full functionality and backward compatibility
 * This class extends ResidueComposite with additional utilities and maintains API compatibility
 */
class Residue : public ResidueComposite {
public:
    // Inherit constructors
    using ResidueComposite::ResidueComposite;

    // PDB format utilities
    std::string get_residue_id() const {
        // Combine residue number and insertion code (e.g., "153A")
        if (inscode == ' ') {
            return std::to_string(ires);
        }
        return std::to_string(ires) + inscode;
    }

    void set_residue_id(const std::string& resid) {
        // Parse residue ID (e.g., "153A" -> ires=153, inscode='A')
        size_t numLen = 0;
        try {
            ires = std::stoi(resid, &numLen);
        } catch (const std::exception&) {
            throw std::invalid_argument("Invalid residue ID format");
        }
        
        if (numLen < resid.length()) {
            inscode = resid[numLen];
        } else {
            inscode = ' ';
        }
    }

    // Enhanced atom lookup methods
    std::shared_ptr<atom::Atom> find_atom_by_pdb_name(const std::string& pdbName) const {
        auto it = std::find_if(atoms.begin(), atoms.end(),
            [&pdbName](const std::shared_ptr<atom::Atom>& atom) {
                return atom && atom->get_formatted_atom_name() == pdbName;
            });
        return (it != atoms.end()) ? *it : nullptr;
    }

    /**
     * @brief Update internal atom mapping for efficient lookups
     */
    void refresh_atom_map() {
        atomMap.clear();
        for (const auto& atom : atoms) {
            if (atom) {
                atomMap[atom->get_type()] = atom;
            }
        }
    }

    /**
     * @brief Update internal atom mapping (for compatibility)
     */
    void update_atom_map() {
        atomMap.clear();
        for (const auto& atom : atoms) {
            if (atom) {
                atomMap[atom->get_type()] = atom;
            }
        }
    }

    // Advanced selection methods
    std::vector<std::shared_ptr<atom::Atom>> select_atoms_by_type(
        const std::vector<std::string>& types) const {
        std::vector<std::shared_ptr<atom::Atom>> selected;
        for (const std::string& type : types) {
            auto atom = find_atom(type);
            if (atom) selected.push_back(atom);
        }
        return selected;
    }

    std::vector<std::shared_ptr<atom::Atom>> select_backbone_atoms() const {
        static const std::vector<std::string> backbone_types = {"N", "CA", "C", "O"};
        return select_atoms_by_type(backbone_types);
    }

    std::vector<std::shared_ptr<atom::Atom>> select_sidechain_atoms() const {
        std::vector<std::shared_ptr<atom::Atom>> sidechain;
        static const std::set<std::string> backbone_set = {"N", "CA", "C", "O", "H", "HA"};
        
        for (const auto& atom : atoms) {
            if (atom && backbone_set.find(atom->get_type()) == backbone_set.end()) {
                sidechain.push_back(atom);
            }
        }
        return sidechain;
    }

    std::vector<std::shared_ptr<atom::Atom>> select_heavy_atoms() const {
        return select_atoms([](const atom::Atom& atom) {
            const std::string& type = atom.get_type();
            return !type.empty() && type[0] != 'H';
        });
    }

    std::vector<std::shared_ptr<atom::Atom>> select_hydrogen_atoms() const {
        return select_atoms([](const atom::Atom& atom) {
            const std::string& type = atom.get_type();
            return !type.empty() && type[0] == 'H';
        });
    }

    // Enhanced properties calculation
    double get_total_mass() const {
        return ResidueValidator::calculate_total_mass(*this);
    }

    double get_total_charge() const {
        return ResidueValidator::calculate_total_charge(*this);
    }

    bool has_reasonable_geometry() const {
        return ResidueValidator::has_reasonable_geometry(*this);
    }

    bool has_reasonable_charge(double tolerance = 0.01) const {
        return ResidueValidator::is_charge_reasonable(*this, tolerance);
    }

    bool has_reasonable_mass() const {
        return ResidueValidator::is_mass_reasonable(*this);
    }

    // Residue type checking
    bool is_protein_residue() const {
        return ResidueValidator::is_protein_residue(resname);
    }

    bool is_nucleic_acid_residue() const {
        return ResidueValidator::is_nucleic_acid_residue(resname);
    }

    bool has_complete_backbone() const {
        if (is_protein_residue()) {
            return ResidueValidator::has_backbone_atoms(*this);
        } else if (is_nucleic_acid_residue()) {
            return ResidueValidator::has_nucleic_backbone_atoms(*this);
        }
        return true; // Unknown residue types are assumed complete
    }

    // Comprehensive validation
    bool is_valid() const override {
        return ResidueValidator::is_valid(*this);
    }

    std::string get_validation_report() const {
        return ResidueValidator::get_validation_report(*this);
    }

    // Enhanced atom access with bounds checking
    std::shared_ptr<atom::Atom> get_atom_by_index(size_t index) const {
        return (index < atoms.size()) ? atoms[index] : nullptr;
    }

    // Distance calculations
    double distance_to(const Residue& other) const {
        // Distance between centers of mass
        const auto& com1 = get_center_of_mass();
        const auto& com2 = other.get_center_of_mass();
        
        double dx = com1[0] - com2[0];
        double dy = com1[1] - com2[1];
        double dz = com1[2] - com2[2];
        
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    double min_distance_to(const Residue& other) const {
        // Minimum distance between any atoms
        double min_dist = std::numeric_limits<double>::max();
        
        for (const auto& atom1 : atoms) {
            if (!atom1) continue;
            for (const auto& atom2 : other.atoms) {
                if (!atom2) continue;
                double dist = atom1->distance_to(*atom2);
                min_dist = std::min(min_dist, dist);
            }
        }
        
        return (min_dist == std::numeric_limits<double>::max()) ? 0.0 : min_dist;
    }

    // Comparison operators for sorting/searching
    bool operator<(const Residue& other) const {
        if (segid != other.segid) return segid < other.segid;
        if (chain != other.chain) return chain < other.chain;
        return ires < other.ires;
    }

    bool operator==(const Residue& other) const {
        return resname == other.resname && 
               ires == other.ires && 
               segid == other.segid &&
               chain == other.chain &&
               inscode == other.inscode;
    }

    bool operator!=(const Residue& other) const {
        return !(*this == other);
    }

    // Hash function for use in unordered containers
    struct Hash {
        std::size_t operator()(const Residue& residue) const {
            std::size_t seed = 0;
            common::utils::hash_combine(seed, residue.resname);
            common::utils::hash_combine(seed, residue.ires);
            common::utils::hash_combine(seed, residue.segid);
            common::utils::hash_combine(seed, residue.chain);
            common::utils::hash_combine(seed, residue.inscode);
            return seed;
        }
    };

    // Clone method
    std::unique_ptr<Residue> clone() const {
        return std::make_unique<Residue>(*this);
    }

    // String representation
    std::string to_string() const {
        std::stringstream ss;
        ss << resname << " " << ires;
        if (inscode != ' ') ss << inscode;
        if (chain != ' ') ss << " (chain " << chain << ")";
        ss << " [" << segid << "]";
        ss << " (" << atoms.size() << " atoms)";
        return ss.str();
    }

    // PDB format output for all atoms
    std::string to_pdb_string() const {
        std::stringstream ss;
        for (const auto& atom : atoms) {
            if (atom) {
                ss << atom->get_pdb_record() << "\n";
            }
        }
        return ss.str();
    }

    // Summary statistics
    struct Statistics {
        size_t total_atoms;
        size_t heavy_atoms;
        size_t hydrogen_atoms;
        double total_mass;
        double total_charge;
        std::array<double, 3> center_of_mass;
        bool is_complete;
        bool is_valid;
    };

    Statistics get_statistics() const {
        Statistics stats;
        stats.total_atoms = atoms.size();
        stats.heavy_atoms = select_heavy_atoms().size();
        stats.hydrogen_atoms = select_hydrogen_atoms().size();
        stats.total_mass = get_total_mass();
        stats.total_charge = get_total_charge();
        stats.center_of_mass = get_center_of_mass();
        stats.is_complete = has_complete_backbone();
        stats.is_valid = is_valid();
        return stats;
    }

    // Utility methods for CHARMM compatibility
    std::string get_unique_identifier() const {
        std::stringstream ss;
        ss << segid << ":" << resname << ":" << ires;
        if (inscode != ' ') ss << inscode;
        if (chain != ' ') ss << ":" << chain;
        return ss.str();
    }

private:
};

// Utility functions for residue collections
namespace utils {

/**
 * @brief Find residues by name in a collection
 */
template<typename Container>
auto find_residues_by_name(const Container& residues, const std::string& resname) {
    std::vector<typename Container::value_type> result;
    std::copy_if(residues.begin(), residues.end(), std::back_inserter(result),
                [&resname](const auto& residue) { 
                    return residue.get_resname() == resname; 
                });
    return result;
}

/**
 * @brief Find residues by chain in a collection
 */
template<typename Container>
auto find_residues_by_chain(const Container& residues, char chain) {
    std::vector<typename Container::value_type> result;
    std::copy_if(residues.begin(), residues.end(), std::back_inserter(result),
                [chain](const auto& residue) { 
                    return residue.get_chain() == chain; 
                });
    return result;
}

/**
 * @brief Find residues in range
 */
template<typename Container>
auto find_residues_in_range(const Container& residues, int start, int end) {
    std::vector<typename Container::value_type> result;
    std::copy_if(residues.begin(), residues.end(), std::back_inserter(result),
                [start, end](const auto& residue) { 
                    int ires = residue.get_ires();
                    return ires >= start && ires <= end; 
                });
    return result;
}

} // namespace utils

} // namespace residue
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_RESIDUE_MAIN_HPP 