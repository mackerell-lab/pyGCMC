#pragma once

#ifndef PYGCMC_MODEL_MOLECULE_MAIN_HPP
#define PYGCMC_MODEL_MOLECULE_MAIN_HPP

#include "MolecularComposite.hpp"
#include "MolecularUtils.hpp"
#include "../common/ModelUtils.hpp"
#include <string>

namespace pygcmc {
namespace model {
namespace molecule {

/**
 * @brief Complete Molecular class with full functionality and backward compatibility
 * This class maintains the same API as the original molecular.hpp while using the refactored structure
 */
class Molecular : public MolecularComposite {
public:
    Molecular() = default;
    ~Molecular() = default;

    // Backward compatibility accessors - use functions instead of references for Python binding compatibility
    std::vector<std::shared_ptr<atom::Atom>>& get_atoms_ref() { return MolecularComposite::get_atoms(); }
    std::vector<std::shared_ptr<residue::Residue>>& get_residues_ref() { return MolecularComposite::get_residues(); }
    
    // Structure info accessors
    std::vector<StructureInfo::TerminalInfo>& get_terminals_ref() { return get_structure_info().terminals; }
    std::map<std::string, std::vector<StructureInfo::SecondaryStructure>>& get_helices_ref() { return get_structure_info().helices; }
    std::map<std::string, std::vector<std::string>>& get_sheets_ref() { return get_structure_info().sheets; }
    std::vector<std::string>& get_ssbonds_ref() { return get_structure_info().ssbonds; }
    std::vector<double>& get_box_dimensions_ref() { return get_structure_info().box_dimensions; }

    // Topology info accessors
    std::vector<topology::TopologyAtom>& get_topology_atoms_ref() { return get_topology_info().atoms; }
    std::vector<topology::TopologyResidue>& get_topology_residues_ref() { return get_topology_info().residues; }
    std::vector<topology::TopologySegment>& get_segments_ref() { return get_topology_info().segments; }
    std::vector<topology::TopologyBond>& get_bonds_ref() { return get_topology_info().bonds; }
    std::vector<topology::TopologyAngle>& get_angles_ref() { return get_topology_info().angles; }
    std::vector<topology::TopologyDihedral>& get_dihedrals_ref() { return get_topology_info().dihedrals; }
    std::vector<topology::TopologyDonor>& get_donors_ref() { return get_topology_info().donors; }
    std::vector<topology::TopologyAcceptor>& get_acceptors_ref() { return get_topology_info().acceptors; }
    std::map<int, std::set<int>>& get_exclusions_ref() { return get_topology_info().exclusions; }
    std::vector<topology::TopologyGroup>& get_groups_ref() { return get_topology_info().groups; }
    std::vector<topology::TopologyCmap>& get_cmaps_ref() { return get_topology_info().cmaps; }
    std::vector<std::string>& get_titles_ref() { return get_topology_info().titles; }
    std::vector<StandardCmap>& get_standard_cmaps_ref() { return get_standard_cmaps(); }

    // Lookup mapping accessors 
    std::unordered_map<std::string, int>& get_segment_map_ref() { return const_cast<std::unordered_map<std::string, int>&>(get_segment_map()); }
    std::map<std::pair<std::string, int>, int>& get_residue_map_ref() { return const_cast<std::map<std::pair<std::string, int>, int>&>(get_residue_map()); }
    std::map<std::tuple<std::string, int, std::string>, int>& get_atom_map_ref() { return const_cast<std::map<std::tuple<std::string, int, std::string>, int>&>(get_atom_map()); }

    // Backward compatibility properties - kept as member variables for legacy code compatibility
    std::vector<std::shared_ptr<atom::Atom>>& atoms = get_atoms_ref();
    std::vector<std::shared_ptr<residue::Residue>>& residues = get_residues_ref();
    std::vector<StructureInfo::TerminalInfo>& terminals = get_terminals_ref();
    std::map<std::string, std::vector<StructureInfo::SecondaryStructure>>& helices = get_helices_ref();
    std::map<std::string, std::vector<std::string>>& sheets = get_sheets_ref();
    std::vector<std::string>& ssbonds = get_ssbonds_ref();
    std::vector<double>& box_dimensions = get_box_dimensions_ref();
    std::vector<topology::TopologyAtom>& topology_atoms = get_topology_atoms_ref();
    std::vector<topology::TopologyResidue>& topology_residues = get_topology_residues_ref();
    std::vector<topology::TopologySegment>& segments = get_segments_ref();
    std::vector<topology::TopologyBond>& bonds = get_bonds_ref();
    std::vector<topology::TopologyAngle>& angles = get_angles_ref();
    std::vector<topology::TopologyDihedral>& dihedrals = get_dihedrals_ref();
    std::vector<topology::TopologyDonor>& donors = get_donors_ref();
    std::vector<topology::TopologyAcceptor>& acceptors = get_acceptors_ref();
    std::map<int, std::set<int>>& exclusions = get_exclusions_ref();
    std::vector<topology::TopologyGroup>& groups = get_groups_ref();
    std::vector<topology::TopologyCmap>& cmaps = get_cmaps_ref();
    std::vector<std::string>& titles = get_titles_ref();
    std::vector<StandardCmap>& standard_cmaps = get_standard_cmaps_ref();
    std::unordered_map<std::string, int>& segment_map = get_segment_map_ref();
    std::map<std::pair<std::string, int>, int>& residue_map = get_residue_map_ref();
    std::map<std::tuple<std::string, int, std::string>, int>& atom_map = get_atom_map_ref();

    // Enhanced selection and analysis methods
    template<typename Predicate>
    std::vector<std::shared_ptr<atom::Atom>> select_atoms(Predicate&& predicate) const {
        return utils::selection::find_atoms_if(*this, std::forward<Predicate>(predicate));
    }

    template<typename Predicate>
    std::vector<std::shared_ptr<residue::Residue>> select_residues(Predicate&& predicate) const {
        return utils::selection::find_residues_if(*this, std::forward<Predicate>(predicate));
    }

    // Convenient selection methods
    std::vector<std::shared_ptr<atom::Atom>> get_atoms_by_segment(const std::string& segment_id) const {
        return utils::selection::find_atoms_by_segment(*this, segment_id);
    }

    std::vector<std::shared_ptr<atom::Atom>> get_atoms_by_resname(const std::string& resname) const {
        return utils::selection::find_atoms_by_resname(*this, resname);
    }

    std::vector<std::shared_ptr<atom::Atom>> get_atoms_by_type(const std::string& atom_type) const {
        return utils::selection::find_atoms_by_type(*this, atom_type);
    }

    std::vector<std::shared_ptr<residue::Residue>> get_residues_by_segment(const std::string& segment_id) const {
        return utils::selection::find_residues_by_segment(*this, segment_id);
    }

    std::vector<std::shared_ptr<residue::Residue>> get_residues_by_chain(char chain_id) const {
        return utils::selection::find_residues_by_chain(*this, chain_id);
    }

    std::vector<std::shared_ptr<residue::Residue>> get_residues_by_name(const std::string& resname) const {
        return utils::selection::find_residues_by_name(*this, resname);
    }

    std::vector<std::shared_ptr<residue::Residue>> get_protein_residues() const {
        return utils::selection::find_protein_residues(*this);
    }

    std::vector<std::shared_ptr<residue::Residue>> get_nucleic_residues() const {
        return utils::selection::find_nucleic_residues(*this);
    }

    std::vector<std::shared_ptr<atom::Atom>> get_heavy_atoms() const {
        return utils::selection::find_heavy_atoms(*this);
    }

    std::vector<std::shared_ptr<atom::Atom>> get_hydrogen_atoms() const {
        return utils::selection::find_hydrogen_atoms(*this);
    }

    // Analysis methods
    double get_total_mass() const {
        return utils::analysis::calculate_total_mass(*this);
    }

    double get_total_charge() const {
        return utils::analysis::calculate_total_charge(*this);
    }

    std::array<double, 3> get_center_of_mass() const {
        return utils::analysis::calculate_center_of_mass(*this);
    }

    std::array<double, 3> get_geometric_center() const {
        return utils::analysis::calculate_geometric_center(*this);
    }

    std::pair<std::array<double, 3>, std::array<double, 3>> get_bounding_box() const {
        return utils::analysis::calculate_bounding_box(*this);
    }

    std::set<char> get_chain_ids() const {
        return utils::analysis::get_chain_ids(*this);
    }

    std::set<std::string> get_segment_ids() const {
        return utils::analysis::get_segment_ids(*this);
    }

    utils::analysis::SystemStatistics get_system_statistics() const {
        return utils::analysis::calculate_system_statistics(*this);
    }

    // Distance utilities
    std::vector<std::shared_ptr<atom::Atom>> find_atoms_within_distance(
        const std::array<double, 3>& point, double max_distance) const {
        return utils::distance::find_atoms_within_distance(*this, point, max_distance);
    }

    std::vector<std::shared_ptr<atom::Atom>> find_atoms_within_distance(
        const atom::Atom& reference_atom, double max_distance) const {
        return utils::distance::find_atoms_within_distance(*this, reference_atom, max_distance);
    }

    // Grouping utilities
    std::map<std::string, std::vector<std::shared_ptr<atom::Atom>>> group_atoms_by_segment() const {
        return utils::grouping::group_atoms_by_segment(*this);
    }

    std::map<std::string, std::vector<std::shared_ptr<atom::Atom>>> group_atoms_by_resname() const {
        return utils::grouping::group_atoms_by_resname(*this);
    }

    std::map<std::string, std::vector<std::shared_ptr<atom::Atom>>> group_atoms_by_type() const {
        return utils::grouping::group_atoms_by_type(*this);
    }

    std::map<std::string, std::vector<std::shared_ptr<residue::Residue>>> group_residues_by_segment() const {
        return utils::grouping::group_residues_by_segment(*this);
    }

    std::map<char, std::vector<std::shared_ptr<residue::Residue>>> group_residues_by_chain() const {
        return utils::grouping::group_residues_by_chain(*this);
    }

    std::map<std::string, std::vector<std::shared_ptr<residue::Residue>>> group_residues_by_name() const {
        return utils::grouping::group_residues_by_name(*this);
    }

    // String representation
    std::string to_string() const {
        std::stringstream ss;
        ss << "Molecular System: " << get_num_atoms() << " atoms, " 
           << get_num_residues() << " residues";
        
        auto segments = get_segment_ids();
        if (!segments.empty()) {
            ss << ", segments: ";
            bool first = true;
            for (const auto& seg : segments) {
                if (!first) ss << ", ";
                ss << seg;
                first = false;
            }
        }
        
        auto chains = get_chain_ids();
        if (chains.size() > 1) {
            ss << ", chains: ";
            bool first = true;
            for (char chain : chains) {
                if (chain != ' ') {
                    if (!first) ss << ", ";
                    ss << chain;
                    first = false;
                }
            }
        }
        
        return ss.str();
    }

    // PDB format output
    std::string to_pdb_string() const {
        std::stringstream ss;
        
        // Write HEADER
        ss << "HEADER    MOLECULAR SYSTEM                        " 
           << std::setfill('0') << std::setw(2) << 1  // day
           << "-" << std::setw(3) << "JAN"             // month  
           << "-" << std::setw(2) << 24                // year
           << "   PYGE\n";  // PDB ID
        
        // Write TITLE
        if (!titles.empty()) {
            for (size_t i = 0; i < titles.size(); ++i) {
                ss << "TITLE    ";
                if (i > 0) ss << std::setw(2) << (i + 1) << " ";
                ss << titles[i] << "\n";
            }
        }
        
        // Write atoms
        for (const auto& atom : atoms) {
            if (atom) {
                ss << atom->get_pdb_record() << "\n";
            }
        }
        
        ss << "END\n";
        return ss.str();
    }

    // Comparison operators
    bool operator==(const Molecular& other) const {
        return get_num_atoms() == other.get_num_atoms() &&
               get_num_residues() == other.get_num_residues() &&
               common::utils::double_equals(get_total_mass(), other.get_total_mass()) &&
               common::utils::double_equals(get_total_charge(), other.get_total_charge());
    }

    bool operator!=(const Molecular& other) const {
        return !(*this == other);
    }

    // Clone method
    std::unique_ptr<Molecular> clone() const {
        auto cloned = std::make_unique<Molecular>();
        
        // Copy atoms
        for (const auto& atom : atoms) {
            if (atom) {
                cloned->add_atom(atom->clone());
            }
        }
        
        // Copy residues
        for (const auto& residue : residues) {
            if (residue) {
                cloned->add_residue(residue->clone());
            }
        }
        
        // Copy structure and topology info
        cloned->get_structure_info() = get_structure_info();
        cloned->get_topology_info() = get_topology_info();
        
        return cloned;
    }
};

/**
 * @brief Validate a complete molecular system
 * @param system The molecular system to validate
 * @return true if the system is valid, false otherwise
 */
inline bool validate_molecular_system(const Molecular& system) {
    // Check basic validity
    if (!system.is_valid()) return false;
    
    // Check all residues
    for (const auto& residue : system.get_residues()) {
        if (!residue || !residue->is_valid()) return false;
    }
    
    // Check all atoms
    for (const auto& atom : system.get_atoms()) {
        if (!atom || !atom->is_valid()) return false;
    }
    
    return true;
}

} // namespace molecule

// Backward compatibility: provide the Molecular class in the model namespace
using Molecular = molecule::Molecular;

// Backward compatibility: provide the StandardCmap struct in the model namespace
using StandardCmap = molecule::StandardCmap;

/**
 * @brief Validate a complete molecular system (backward compatibility)
 */
inline bool validate_molecular_system(const Molecular& system) {
    return molecule::validate_molecular_system(system);
}

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_MOLECULE_MAIN_HPP 