#pragma once

#ifndef PYGCMC_MODEL_TOPOLOGY_UTILS_HPP
#define PYGCMC_MODEL_TOPOLOGY_UTILS_HPP

#include "TopologyCore.hpp"
#include "TopologyStats.hpp"
#include "TopologyValidation.hpp"
#include "TopologySearch.hpp"
#include "TopologyAnalysis.hpp"
#include <algorithm>

namespace pygcmc {
namespace model {
namespace topology {

/**
 * @brief Main topology utilities class combining all functionality
 */
class TopologyUtils {
public:
    TopologyUtils(const TopologyStorage& storage) 
        : storage_(storage),
          stats_generator_(storage),
          validator_(storage),
          searcher_(storage),
          analyzer_(storage) {}

    // Delegate to specialized components
    TopologyStats get_statistics() const { return stats_generator_.get_statistics(); }
    std::string get_summary() const { return stats_generator_.get_summary(); }
    
    bool validate_topology() const { return validator_.validate_topology(); }
    bool has_atom(int index) const { return validator_.has_atom(index); }
    bool has_residue(int index) const { return validator_.has_residue(index); }
    bool has_segment(int index) const { return validator_.has_segment(index); }
    
    std::vector<int> get_residue_atoms(int residue_id) const { return searcher_.get_residue_atoms(residue_id); }
    std::vector<int> get_segment_residues(int segment_id) const { return searcher_.get_segment_residues(segment_id); }
    std::vector<int> get_bonded_atoms(int atom_id) const { return searcher_.get_bonded_atoms(atom_id); }
    std::optional<int> find_atom(const std::string& residue_name, int residue_number, const std::string& atom_name) const {
        return searcher_.find_atom(residue_name, residue_number, atom_name);
    }
    std::optional<int> find_residue(const std::string& name, int number) const { return searcher_.find_residue(name, number); }
    std::optional<int> find_segment(const std::string& name) const { return searcher_.find_segment(name); }
    bool are_bonded(int atom1, int atom2) const { return searcher_.are_bonded(atom1, atom2); }
    const std::vector<std::string>& get_titles() const { return searcher_.get_titles(); }
    
    std::map<std::string, int> count_atom_types() const { return analyzer_.count_atom_types(); }
    std::map<std::string, int> count_residue_types() const { return analyzer_.count_residue_types(); }
    int get_bond_distance(int atom1, int atom2) const { return analyzer_.get_bond_distance(atom1, atom2); }
    std::vector<std::vector<int>> get_connected_components() const { return analyzer_.get_connected_components(); }

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

    // Direct access to specialized components
    const TopologyStatsGenerator& get_stats_generator() const { return stats_generator_; }
    const TopologyValidator& get_validator() const { return validator_; }
    const TopologySearcher& get_searcher() const { return searcher_; }
    const TopologyAnalyzer& get_analyzer() const { return analyzer_; }

private:
    const TopologyStorage& storage_;
    TopologyStatsGenerator stats_generator_;
    TopologyValidator validator_;
    TopologySearcher searcher_;
    TopologyAnalyzer analyzer_;
};

} // namespace topology
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_TOPOLOGY_UTILS_HPP 