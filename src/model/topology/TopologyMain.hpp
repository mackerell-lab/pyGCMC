#pragma once

#include "TopologyStructures.hpp"
#include "TopologyOperations.hpp"
#include "TopologyQueries.hpp"

namespace pygcmc {
namespace model {
namespace topology {

/**
 * @brief Main topology class following system module design patterns
 * AI-friendly: clear separation of structures, operations, and queries
 */
class Topology {
public:
    Topology() = default;
    ~Topology() = default;

    // Add elements to topology - delegate to TopologyOperations
    inline int add_atom(const std::string& name, const std::string& type, double charge, double mass,
                const std::string& residue_name, int residue_number, const std::string& segment_name) {
        return TopologyOperations::addAtom(atoms_, residues_, segments_, segment_map_, 
                                          residue_map_, atom_map_, name, type, charge, mass, 
                                          residue_name, residue_number, segment_name);
    }

    inline void add_bond(int atom1, int atom2, double length = 0.0, double force_constant = 0.0, int function_type = 1) {
        TopologyOperations::addBond(bonds_, atoms_, atom1, atom2, length, force_constant, function_type);
    }

    inline void add_angle(int atom1, int atom2, int atom3, double angle = 0.0, double force_constant = 0.0, int function_type = 1) {
        TopologyOperations::addAngle(angles_, atoms_, atom1, atom2, atom3, angle, force_constant, function_type);
    }

    inline void add_dihedral(int atom1, int atom2, int atom3, int atom4, int multiplicity = 1,
                     double angle = 0.0, double force_constant = 0.0, bool improper = false, int function_type = 1) {
        TopologyOperations::addDihedral(dihedrals_, atoms_, atom1, atom2, atom3, atom4, 
                                       multiplicity, angle, force_constant, improper, function_type);
    }

    inline void add_improper(int atom1, int atom2, int atom3, int atom4, double angle = 0.0, double force_constant = 0.0) {
        TopologyOperations::addDihedral(dihedrals_, atoms_, atom1, atom2, atom3, atom4, 
                                       0, angle, force_constant, true, 1);
    }

    inline int add_residue(const std::string& name, int number, const std::string& segment) {
        return TopologyOperations::addResidue(residues_, segments_, segment_map_, residue_map_, name, number, segment);
    }

    inline int add_segment(const std::string& name) {
        return TopologyOperations::addSegment(segments_, segment_map_, name);
    }

    inline void add_donor(int donor, int hydrogen) {
        TopologyOperations::addDonor(donors_, atoms_, donor, hydrogen);
    }

    inline void add_acceptor(int acceptor) {
        TopologyOperations::addAcceptor(acceptors_, atoms_, acceptor);
    }

    inline void add_nonbonded_exclusion(int atom1, int atom2) {
        TopologyOperations::addExclusion(exclusions_, atoms_, atom1, atom2);
    }

    inline void add_group(int id, const std::vector<int>& atoms, const std::string& type = "") {
        TopologyOperations::addGroup(groups_, atoms_, id, atoms, type);
    }

    inline void add_cmap(const std::array<int, 8>& atoms) {
        TopologyOperations::addCmap(cmaps_, atoms_, atoms);
    }

    inline void add_cmap(const std::array<int, 5>& atoms, int function_type = 1) {
        TopologyOperations::addCmap(cmaps_, atoms_, atoms, function_type);
    }

    // Getters
    inline const TopologyAtom& get_atom(int index) const {
        if (!TopologyQueries::hasAtom(atoms_, index)) throw std::out_of_range("Invalid atom index");
        return atoms_[index];
    }

    inline const TopologyResidue& get_residue(int index) const {
        if (!TopologyQueries::hasResidue(residues_, index)) throw std::out_of_range("Invalid residue index");
        return residues_[index];
    }

    inline const TopologySegment& get_segment(int index) const {
        if (!TopologyQueries::hasSegment(segments_, index)) throw std::out_of_range("Invalid segment index");
        return segments_[index];
    }

    // Collection getters
    inline const std::vector<TopologyBond>& get_bonds() const { return bonds_; }
    inline const std::vector<TopologyAngle>& get_angles() const { return angles_; }
    inline const std::vector<TopologyDihedral>& get_dihedrals() const { return dihedrals_; }
    inline const std::vector<TopologyDonor>& get_donors() const { return donors_; }
    inline const std::vector<TopologyAcceptor>& get_acceptors() const { return acceptors_; }
    inline const std::vector<TopologyGroup>& get_groups() const { return groups_; }
    inline const std::vector<TopologyCmap>& get_cmaps() const { return cmaps_; }
    inline const std::map<int, std::set<int>>& get_exclusions() const { return exclusions_; }

    // Utility functions - delegate to TopologyQueries
    inline bool has_atom(int index) const { return TopologyQueries::hasAtom(atoms_, index); }
    inline bool has_residue(int index) const { return TopologyQueries::hasResidue(residues_, index); }
    inline bool has_segment(int index) const { return TopologyQueries::hasSegment(segments_, index); }

    // Count functions
    inline int get_num_atoms() const { return atoms_.size(); }
    inline int get_num_residues() const { return residues_.size(); }
    inline int get_num_segments() const { return segments_.size(); }
    inline size_t get_num_bonds() const { return bonds_.size(); }
    inline size_t get_num_angles() const { return angles_.size(); }
    inline size_t get_num_dihedrals() const { return TopologyQueries::countDihedrals(dihedrals_); }
    inline size_t get_num_impropers() const { return TopologyQueries::countImpropers(dihedrals_); }
    inline size_t get_num_donors() const { return donors_.size(); }
    inline size_t get_num_acceptors() const { return acceptors_.size(); }
    inline size_t get_num_cmaps() const { return cmaps_.size(); }
    inline size_t get_num_groups() const { return groups_.size(); }

    // Find operations - delegate to TopologyQueries
    inline std::optional<int> find_atom(const std::string& residue_name, int residue_number, const std::string& atom_name) const {
        return TopologyQueries::findAtom(segments_, atom_map_, residue_name, residue_number, atom_name);
    }

    inline std::optional<int> find_residue(const std::string& name, int number) const {
        return TopologyQueries::findResidue(segments_, residue_map_, name, number);
    }

    inline std::optional<int> find_segment(const std::string& name) const {
        return TopologyQueries::findSegment(segment_map_, name);
    }

    // Check operations - delegate to TopologyQueries
    inline bool has_bond(int atom1, int atom2) const {
        return TopologyQueries::hasBond(bonds_, atom1, atom2);
    }

    inline bool has_angle(int atom1, int atom2, int atom3) const {
        return TopologyQueries::hasAngle(angles_, atom1, atom2, atom3);
    }

    inline bool has_dihedral(int atom1, int atom2, int atom3, int atom4) const {
        return TopologyQueries::hasDihedral(dihedrals_, atom1, atom2, atom3, atom4);
    }

    inline bool has_improper(int atom1, int atom2, int atom3, int atom4) const {
        return TopologyQueries::hasImproper(dihedrals_, atom1, atom2, atom3, atom4);
    }

    inline bool has_donor(int donor_atom) const {
        return TopologyQueries::hasDonor(donors_, donor_atom);
    }

    inline bool has_donor(int donor_atom, int hydrogen_atom) const {
        return TopologyQueries::hasDonor(donors_, donor_atom, hydrogen_atom);
    }

    inline bool has_acceptor(int acceptor_atom) const {
        return TopologyQueries::hasAcceptor(acceptors_, acceptor_atom);
    }

    inline bool has_cmap() const {
        return TopologyQueries::hasCmap(cmaps_);
    }

    inline bool has_cmap(const std::vector<int>& atoms) const {
        return TopologyQueries::hasCmap(cmaps_, atoms);
    }

    inline bool has_group(int group_id) const {
        return TopologyQueries::hasGroup(groups_, group_id);
    }

    inline const TopologyGroup& get_group(int index) const {
        if (index < 0 || index >= static_cast<int>(groups_.size())) {
            throw std::out_of_range("Invalid group index");
        }
        return groups_[index];
    }

    // Additional compatibility methods
    inline void add_title(const std::string& title) { titles_.push_back(title); }
    inline const std::vector<std::string>& get_titles() const { return titles_; }
    inline void reserve_atoms(size_t n) { atoms_.reserve(n); }

private:
    std::vector<TopologyAtom> atoms_;
    std::vector<TopologyResidue> residues_;
    std::vector<TopologySegment> segments_;
    std::vector<TopologyBond> bonds_;
    std::vector<TopologyAngle> angles_;
    std::vector<TopologyDihedral> dihedrals_;
    std::vector<TopologyDonor> donors_;
    std::vector<TopologyAcceptor> acceptors_;
    std::map<int, std::set<int>> exclusions_;
    std::vector<TopologyGroup> groups_;
    std::vector<TopologyCmap> cmaps_;

    // Maps for fast lookup
    std::unordered_map<std::string, int> segment_map_;
    std::map<std::tuple<std::string, int, std::string>, int> residue_map_;
    std::map<std::tuple<std::string, int, std::string, std::string>, int> atom_map_;

    std::vector<std::string> titles_;
};

} // namespace topology

// Backward compatibility type aliases
using Topology = topology::Topology;
using TopologyAtom = topology::TopologyAtom;
using TopologyResidue = topology::TopologyResidue;
using TopologySegment = topology::TopologySegment;
using TopologyBond = topology::TopologyBond;
using TopologyAngle = topology::TopologyAngle;
using TopologyDihedral = topology::TopologyDihedral;
using TopologyDonor = topology::TopologyDonor;
using TopologyAcceptor = topology::TopologyAcceptor;
using TopologyGroup = topology::TopologyGroup;
using TopologyCmap = topology::TopologyCmap;

} // namespace model
} // namespace pygcmc