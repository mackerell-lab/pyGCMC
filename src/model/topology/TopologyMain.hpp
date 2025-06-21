#pragma once

#ifndef PYGCMC_MODEL_TOPOLOGY_MAIN_HPP
#define PYGCMC_MODEL_TOPOLOGY_MAIN_HPP

#include "../common/ModelInterface.hpp"
#include "../common/ModelConstants.hpp"
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <memory>
#include <array>
#include <optional>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <tuple>

namespace pygcmc {
namespace model {
namespace topology {

/**
 * @brief Represents an atom in the topology
 */
struct TopologyAtom {
    int id;                 ///< Atom ID (0-based)
    std::string name;       ///< Atom name (e.g., "CA", "N", "O")
    std::string type;       ///< Atom type (e.g., "CT1", "NH1")
    double charge;          ///< Atomic charge
    double mass;            ///< Atomic mass
    int residue_id;        ///< ID of the residue this atom belongs to
    int segment_id;        ///< ID of the segment this atom belongs to
    // B-state parameters for free energy calculations
    std::string typeB;
    double chargeB = 0.0;
    double massB = 0.0;
    bool has_b_state = false;
};

/**
 * @brief Represents a residue in the topology
 */
struct TopologyResidue {
    int id;                 ///< Residue ID (0-based)
    std::string name;       ///< Residue name (e.g., "ALA", "GLY")
    int number;             ///< Residue number from PDB/PSF
    std::vector<int> atoms; ///< Indices of atoms in this residue
    std::string segment;    ///< Segment name this residue belongs to
};

/**
 * @brief Represents a segment in the topology
 */
struct TopologySegment {
    int id;                 ///< Segment ID (0-based)
    std::string name;       ///< Segment name (e.g., "PROT", "MEMB")
    std::vector<int> residues; ///< Indices of residues in this segment
};

/**
 * @brief Represents a bond between two atoms
 */
struct TopologyBond {
    int atom1;              ///< Index of first atom
    int atom2;              ///< Index of second atom
    double length;          ///< Equilibrium bond length (optional)
    double force_constant;  ///< Bond force constant (optional)
    int function_type = 1;    // Default to GROMACS function type 1
};

/**
 * @brief Represents an angle between three atoms
 */
struct TopologyAngle {
    int atom1;              ///< Index of first atom
    int atom2;              ///< Index of central atom
    int atom3;              ///< Index of third atom
    double angle;           ///< Equilibrium angle in degrees (optional)
    double force_constant;  ///< Angle force constant (optional)
    int function_type = 1;   // Default to GROMACS function type 1
    // Urey-Bradley terms
    double ub_length = 0.0;
    double ub_constant = 0.0;
    bool has_ub = false;
};

/**
 * @brief Represents a dihedral angle between four atoms
 */
struct TopologyDihedral {
    int atom1;              ///< Index of first atom
    int atom2;              ///< Index of second atom
    int atom3;              ///< Index of third atom
    int atom4;              ///< Index of fourth atom
    int multiplicity;       ///< Dihedral multiplicity
    double angle;           ///< Equilibrium angle in degrees
    double force_constant;  ///< Dihedral force constant
    bool improper;          ///< Whether this is an improper dihedral
    int function_type = 1;   // Default to GROMACS function type 1
};

/**
 * @brief Represents a hydrogen bond donor
 */
struct TopologyDonor {
    int donor_atom;     ///< Index of the donor atom
    int hydrogen_atom;  ///< Index of the hydrogen atom
};

/**
 * @brief Represents a hydrogen bond acceptor
 */
struct TopologyAcceptor {
    int acceptor_atom;  ///< Index of the acceptor atom
};

/**
 * @brief Represents a group in the topology
 */
struct TopologyGroup {
    int id;                    ///< Group ID
    std::vector<int> atoms;    ///< Indices of atoms in this group
    std::string type;          ///< Group type (e.g., "WATER")
};

/**
 * @brief Represents a CMAP (correction map) term
 */
struct TopologyCmap {
    std::array<int, 8> atoms;  ///< 8 atoms involved in CHARMM CMAP term
    int function_type = 1;     ///< CMAP function type (default: 1)
};

/**
 * @brief Main topology class that holds all molecular topology information
 * This class maintains backward compatibility with the original topology.hpp
 */
class Topology : public common::IValidatable {
public:
    Topology() = default;
    ~Topology() = default;

    // IValidatable interface
    bool is_valid() const override {
        // Check basic consistency
        for (const auto& atom : atoms_) {
            if (atom.residue_id < 0 || atom.residue_id >= static_cast<int>(residues_.size()) ||
                atom.segment_id < 0 || atom.segment_id >= static_cast<int>(segments_.size())) {
                return false;
            }
        }
        return true;
    }

    // Add elements to topology
    inline int add_atom(const std::string& name, const std::string& type, double charge, double mass,
                const std::string& residue_name, int residue_number, const std::string& segment_name) {
        try {
            // First ensure we have the segment
            int segment_id;
            auto segment_it = segment_map_.find(segment_name);
            if (segment_it == segment_map_.end()) {
                // Create new segment
                TopologySegment segment;
                segment.id = static_cast<int>(segments_.size());
                segment.name = segment_name;
                segment_id = segment.id;
                segments_.push_back(segment);
                segment_map_[segment_name] = segment_id;
            } else {
                segment_id = segment_it->second;
            }

            // Then ensure we have the residue
            int residue_id;
            auto residue_key = std::make_tuple(residue_name, residue_number, segment_name);
            auto residue_it = residue_map_.find(residue_key);
            if (residue_it == residue_map_.end()) {
                // Create new residue
                TopologyResidue residue;
                residue.id = static_cast<int>(residues_.size());
                residue.name = residue_name;
                residue.number = residue_number;
                residue.segment = segment_name;
                residue_id = residue.id;
                residues_.push_back(residue);
                residue_map_[residue_key] = residue_id;
                
                // Add residue to its segment
                if (segment_id >= 0 && segment_id < static_cast<int>(segments_.size())) {
                    segments_[segment_id].residues.push_back(residue_id);
                } else {
                    throw std::runtime_error("Invalid segment ID");
                }
            } else {
                residue_id = residue_it->second;
            }

            // Create and add the atom
            TopologyAtom atom;
            atom.id = static_cast<int>(atoms_.size());
            atom.name = name;
            atom.type = type;
            atom.charge = charge;
            atom.mass = mass;
            atom.residue_id = residue_id;
            atom.segment_id = segment_id;

            // Add atom to the vectors and maps
            int atom_id = static_cast<int>(atoms_.size());
            atoms_.push_back(atom);
            atom_map_[std::make_tuple(residue_name, residue_number, segment_name, name)] = atom_id;

            // Add atom to its residue
            if (residue_id >= 0 && residue_id < static_cast<int>(residues_.size())) {
                residues_[residue_id].atoms.push_back(atom_id);
            } else {
                throw std::runtime_error("Invalid residue ID");
            }

            return atom_id;
        } catch (const std::exception& e) {
            std::cerr << "Error in add_atom: " << e.what() << std::endl;
            throw;
        }
    }

    inline void add_bond(int atom1, int atom2, double length = 0.0, double force_constant = 0.0, int function_type = 1) {
        if (!has_atom(atom1) || !has_atom(atom2)) {
            throw std::out_of_range("Invalid atom indices in add_bond");
        }

        TopologyBond bond;
        bond.atom1 = atom1;
        bond.atom2 = atom2;
        bond.length = length;
        bond.force_constant = force_constant;
        bond.function_type = function_type;
        bonds_.push_back(bond);
    }

    inline void add_angle(int atom1, int atom2, int atom3, double angle = 0.0, double force_constant = 0.0, int function_type = 1) {
        if (!has_atom(atom1) || !has_atom(atom2) || !has_atom(atom3)) {
            throw std::out_of_range("Invalid atom indices in add_angle");
        }

        TopologyAngle ang;
        ang.atom1 = atom1;
        ang.atom2 = atom2;
        ang.atom3 = atom3;
        ang.angle = angle;
        ang.force_constant = force_constant;
        ang.function_type = function_type;
        ang.ub_length = 0.0;
        ang.ub_constant = 0.0;
        angles_.push_back(ang);
    }

    inline void add_dihedral(int atom1, int atom2, int atom3, int atom4, int multiplicity = 1,
                     double angle = 0.0, double force_constant = 0.0, bool improper = false, int function_type = 1) {
        if (!has_atom(atom1) || !has_atom(atom2) || !has_atom(atom3) || !has_atom(atom4)) {
            throw std::out_of_range("Invalid atom indices in add_dihedral");
        }

        TopologyDihedral dihedral;
        dihedral.atom1 = atom1;
        dihedral.atom2 = atom2;
        dihedral.atom3 = atom3;
        dihedral.atom4 = atom4;
        dihedral.multiplicity = multiplicity;
        dihedral.angle = angle;
        dihedral.force_constant = force_constant;
        dihedral.improper = improper;
        dihedral.function_type = function_type;
        dihedrals_.push_back(dihedral);
    }

    inline int add_residue(const std::string& name, int number, const std::string& segment) {
        // First ensure we have the segment
        int segment_id;
        auto segment_it = segment_map_.find(segment);
        if (segment_it == segment_map_.end()) {
            segment_id = add_segment(segment);
        } else {
            segment_id = segment_it->second;
        }

        // Create and add the residue
        TopologyResidue residue;
        residue.id = residues_.size();
        residue.name = name;
        residue.number = number;
        residue.segment = segment;

        // Add residue to the vectors and maps
        int residue_id = residues_.size();
        residues_.push_back(residue);
        residue_map_[std::make_tuple(name, number, segment)] = residue_id;

        // Add residue to its segment
        segments_[segment_id].residues.push_back(residue_id);

        return residue_id;
    }

    inline int add_segment(const std::string& name) {
        TopologySegment segment;
        segment.id = segments_.size();
        segment.name = name;

        // Add segment to the vectors and maps
        int segment_id = segments_.size();
        segments_.push_back(segment);
        segment_map_[name] = segment_id;

        return segment_id;
    }

    // Getters
    inline const TopologyAtom& get_atom(int index) const {
        if (!has_atom(index)) {
            throw std::out_of_range("Invalid atom index");
        }
        return atoms_[index];
    }

    inline const TopologyResidue& get_residue(int index) const {
        if (!has_residue(index)) {
            throw std::out_of_range("Invalid residue index");
        }
        return residues_[index];
    }

    inline const TopologySegment& get_segment(int index) const {
        if (!has_segment(index)) {
            throw std::out_of_range("Invalid segment index");
        }
        return segments_[index];
    }

    inline const std::vector<TopologyBond>& get_bonds() const { return bonds_; }
    inline const std::vector<TopologyAngle>& get_angles() const { return angles_; }
    inline const std::vector<TopologyDihedral>& get_dihedrals() const { return dihedrals_; }

    // Utility functions
    inline bool has_atom(int index) const {
        return index >= 0 && index < static_cast<int>(atoms_.size());
    }

    inline bool has_residue(int index) const {
        return index >= 0 && index < static_cast<int>(residues_.size());
    }

    inline bool has_segment(int index) const {
        return index >= 0 && index < static_cast<int>(segments_.size());
    }

    inline int get_num_atoms() const { return atoms_.size(); }
    inline int get_num_residues() const { return residues_.size(); }
    inline int get_num_segments() const { return segments_.size(); }
    
    // Find elements
    inline std::optional<int> find_atom(const std::string& residue_name, int residue_number,
                                const std::string& atom_name) const {
        // Try to find the atom in any segment
        for (const auto& segment : segments_) {
            auto it = atom_map_.find(std::make_tuple(residue_name, residue_number, segment.name, atom_name));
            if (it != atom_map_.end()) {
                return it->second;
            }
        }
        return std::nullopt;
    }

    inline std::optional<int> find_residue(const std::string& name, int number) const {
        // Try to find the residue in any segment
        for (const auto& segment : segments_) {
            auto it = residue_map_.find(std::make_tuple(name, number, segment.name));
            if (it != residue_map_.end()) {
                return it->second;
            }
        }
        return std::nullopt;
    }

    inline std::optional<int> find_segment(const std::string& name) const {
        auto it = segment_map_.find(name);
        if (it != segment_map_.end()) {
            return it->second;
        }
        return std::nullopt;
    }

    // Additional methods for PSF sections
    void add_title(const std::string& title) { titles_.push_back(title); }
    
    void add_improper(int atom1, int atom2, int atom3, int atom4,
                     double angle = 0.0, double force_constant = 0.0) {
        add_dihedral(atom1, atom2, atom3, atom4, 0, angle, force_constant, true);
    }

    void add_donor(int donor, int hydrogen) {
        if (!has_atom(donor) || !has_atom(hydrogen)) {
            throw std::out_of_range("Invalid atom indices in add_donor");
        }
        TopologyDonor d;
        d.donor_atom = donor;
        d.hydrogen_atom = hydrogen;
        donors_.push_back(d);
    }

    void add_acceptor(int acceptor) {
        if (!has_atom(acceptor)) {
            throw std::out_of_range("Invalid atom index in add_acceptor");
        }
        TopologyAcceptor a;
        a.acceptor_atom = acceptor;
        acceptors_.push_back(a);
    }

    void add_nonbonded_exclusion(int atom1, int atom2) {
        if (!has_atom(atom1) || !has_atom(atom2)) {
            throw std::out_of_range("Invalid atom indices in add_nonbonded_exclusion");
        }
        exclusions_[atom1].insert(atom2);
        exclusions_[atom2].insert(atom1);
    }

    void add_group(int id, const std::vector<int>& atoms, const std::string& type = "") {
        for (int atom : atoms) {
            if (!has_atom(atom)) {
                throw std::out_of_range("Invalid atom index in add_group");
            }
        }
        TopologyGroup group;
        group.id = id;
        group.atoms = atoms;
        group.type = type;
        groups_.push_back(group);
    }

    void add_cmap(const std::array<int, 8>& atoms) {
        for (int atom : atoms) {
            if (!has_atom(atom)) {
                throw std::out_of_range("Invalid atom index in add_cmap");
            }
        }
        TopologyCmap cmap;
        cmap.atoms = atoms;
        cmap.function_type = 1;
        cmaps_.push_back(cmap);
    }

    void add_cmap(const std::array<int, 5>& atoms, int function_type = 1) {
        for (int atom : atoms) {
            if (!has_atom(atom)) {
                throw std::out_of_range("Invalid atom index in add_cmap");
            }
        }
        // Convert 5-atom GROMACS format to 8-atom CHARMM format
        std::array<int, 8> charmm_atoms;
        for (int i = 0; i < 5; ++i) {
            charmm_atoms[i] = atoms[i];
        }
        for (int i = 5; i < 8; ++i) {
            charmm_atoms[i] = -1;
        }
        
        TopologyCmap cmap;
        cmap.atoms = charmm_atoms;
        cmap.function_type = function_type;
        cmaps_.push_back(cmap);
    }

    // Getters for additional sections
    const std::vector<std::string>& get_titles() const { return titles_; }
    const std::vector<TopologyDonor>& get_donors() const { return donors_; }
    const std::vector<TopologyAcceptor>& get_acceptors() const { return acceptors_; }
    const std::vector<TopologyGroup>& get_groups() const { return groups_; }
    const std::vector<TopologyCmap>& get_cmaps() const { return cmaps_; }
    const std::map<int, std::set<int>>& get_exclusions() const { return exclusions_; }

    // Count methods
    inline size_t get_num_bonds() const { return bonds_.size(); }
    inline size_t get_num_angles() const { return angles_.size(); }
    inline size_t get_num_dihedrals() const { 
        size_t count = 0;
        for (const auto& dih : dihedrals_) {
            if (!dih.improper) count++;
        }
        return count;
    }
    inline size_t get_num_impropers() const { 
        size_t count = 0;
        for (const auto& dih : dihedrals_) {
            if (dih.improper) count++;
        }
        return count;
    }
    inline size_t get_num_donors() const { return donors_.size(); }
    inline size_t get_num_acceptors() const { return acceptors_.size(); }
    inline size_t get_num_cmaps() const { return cmaps_.size(); }
    inline size_t get_num_groups() const { return groups_.size(); }

    // Check methods for different topology elements
    inline bool has_donor(int donor_atom) const {
        return std::any_of(donors_.begin(), donors_.end(),
            [donor_atom](const TopologyDonor& d) { return d.donor_atom == donor_atom; });
    }
    
    inline bool has_acceptor(int acceptor_atom) const {
        return std::any_of(acceptors_.begin(), acceptors_.end(),
            [acceptor_atom](const TopologyAcceptor& a) { return a.acceptor_atom == acceptor_atom; });
    }
    
    inline bool has_cmap() const {
        return !cmaps_.empty();
    }
    
    inline bool has_group(int group_id) const {
        return std::any_of(groups_.begin(), groups_.end(),
            [group_id](const TopologyGroup& g) { return g.id == group_id; });
    }
    
    inline const TopologyGroup& get_group(int index) const {
        if (index < 0 || index >= static_cast<int>(groups_.size())) {
            throw std::out_of_range("Invalid group index");
        }
        return groups_[index];
    }

    // Existence check methods
    inline bool has_bond(int atom1, int atom2) const {
        for (const auto& bond : bonds_) {
            if ((bond.atom1 == atom1 && bond.atom2 == atom2) ||
                (bond.atom1 == atom2 && bond.atom2 == atom1)) {
                return true;
            }
        }
        return false;
    }

    inline bool has_angle(int atom1, int atom2, int atom3) const {
        for (const auto& angle : angles_) {
            if ((angle.atom1 == atom1 && angle.atom2 == atom2 && angle.atom3 == atom3) ||
                (angle.atom1 == atom3 && angle.atom2 == atom2 && angle.atom3 == atom1)) {
                return true;
            }
        }
        return false;
    }

    inline bool has_dihedral(int atom1, int atom2, int atom3, int atom4) const {
        for (const auto& dihedral : dihedrals_) {
            if (!dihedral.improper &&
                ((dihedral.atom1 == atom1 && dihedral.atom2 == atom2 && 
                  dihedral.atom3 == atom3 && dihedral.atom4 == atom4) ||
                 (dihedral.atom1 == atom4 && dihedral.atom2 == atom3 && 
                  dihedral.atom3 == atom2 && dihedral.atom4 == atom1))) {
                return true;
            }
        }
        return false;
    }

    inline bool has_improper(int atom1, int atom2, int atom3, int atom4) const {
        for (const auto& dihedral : dihedrals_) {
            if (dihedral.improper &&
                dihedral.atom1 == atom1 && dihedral.atom2 == atom2 && 
                dihedral.atom3 == atom3 && dihedral.atom4 == atom4) {
                return true;
            }
        }
        return false;
    }

    void reserve_atoms(size_t n) {
        atoms_.reserve(n);
    }

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

    // Lookup maps
    std::unordered_map<std::string, int> segment_map_;
    std::map<std::tuple<std::string, int, std::string>, int> residue_map_;
    std::map<std::tuple<std::string, int, std::string, std::string>, int> atom_map_;

    std::vector<std::string> titles_;
};

} // namespace topology

// Backward compatibility: provide the Topology classes in the model namespace
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

#endif // PYGCMC_MODEL_TOPOLOGY_MAIN_HPP 