// src/model/topology.hpp

#pragma once

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

namespace pygcmc {
namespace model {

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
    std::array<int, 8> atoms;  ///< 8 atoms involved in CMAP term
};

/**
 * @brief Main topology class that holds all molecular topology information
 */
class Topology {
public:
    Topology() = default;
    ~Topology() = default;

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
                segment.id = static_cast<int>(segments_.size());  // 0-based indexing
                segment.name = segment_name;
                segment_id = segment.id;
                segments_.push_back(segment);
                segment_map_[segment_name] = segment_id;
            } else {
                segment_id = segment_it->second;
            }

            // Then ensure we have the residue
            int residue_id;
            auto residue_key = std::make_pair(residue_name, residue_number);
            auto residue_it = residue_map_.find(residue_key);
            if (residue_it == residue_map_.end()) {
                // Create new residue
                TopologyResidue residue;
                residue.id = static_cast<int>(residues_.size());  // 0-based indexing
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
            atom.id = static_cast<int>(atoms_.size());  // 0-based indexing
            atom.name = name;
            atom.type = type;
            atom.charge = charge;
            atom.mass = mass;
            atom.residue_id = residue_id;
            atom.segment_id = segment_id;

            // Add atom to the vectors and maps
            int atom_id = static_cast<int>(atoms_.size());
            atoms_.push_back(atom);
            atom_map_[std::make_tuple(residue_name, residue_number, name)] = atom_id;

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
        residue.id = residues_.size();  // 0-based indexing
        residue.name = name;
        residue.number = number;
        residue.segment = segment;

        // Add residue to the vectors and maps
        int residue_id = residues_.size();
        residues_.push_back(residue);
        residue_map_[std::make_pair(name, number)] = residue_id;

        // Add residue to its segment
        segments_[segment_id].residues.push_back(residue_id);

        return residue_id;
    }

    inline int add_segment(const std::string& name) {
        TopologySegment segment;
        segment.id = segments_.size();  // 0-based indexing
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
        auto it = atom_map_.find(std::make_tuple(residue_name, residue_number, atom_name));
        if (it != atom_map_.end()) {
            return it->second;
        }
        return std::nullopt;
    }

    inline std::optional<int> find_residue(const std::string& name, int number) const {
        auto it = residue_map_.find(std::make_pair(name, number));
        if (it != residue_map_.end()) {
            return it->second;
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

    // New methods for additional PSF sections
    void add_title(const std::string& title) { titles_.push_back(title); }
    void add_improper(int atom1, int atom2, int atom3, int atom4,
                     double angle = 0.0, double force_constant = 0.0) {
        if (!has_atom(atom1) || !has_atom(atom2) || !has_atom(atom3) || !has_atom(atom4)) {
            throw std::out_of_range("Invalid atom indices in add_improper");
        }

        TopologyDihedral improper;
        improper.atom1 = atom1;
        improper.atom2 = atom2;
        improper.atom3 = atom3;
        improper.atom4 = atom4;
        improper.angle = angle;
        improper.force_constant = force_constant;
        improper.multiplicity = 0;  // Not used for impropers
        improper.improper = true;
        dihedrals_.push_back(improper);
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
        exclusions_[atom2].insert(atom1);  // Exclusions are symmetric
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
        cmaps_.push_back(cmap);
    }

    // Getters for new sections
    const std::vector<std::string>& get_titles() const { return titles_; }
    const std::vector<TopologyDonor>& get_donors() const { return donors_; }
    const std::vector<TopologyAcceptor>& get_acceptors() const { return acceptors_; }
    const std::vector<TopologyGroup>& get_groups() const { return groups_; }
    const std::vector<TopologyCmap>& get_cmaps() const { return cmaps_; }
    const std::map<int, std::set<int>>& get_exclusions() const { return exclusions_; }

    // Bond methods
    inline size_t get_num_bonds() const { return bonds_.size(); }
    inline bool has_bond(int atom1, int atom2) const {
        for (const auto& bond : bonds_) {
            if ((bond.atom1 == atom1 && bond.atom2 == atom2) ||
                (bond.atom1 == atom2 && bond.atom2 == atom1)) {
                return true;
            }
        }
        return false;
    }

    // Angle methods
    inline size_t get_num_angles() const { return angles_.size(); }
    inline bool has_angle(int atom1, int atom2, int atom3) const {
        for (const auto& angle : angles_) {
            if ((angle.atom1 == atom1 && angle.atom2 == atom2 && angle.atom3 == atom3) ||
                (angle.atom1 == atom3 && angle.atom2 == atom2 && angle.atom3 == atom1)) {
                return true;
            }
        }
        return false;
    }

    // Dihedral methods
    inline size_t get_num_dihedrals() const { 
        size_t count = 0;
        for (const auto& dihedral : dihedrals_) {
            if (!dihedral.improper) count++;
        }
        return count;
    }

    inline bool has_dihedral(int atom1, int atom2, int atom3, int atom4) const {
        for (const auto& dihedral : dihedrals_) {
            if (dihedral.improper) continue;  // Skip impropers
            if ((dihedral.atom1 == atom1 && dihedral.atom2 == atom2 &&
                 dihedral.atom3 == atom3 && dihedral.atom4 == atom4) ||
                (dihedral.atom1 == atom4 && dihedral.atom2 == atom3 &&
                 dihedral.atom3 == atom2 && dihedral.atom4 == atom1)) {
                return true;
            }
        }
        return false;
    }

    // Improper methods
    inline size_t get_num_impropers() const { 
        size_t count = 0;
        for (const auto& dihedral : dihedrals_) {
            if (dihedral.improper) count++;
        }
        return count;
    }

    inline bool has_improper(int atom1, int atom2, int atom3, int atom4) const {
        // atom1 is assumed to be the central atom in the query
        std::vector<int> query_others = {atom2, atom3, atom4};
        std::sort(query_others.begin(), query_others.end());

        for (const auto& dihedral : dihedrals_) {
            if (!dihedral.improper) continue;  // Skip regular dihedrals
            
            // Try each position as the potential central atom
            std::array<int, 4> atoms = {dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4};
            for (int i = 0; i < 4; ++i) {
                if (atoms[i] == atom1) {  // Found potential central atom match
                    std::vector<int> others;
                    for (int j = 0; j < 4; ++j) {
                        if (j != i) {
                            others.push_back(atoms[j]);
                        }
                    }
                    std::sort(others.begin(), others.end());
                    if (others == query_others) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    // Donor/Acceptor methods
    inline size_t get_num_donors() const { return donors_.size(); }
    inline size_t get_num_acceptors() const { return acceptors_.size(); }
    inline bool has_donor(int donor_atom, int hydrogen_atom) const {
        for (const auto& donor : donors_) {
            if (donor.donor_atom == donor_atom && donor.hydrogen_atom == hydrogen_atom) {
                return true;
            }
        }
        return false;
    }
    inline bool has_acceptor(int acceptor_atom) const {
        for (const auto& acceptor : acceptors_) {
            if (acceptor.acceptor_atom == acceptor_atom) {
                return true;
            }
        }
        return false;
    }

    // CMAP methods
    inline size_t get_num_cmaps() const { return cmaps_.size(); }
    inline bool has_cmap(const std::vector<int>& atoms) const {
        // CHARMM format requires 8 atoms for CMAP terms:
        // C(i-1), N(i), CA(i), C(i), N(i+1), CA(i+1), C(i+1), N(i+2)
        if (atoms.size() != 8) return false;
        for (const auto& cmap : cmaps_) {
            bool match = true;
            for (size_t i = 0; i < 8; ++i) {
                if (cmap.atoms[i] != atoms[i]) {
                    match = false;
                    break;
                }
            }
            if (match) return true;
        }
        return false;
    }

    // Group methods
    inline size_t get_num_groups() const { return groups_.size(); }
    inline bool has_group(const std::vector<int>& atoms) const {
        for (const auto& group : groups_) {
            if (group.atoms == atoms) {
                return true;
            }
        }
        return false;
    }

    inline const TopologyGroup& get_group(int index) const {
        if (index < 0 || index >= static_cast<int>(groups_.size())) {
            throw std::out_of_range("Invalid group index");
        }
        return groups_[index];
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
    std::vector<TopologyDihedral> dihedrals_;  // Contains both dihedrals and impropers
    std::vector<TopologyDonor> donors_;
    std::vector<TopologyAcceptor> acceptors_;
    std::map<int, std::set<int>> exclusions_;  // atom_index -> set of excluded atom indices
    std::vector<TopologyGroup> groups_;
    std::vector<TopologyCmap> cmaps_;

    // Lookup maps for efficient searching
    std::unordered_map<std::string, int> segment_map_;  // segment_name -> index
    std::map<std::pair<std::string, int>, int> residue_map_;  // (residue_name, number) -> index
    std::map<std::tuple<std::string, int, std::string>, int> atom_map_;  // (residue_name, number, atom_name) -> index

    // Additional PSF sections
    std::vector<std::string> titles_;
};

} // namespace model
} // namespace pygcmc
