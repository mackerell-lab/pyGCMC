#pragma once

#include "TopologyMain.hpp"
#include <stdexcept>
#include <iostream>
#include <algorithm>

namespace pygcmc {
namespace model {
namespace topology {

// Implementation of inline functions for Topology class

inline int Topology::add_atom(const std::string& name, const std::string& type, double charge, double mass,
            const std::string& residue_name, int residue_number, const std::string& segment_name) {
    // Ensure we have the segment
    int segment_id;
    auto segment_it = segment_map_.find(segment_name);
    if (segment_it == segment_map_.end()) {
        TopologySegment segment;
        segment.id = static_cast<int>(segments_.size());
        segment.name = segment_name;
        segment_id = segment.id;
        segments_.push_back(segment);
        segment_map_[segment_name] = segment_id;
    } else {
        segment_id = segment_it->second;
    }

    // Ensure we have the residue
    int residue_id;
    auto residue_key = std::make_tuple(residue_name, residue_number, segment_name);
    auto residue_it = residue_map_.find(residue_key);
    if (residue_it == residue_map_.end()) {
        TopologyResidue residue;
        residue.id = static_cast<int>(residues_.size());
        residue.name = residue_name;
        residue.number = residue_number;
        residue.segment = segment_name;
        residue_id = residue.id;
        residues_.push_back(residue);
        residue_map_[residue_key] = residue_id;
        segments_[segment_id].residues.push_back(residue_id);
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

    int atom_id = static_cast<int>(atoms_.size());
    atoms_.push_back(atom);
    atom_map_[std::make_tuple(residue_name, residue_number, segment_name, name)] = atom_id;
    residues_[residue_id].atoms.push_back(atom_id);

    return atom_id;
}

inline void Topology::add_bond(int atom1, int atom2, double length, double force_constant, int function_type) {
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

inline void Topology::add_angle(int atom1, int atom2, int atom3, double angle, double force_constant, int function_type) {
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
    angles_.push_back(ang);
}

inline void Topology::add_dihedral(int atom1, int atom2, int atom3, int atom4, int multiplicity,
                 double angle, double force_constant, bool improper, int function_type) {
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

inline int Topology::add_residue(const std::string& name, int number, const std::string& segment) {
    // Ensure we have the segment
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

    int residue_id = residues_.size();
    residues_.push_back(residue);
    residue_map_[std::make_tuple(name, number, segment)] = residue_id;
    segments_[segment_id].residues.push_back(residue_id);

    return residue_id;
}

inline int Topology::add_segment(const std::string& name) {
    TopologySegment segment;
    segment.id = segments_.size();
    segment.name = name;

    int segment_id = segments_.size();
    segments_.push_back(segment);
    segment_map_[name] = segment_id;

    return segment_id;
}

inline void Topology::add_donor(int donor, int hydrogen) {
    if (!has_atom(donor) || !has_atom(hydrogen)) {
        throw std::out_of_range("Invalid atom indices in add_donor");
    }
    TopologyDonor d;
    d.donor_atom = donor;
    d.hydrogen_atom = hydrogen;
    donors_.push_back(d);
}

inline void Topology::add_acceptor(int acceptor) {
    if (!has_atom(acceptor)) {
        throw std::out_of_range("Invalid atom index in add_acceptor");
    }
    TopologyAcceptor a;
    a.acceptor_atom = acceptor;
    acceptors_.push_back(a);
}

inline void Topology::add_nonbonded_exclusion(int atom1, int atom2) {
    if (!has_atom(atom1) || !has_atom(atom2)) {
        throw std::out_of_range("Invalid atom indices in add_nonbonded_exclusion");
    }
    exclusions_[atom1].insert(atom2);
    exclusions_[atom2].insert(atom1);
}

inline void Topology::add_group(int id, const std::vector<int>& atoms, const std::string& type) {
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

inline void Topology::add_cmap(const std::array<int, 8>& atoms) {
    for (int atom : atoms) {
        if (atom >= 0 && !has_atom(atom)) {
            throw std::out_of_range("Invalid atom index in add_cmap");
        }
    }
    TopologyCmap cmap;
    cmap.atoms = atoms;
    cmap.function_type = 1;
    cmaps_.push_back(cmap);
}

inline void Topology::add_cmap(const std::array<int, 5>& atoms, int function_type) {
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

inline std::optional<int> Topology::find_atom(const std::string& residue_name, int residue_number, const std::string& atom_name) const {
    for (const auto& segment : segments_) {
        auto it = atom_map_.find(std::make_tuple(residue_name, residue_number, segment.name, atom_name));
        if (it != atom_map_.end()) {
            return it->second;
        }
    }
    return std::nullopt;
}

inline std::optional<int> Topology::find_residue(const std::string& name, int number) const {
    for (const auto& segment : segments_) {
        auto it = residue_map_.find(std::make_tuple(name, number, segment.name));
        if (it != residue_map_.end()) {
            return it->second;
        }
    }
    return std::nullopt;
}

inline std::optional<int> Topology::find_segment(const std::string& name) const {
    auto it = segment_map_.find(name);
    return (it != segment_map_.end()) ? std::optional<int>(it->second) : std::nullopt;
}

inline bool Topology::has_bond(int atom1, int atom2) const {
    for (const auto& bond : bonds_) {
        if ((bond.atom1 == atom1 && bond.atom2 == atom2) ||
            (bond.atom1 == atom2 && bond.atom2 == atom1)) {
            return true;
        }
    }
    return false;
}

inline bool Topology::has_angle(int atom1, int atom2, int atom3) const {
    for (const auto& angle : angles_) {
        if ((angle.atom1 == atom1 && angle.atom2 == atom2 && angle.atom3 == atom3) ||
            (angle.atom1 == atom3 && angle.atom2 == atom2 && angle.atom3 == atom1)) {
            return true;
        }
    }
    return false;
}

inline bool Topology::has_dihedral(int atom1, int atom2, int atom3, int atom4) const {
    for (const auto& dihedral : dihedrals_) {
        if (dihedral.improper) continue;
        if ((dihedral.atom1 == atom1 && dihedral.atom2 == atom2 &&
             dihedral.atom3 == atom3 && dihedral.atom4 == atom4) ||
            (dihedral.atom1 == atom4 && dihedral.atom2 == atom3 &&
             dihedral.atom3 == atom2 && dihedral.atom4 == atom1)) {
            return true;
        }
    }
    return false;
}

inline bool Topology::has_improper(int atom1, int atom2, int atom3, int atom4) const {
    std::vector<int> query_others = {atom2, atom3, atom4};
    std::sort(query_others.begin(), query_others.end());

    for (const auto& dihedral : dihedrals_) {
        if (!dihedral.improper) continue;
        
        std::array<int, 4> atoms = {dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4};
        for (int i = 0; i < 4; ++i) {
            if (atoms[i] == atom1) {
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

inline bool Topology::has_donor(int donor_atom) const {
    for (const auto& donor : donors_) {
        if (donor.donor_atom == donor_atom) {
            return true;
        }
    }
    return false;
}

inline bool Topology::has_donor(int donor_atom, int hydrogen_atom) const {
    for (const auto& donor : donors_) {
        if (donor.donor_atom == donor_atom && donor.hydrogen_atom == hydrogen_atom) {
            return true;
        }
    }
    return false;
}

inline bool Topology::has_acceptor(int acceptor_atom) const {
    for (const auto& acceptor : acceptors_) {
        if (acceptor.acceptor_atom == acceptor_atom) {
            return true;
        }
    }
    return false;
}

inline bool Topology::has_cmap(const std::vector<int>& atoms) const {
    if (atoms.size() == 8) {
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
    }
    else if (atoms.size() == 5) {
        for (const auto& cmap : cmaps_) {
            bool match = true;
            for (size_t i = 0; i < 5; ++i) {
                if (cmap.atoms[i] != atoms[i]) {
                    match = false;
                    break;
                }
            }
            if (match) return true;
        }
    }
    return false;
}

inline bool Topology::has_group(int group_id) const {
    for (const auto& group : groups_) {
        if (group.id == group_id) {
            return true;
        }
    }
    return false;
}

inline const TopologyGroup& Topology::get_group(int index) const {
    if (index < 0 || index >= static_cast<int>(groups_.size())) {
        throw std::out_of_range("Invalid group index");
    }
    return groups_[index];
}

} // namespace topology
} // namespace model
} // namespace pygcmc