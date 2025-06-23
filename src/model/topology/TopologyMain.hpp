#pragma once

#ifndef PYGCMC_MODEL_TOPOLOGY_MAIN_HPP
#define PYGCMC_MODEL_TOPOLOGY_MAIN_HPP

#include "TopologyCore.hpp"
#include "../common/ModelInterface.hpp"
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
 * @brief Main topology class - compact version under 300 lines
 */
class Topology : public common::IValidatable {
public:
    Topology() = default;
    ~Topology() = default;

    // IValidatable interface
    bool is_valid() const override {
        for (const auto& atom : storage_.atoms) {
            if (atom.residue_id < 0 || atom.residue_id >= static_cast<int>(storage_.residues.size()) ||
                atom.segment_id < 0 || atom.segment_id >= static_cast<int>(storage_.segments.size())) {
                return false;
            }
        }
        return true;
    }

    // Core add methods
    inline int add_atom(const std::string& name, const std::string& type, double charge, double mass,
                const std::string& residue_name, int residue_number, const std::string& segment_name) {
        try {
            int segment_id = ensure_segment(segment_name);
            int residue_id = ensure_residue(residue_name, residue_number, segment_name, segment_id);

            TopologyAtom atom;
            atom.id = static_cast<int>(storage_.atoms.size());
            atom.name = name;
            atom.type = type;
            atom.charge = charge;
            atom.mass = mass;
            atom.residue_id = residue_id;
            atom.segment_id = segment_id;

            int atom_id = static_cast<int>(storage_.atoms.size());
            storage_.atoms.push_back(atom);
            storage_.atom_map[std::make_tuple(residue_name, residue_number, segment_name, name)] = atom_id;

            if (residue_id >= 0 && residue_id < static_cast<int>(storage_.residues.size())) {
                storage_.residues[residue_id].atoms.push_back(atom_id);
            }

            return atom_id;
        } catch (const std::exception& e) {
            std::cerr << "Error in add_atom: " << e.what() << std::endl;
            throw;
        }
    }

    inline void add_bond(int atom1, int atom2, double length = 0.0, double force_constant = 0.0, int function_type = 1) {
        TopologyBond bond;
        bond.atom1 = atom1;
        bond.atom2 = atom2;
        bond.length = length;
        bond.force_constant = force_constant;
        bond.function_type = function_type;
        storage_.bonds.push_back(bond);
    }

    inline void add_angle(int atom1, int atom2, int atom3, double angle = 0.0, double force_constant = 0.0, int function_type = 1) {
        TopologyAngle ang;
        ang.atom1 = atom1; ang.atom2 = atom2; ang.atom3 = atom3;
        ang.angle = angle; ang.force_constant = force_constant; ang.function_type = function_type;
        ang.ub_length = 0.0; ang.ub_constant = 0.0;
        storage_.angles.push_back(ang);
    }

    inline void add_dihedral(int atom1, int atom2, int atom3, int atom4, int multiplicity = 1,
                     double angle = 0.0, double force_constant = 0.0, bool improper = false, int function_type = 1) {
        TopologyDihedral dihedral;
        dihedral.atom1 = atom1; dihedral.atom2 = atom2; dihedral.atom3 = atom3; dihedral.atom4 = atom4;
        dihedral.multiplicity = multiplicity; dihedral.angle = angle; dihedral.force_constant = force_constant;
        dihedral.improper = improper; dihedral.function_type = function_type;
        storage_.dihedrals.push_back(dihedral);
    }

    inline int add_residue(const std::string& name, int number, const std::string& segment) {
        int segment_id = ensure_segment(segment);
        TopologyResidue residue;
        residue.id = static_cast<int>(storage_.residues.size());
        residue.name = name; residue.number = number; residue.segment = segment;
        int residue_id = static_cast<int>(storage_.residues.size());
        storage_.residues.push_back(residue);
        storage_.residue_map[std::make_tuple(name, number, segment)] = residue_id;
        if (segment_id >= 0 && segment_id < static_cast<int>(storage_.segments.size())) {
            storage_.segments[segment_id].residues.push_back(residue_id);
        }
        return residue_id;
    }

    inline int add_segment(const std::string& name) {
        TopologySegment segment;
        segment.id = static_cast<int>(storage_.segments.size());
        segment.name = name;
        int segment_id = static_cast<int>(storage_.segments.size());
        storage_.segments.push_back(segment);
        storage_.segment_map[name] = segment_id;
        return segment_id;
    }

    // Getters
    inline const TopologyAtom& get_atom(int index) const {
        if (index < 0 || index >= static_cast<int>(storage_.atoms.size())) {
            throw std::out_of_range("Invalid atom index");
        }
        return storage_.atoms[index];
    }

    inline const TopologyResidue& get_residue(int index) const {
        if (index < 0 || index >= static_cast<int>(storage_.residues.size())) {
            throw std::out_of_range("Invalid residue index");
        }
        return storage_.residues[index];
    }

    inline const TopologySegment& get_segment(int index) const {
        if (index < 0 || index >= static_cast<int>(storage_.segments.size())) {
            throw std::out_of_range("Invalid segment index");
        }
        return storage_.segments[index];
    }

    inline const std::vector<TopologyBond>& get_bonds() const { return storage_.bonds; }
    inline const std::vector<TopologyAngle>& get_angles() const { return storage_.angles; }
    inline const std::vector<TopologyDihedral>& get_dihedrals() const { return storage_.dihedrals; }

    // Utility functions
    inline bool has_atom(int index) const { return index >= 0 && index < static_cast<int>(storage_.atoms.size()); }
    inline bool has_residue(int index) const { return index >= 0 && index < static_cast<int>(storage_.residues.size()); }
    inline bool has_segment(int index) const { return index >= 0 && index < static_cast<int>(storage_.segments.size()); }

    inline int get_num_atoms() const { return static_cast<int>(storage_.atoms.size()); }
    inline int get_num_residues() const { return static_cast<int>(storage_.residues.size()); }
    inline int get_num_segments() const { return static_cast<int>(storage_.segments.size()); }

    // Find elements
    std::optional<int> find_atom(const std::string& residue_name, int residue_number, const std::string& atom_name) const;
    std::optional<int> find_residue(const std::string& name, int number) const;
    std::optional<int> find_segment(const std::string& name) const;

    // Additional methods for compatibility
    void add_title(const std::string& title) { storage_.titles.push_back(title); }
    void add_improper(int atom1, int atom2, int atom3, int atom4, double angle = 0.0, double force_constant = 0.0) {
        add_dihedral(atom1, atom2, atom3, atom4, 0, angle, force_constant, true);
    }
    void add_donor(int donor, int hydrogen) {
        TopologyDonor d; d.donor_atom = donor; d.hydrogen_atom = hydrogen;
        storage_.donors.push_back(d);
    }
    void add_acceptor(int acceptor) {
        TopologyAcceptor a; a.acceptor_atom = acceptor;
        storage_.acceptors.push_back(a);
    }
    void add_nonbonded_exclusion(int atom1, int atom2) {
        storage_.exclusions[atom1].insert(atom2);
        storage_.exclusions[atom2].insert(atom1);
    }
    void add_group(int id, const std::vector<int>& atoms, const std::string& type = "") {
        TopologyGroup group; group.id = id; group.atoms = atoms; group.type = type;
        storage_.groups.push_back(group);
    }
    void add_cmap(const std::array<int, 8>& atoms) {
        TopologyCmap cmap; cmap.atoms = atoms; cmap.function_type = 1;
        storage_.cmaps.push_back(cmap);
    }
    void add_cmap(const std::array<int, 5>& atoms, int function_type = 1) {
        std::array<int, 8> charmm_atoms;
        for (int i = 0; i < 5; ++i) charmm_atoms[i] = atoms[i];
        for (int i = 5; i < 8; ++i) charmm_atoms[i] = -1;
        TopologyCmap cmap; cmap.atoms = charmm_atoms; cmap.function_type = function_type;
        storage_.cmaps.push_back(cmap);
    }

    // Additional getters
    const std::vector<std::string>& get_titles() const { return storage_.titles; }
    const std::vector<TopologyDonor>& get_donors() const { return storage_.donors; }
    const std::vector<TopologyAcceptor>& get_acceptors() const { return storage_.acceptors; }
    const std::vector<TopologyGroup>& get_groups() const { return storage_.groups; }
    const std::vector<TopologyCmap>& get_cmaps() const { return storage_.cmaps; }
    const std::map<int, std::set<int>>& get_exclusions() const { return storage_.exclusions; }

    // Count methods
    size_t get_num_bonds() const { return storage_.bonds.size(); }
    size_t get_num_angles() const { return storage_.angles.size(); }
    size_t get_num_dihedrals() const;
    size_t get_num_impropers() const;
    size_t get_num_donors() const { return storage_.donors.size(); }
    size_t get_num_acceptors() const { return storage_.acceptors.size(); }
    size_t get_num_cmaps() const { return storage_.cmaps.size(); }
    size_t get_num_groups() const { return storage_.groups.size(); }

    // Essential check methods
    bool has_cmap() const { return !storage_.cmaps.empty(); }
    bool has_cmap(const std::vector<int>& atoms) const;
    bool has_bond(int atom1, int atom2) const;
    bool has_angle(int atom1, int atom2, int atom3) const;
    bool has_dihedral(int atom1, int atom2, int atom3, int atom4) const;
    bool has_improper(int atom1, int atom2, int atom3, int atom4) const;
    bool has_donor(int donor_atom) const;
    bool has_donor(int donor_atom, int hydrogen_atom) const;
    bool has_acceptor(int acceptor_atom) const;
    bool has_group(int group_id) const;
    
    const TopologyGroup& get_group(int index) const {
        if (index < 0 || index >= static_cast<int>(storage_.groups.size())) {
            throw std::out_of_range("Invalid group index");
        }
        return storage_.groups[index];
    }

    void reserve_atoms(size_t n) { storage_.atoms.reserve(n); }

private:
    TopologyStorage storage_;

    int ensure_segment(const std::string& segment_name) {
        auto it = storage_.segment_map.find(segment_name);
        return (it == storage_.segment_map.end()) ? add_segment(segment_name) : it->second;
    }

    int ensure_residue(const std::string& residue_name, int residue_number, 
                      const std::string& segment_name, int segment_id) {
        auto residue_key = std::make_tuple(residue_name, residue_number, segment_name);
        auto it = storage_.residue_map.find(residue_key);
        if (it == storage_.residue_map.end()) {
            TopologyResidue residue;
            residue.id = static_cast<int>(storage_.residues.size());
            residue.name = residue_name; residue.number = residue_number; residue.segment = segment_name;
            int residue_id = residue.id;
            storage_.residues.push_back(residue);
            storage_.residue_map[residue_key] = residue_id;
            if (segment_id >= 0 && segment_id < static_cast<int>(storage_.segments.size())) {
                storage_.segments[segment_id].residues.push_back(residue_id);
            }
            return residue_id;
        }
        return it->second;
    }
};

} // namespace topology

// Backward compatibility
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

// Include implementation
#include "TopologyMainImpl.hpp"

#endif // PYGCMC_MODEL_TOPOLOGY_MAIN_HPP 