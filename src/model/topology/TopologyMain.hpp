#pragma once

#ifndef PYGCMC_MODEL_TOPOLOGY_MAIN_HPP
#define PYGCMC_MODEL_TOPOLOGY_MAIN_HPP

#include "TopologyCore.hpp"
#include "TopologyAtoms.hpp"
#include "TopologyBonds.hpp"
#include "TopologyHBonds.hpp"
#include "TopologySpecial.hpp"
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
 * @brief Main topology class using composition of specialized managers
 */
class Topology : public common::IValidatable {
public:
    Topology() : 
        atom_manager_(storage_),
        bond_manager_(storage_),
        hbond_manager_(storage_),
        special_manager_(storage_) {}

    ~Topology() = default;

    // Disable copy and move semantics to prevent dangling manager references
    Topology(const Topology&) = delete;
    Topology& operator=(const Topology&) = delete;
    Topology(Topology&& other) noexcept : storage_(std::move(other.storage_)),
        atom_manager_(storage_), bond_manager_(storage_), hbond_manager_(storage_), special_manager_(storage_) {}
    Topology& operator=(Topology&&) = delete;

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

    // Atom management - delegate to AtomManager
    int add_atom(const std::string& name, const std::string& type, double charge, double mass,
                const std::string& residue_name, int residue_number, const std::string& segment_name) {
        return atom_manager_.add_atom(name, type, charge, mass, residue_name, residue_number, segment_name);
    }

    int add_residue(const std::string& name, int number, const std::string& segment) {
        return atom_manager_.add_residue(name, number, segment);
    }

    int add_segment(const std::string& name) {
        return atom_manager_.add_segment(name);
    }

    // Bond management - delegate to BondManager
    void add_bond(int atom1, int atom2, double length = 0.0, double force_constant = 0.0, int function_type = 1) {
        bond_manager_.add_bond(atom1, atom2, length, force_constant, function_type);
    }

    void add_angle(int atom1, int atom2, int atom3, double angle = 0.0, double force_constant = 0.0, int function_type = 1) {
        bond_manager_.add_angle(atom1, atom2, atom3, angle, force_constant, function_type);
    }

    void add_dihedral(int atom1, int atom2, int atom3, int atom4, int multiplicity = 1,
                     double angle = 0.0, double force_constant = 0.0, bool improper = false, int function_type = 1) {
        bond_manager_.add_dihedral(atom1, atom2, atom3, atom4, multiplicity, angle, force_constant, improper, function_type);
    }

    void add_improper(int atom1, int atom2, int atom3, int atom4, double angle = 0.0, double force_constant = 0.0) {
        bond_manager_.add_improper(atom1, atom2, atom3, atom4, angle, force_constant);
    }

    // Hydrogen bond management - delegate to HBondManager
    void add_donor(int donor, int hydrogen) {
        hbond_manager_.add_donor(donor, hydrogen);
    }

    void add_acceptor(int acceptor) {
        hbond_manager_.add_acceptor(acceptor);
    }

    // Special features - delegate to SpecialManager
    void add_nonbonded_exclusion(int atom1, int atom2) {
        special_manager_.add_nonbonded_exclusion(atom1, atom2);
    }

    void add_group(int id, const std::vector<int>& atoms, const std::string& type = "") {
        special_manager_.add_group(id, atoms, type);
    }

    void add_cmap(const std::array<int, 8>& atoms) {
        special_manager_.add_cmap(atoms);
    }

    void add_cmap(const std::array<int, 5>& atoms, int function_type = 1) {
        special_manager_.add_cmap(atoms, function_type);
    }

    // Getters - delegate to appropriate managers
    const TopologyAtom& get_atom(int index) const { return atom_manager_.get_atom(index); }
    const TopologyResidue& get_residue(int index) const { return atom_manager_.get_residue(index); }
    const TopologySegment& get_segment(int index) const { return atom_manager_.get_segment(index); }

    const std::vector<TopologyBond>& get_bonds() const { return bond_manager_.get_bonds(); }
    const std::vector<TopologyAngle>& get_angles() const { return bond_manager_.get_angles(); }
    const std::vector<TopologyDihedral>& get_dihedrals() const { return bond_manager_.get_dihedrals(); }

    const std::vector<TopologyDonor>& get_donors() const { return hbond_manager_.get_donors(); }
    const std::vector<TopologyAcceptor>& get_acceptors() const { return hbond_manager_.get_acceptors(); }

    const std::vector<TopologyGroup>& get_groups() const { return special_manager_.get_groups(); }
    const std::vector<TopologyCmap>& get_cmaps() const { return special_manager_.get_cmaps(); }
    const std::map<int, std::set<int>>& get_exclusions() const { return special_manager_.get_exclusions(); }

    // Utility functions
    bool has_atom(int index) const { return atom_manager_.has_atom(index); }
    bool has_residue(int index) const { return atom_manager_.has_residue(index); }
    bool has_segment(int index) const { return atom_manager_.has_segment(index); }

    int get_num_atoms() const { return atom_manager_.get_num_atoms(); }
    int get_num_residues() const { return atom_manager_.get_num_residues(); }
    int get_num_segments() const { return atom_manager_.get_num_segments(); }

    // Count methods - delegate to appropriate managers
    size_t get_num_bonds() const { return bond_manager_.get_num_bonds(); }
    size_t get_num_angles() const { return bond_manager_.get_num_angles(); }
    size_t get_num_dihedrals() const { return bond_manager_.get_num_dihedrals(); }
    size_t get_num_impropers() const { return bond_manager_.get_num_impropers(); }

    size_t get_num_donors() const { return hbond_manager_.get_num_donors(); }
    size_t get_num_acceptors() const { return hbond_manager_.get_num_acceptors(); }

    size_t get_num_cmaps() const { return special_manager_.get_num_cmaps(); }
    size_t get_num_groups() const { return special_manager_.get_num_groups(); }

    // Check methods - delegate to appropriate managers
    bool has_bond(int atom1, int atom2) const { return bond_manager_.has_bond(atom1, atom2); }
    bool has_angle(int atom1, int atom2, int atom3) const { return bond_manager_.has_angle(atom1, atom2, atom3); }
    bool has_dihedral(int atom1, int atom2, int atom3, int atom4) const { return bond_manager_.has_dihedral(atom1, atom2, atom3, atom4); }
    bool has_improper(int atom1, int atom2, int atom3, int atom4) const { return bond_manager_.has_improper(atom1, atom2, atom3, atom4); }

    bool has_donor(int donor_atom) const { return hbond_manager_.has_donor(donor_atom); }
    bool has_donor(int donor_atom, int hydrogen_atom) const { return hbond_manager_.has_donor(donor_atom, hydrogen_atom); }
    bool has_acceptor(int acceptor_atom) const { return hbond_manager_.has_acceptor(acceptor_atom); }
    
    bool has_cmap() const { return special_manager_.has_cmap(); }
    bool has_group(int group_id) const { return special_manager_.has_group(group_id); }

    // Find elements - delegate to AtomManager
    std::optional<int> find_atom(const std::string& residue_name, int residue_number, const std::string& atom_name) const {
        return atom_manager_.find_atom(residue_name, residue_number, atom_name);
    }
    std::optional<int> find_residue(const std::string& name, int number) const {
        return atom_manager_.find_residue(name, number);
    }
    std::optional<int> find_segment(const std::string& name) const {
        return atom_manager_.find_segment(name);
    }

    // Additional compatibility methods
    void add_title(const std::string& title) { storage_.titles.push_back(title); }
    const std::vector<std::string>& get_titles() const { return storage_.titles; }

    const TopologyGroup& get_group(int index) const { return special_manager_.get_group(index); }
    void reserve_atoms(size_t n) { atom_manager_.reserve_atoms(n); }

    // CMAP check with atoms vector - support both 5-atom and 8-atom formats
    bool has_cmap(const std::vector<int>& atoms) const {
        if (atoms.size() == 8) {
            // CHARMM format with 8 atoms
            return std::any_of(special_manager_.get_cmaps().begin(), special_manager_.get_cmaps().end(),
                [&atoms](const TopologyCmap& cmap) {
                    for (size_t i = 0; i < 8; ++i) {
                        if (cmap.atoms[i] != atoms[i]) return false;
                    }
                    return true;
                });
        }
        else if (atoms.size() == 5) {
            // GROMACS format with 5 atoms
            return std::any_of(special_manager_.get_cmaps().begin(), special_manager_.get_cmaps().end(),
                [&atoms](const TopologyCmap& cmap) {
                    for (size_t i = 0; i < 5; ++i) {
                        if (cmap.atoms[i] != atoms[i]) return false;
                    }
                    return true;
                });
        }
        return false;
    }

    // Access to specialized managers for advanced usage
    const TopologyAtomManager& get_atom_manager() const { return atom_manager_; }
    const TopologyBondManager& get_bond_manager() const { return bond_manager_; }
    const TopologyHBondManager& get_hbond_manager() const { return hbond_manager_; }
    const TopologySpecialManager& get_special_manager() const { return special_manager_; }

private:
    TopologyStorage storage_;
    TopologyAtomManager atom_manager_;
    TopologyBondManager bond_manager_;
    TopologyHBondManager hbond_manager_;
    TopologySpecialManager special_manager_;
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

#endif // PYGCMC_MODEL_TOPOLOGY_MAIN_HPP 