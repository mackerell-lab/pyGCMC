#pragma once

#ifndef PYGCMC_MODEL_TOPOLOGY_BONDS_HPP
#define PYGCMC_MODEL_TOPOLOGY_BONDS_HPP

#include "TopologyCore.hpp"
#include "TopologyValidationMixin.hpp"
#include <stdexcept>
#include <algorithm>

namespace pygcmc {
namespace model {
namespace topology {

/**
 * @brief Bond, angle, and dihedral management operations
 */
class TopologyBondManager : public ValidationMixin {
public:
    TopologyBondManager(TopologyStorage& storage) : ValidationMixin(storage), storage_(storage) {}

    /**
     * @brief Add a bond between two atoms
     */
    void add_bond(int atom1, int atom2, double length = 0.0, double force_constant = 0.0, int function_type = 1) {
        if (!are_valid_bond_atoms(atom1, atom2)) {
            throw std::invalid_argument("Invalid atom indices for bond");
        }
        
        TopologyBond bond;
        bond.atom1 = atom1;
        bond.atom2 = atom2;
        bond.length = length;
        bond.force_constant = force_constant;
        bond.function_type = function_type;
        storage_.bonds.push_back(bond);
    }

    /**
     * @brief Add an angle between three atoms
     */
    void add_angle(int atom1, int atom2, int atom3, double angle = 0.0, double force_constant = 0.0, int function_type = 1) {
        if (!are_valid_angle_atoms(atom1, atom2, atom3)) {
            throw std::invalid_argument("Invalid atom indices for angle");
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
        storage_.angles.push_back(ang);
    }

    /**
     * @brief Add angle with Urey-Bradley terms
     */
    void add_angle_with_ub(int atom1, int atom2, int atom3, double angle = 0.0, double force_constant = 0.0, 
                          double ub_length = 0.0, double ub_constant = 0.0, int function_type = 1) {
        if (!are_valid_angle_atoms(atom1, atom2, atom3)) {
            throw std::invalid_argument("Invalid atom indices for angle");
        }
        
        TopologyAngle ang;
        ang.atom1 = atom1;
        ang.atom2 = atom2;
        ang.atom3 = atom3;
        ang.angle = angle;
        ang.force_constant = force_constant;
        ang.function_type = function_type;
        ang.ub_length = ub_length;
        ang.ub_constant = ub_constant;
        ang.has_ub = (ub_constant != 0.0);
        storage_.angles.push_back(ang);
    }

    /**
     * @brief Add a dihedral between four atoms
     */
    void add_dihedral(int atom1, int atom2, int atom3, int atom4, int multiplicity = 1,
                     double angle = 0.0, double force_constant = 0.0, bool improper = false, int function_type = 1) {
        if (!are_valid_dihedral_atoms(atom1, atom2, atom3, atom4)) {
            throw std::invalid_argument("Invalid atom indices for dihedral");
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
        storage_.dihedrals.push_back(dihedral);
    }

    /**
     * @brief Add an improper dihedral
     */
    void add_improper(int atom1, int atom2, int atom3, int atom4,
                     double angle = 0.0, double force_constant = 0.0) {
        add_dihedral(atom1, atom2, atom3, atom4, 0, angle, force_constant, true);
    }

    /**
     * @brief Add a hydrogen bond donor
     */
    void add_donor(int donor, int hydrogen) {
        if (!are_valid_bond_atoms(donor, hydrogen)) {
            throw std::invalid_argument("Invalid atom indices for donor");
        }
        
        TopologyDonor d;
        d.donor_atom = donor;
        d.hydrogen_atom = hydrogen;
        storage_.donors.push_back(d);
    }

    /**
     * @brief Add a hydrogen bond acceptor
     */
    void add_acceptor(int acceptor) {
        if (!is_valid_atom(acceptor)) {
            throw std::invalid_argument("Invalid atom index for acceptor");
        }
        
        TopologyAcceptor a;
        a.acceptor_atom = acceptor;
        storage_.acceptors.push_back(a);
    }

    /**
     * @brief Add a nonbonded exclusion between two atoms
     */
    void add_nonbonded_exclusion(int atom1, int atom2) {
        if (!are_valid_bond_atoms(atom1, atom2)) {
            throw std::invalid_argument("Invalid atom indices for exclusion");
        }
        
        storage_.exclusions[atom1].insert(atom2);
        storage_.exclusions[atom2].insert(atom1);
    }

    /**
     * @brief Add a group of atoms
     */
    void add_group(int id, const std::vector<int>& atoms, const std::string& type = "") {
        // Validate all atoms in the group
        for (int atom : atoms) {
            if (!is_valid_atom(atom)) {
                throw std::invalid_argument("Invalid atom index in group");
            }
        }
        
        TopologyGroup group;
        group.id = id;
        group.atoms = atoms;
        group.type = type;
        storage_.groups.push_back(group);
    }

    /**
     * @brief Add a CMAP term (8-atom CHARMM format)
     */
    void add_cmap(const std::array<int, 8>& atoms) {
        // Validate atoms (only first 5 are required to be valid for CHARMM format)
        for (int i = 0; i < 5; ++i) {
            if (!is_valid_atom(atoms[i])) {
                throw std::invalid_argument("Invalid atom index in CMAP");
            }
        }
        
        TopologyCmap cmap;
        cmap.atoms = atoms;
        cmap.function_type = 1;
        storage_.cmaps.push_back(cmap);
    }

    /**
     * @brief Add a CMAP term (5-atom GROMACS format)
     */
    void add_cmap(const std::array<int, 5>& atoms, int function_type = 1) {
        // Validate all atoms
        for (int atom : atoms) {
            if (!is_valid_atom(atom)) {
                throw std::invalid_argument("Invalid atom index in CMAP");
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
        storage_.cmaps.push_back(cmap);
    }

    // Getters
    const std::vector<TopologyBond>& get_bonds() const { return storage_.bonds; }
    const std::vector<TopologyAngle>& get_angles() const { return storage_.angles; }
    const std::vector<TopologyDihedral>& get_dihedrals() const { return storage_.dihedrals; }
    const std::vector<TopologyDonor>& get_donors() const { return storage_.donors; }
    const std::vector<TopologyAcceptor>& get_acceptors() const { return storage_.acceptors; }
    const std::vector<TopologyGroup>& get_groups() const { return storage_.groups; }
    const std::vector<TopologyCmap>& get_cmaps() const { return storage_.cmaps; }
    const std::map<int, std::set<int>>& get_exclusions() const { return storage_.exclusions; }

    // Count methods
    size_t get_num_bonds() const { return storage_.bonds.size(); }
    size_t get_num_angles() const { return storage_.angles.size(); }
    size_t get_num_dihedrals() const { 
        size_t count = 0;
        for (const auto& dih : storage_.dihedrals) {
            if (!dih.improper) count++;
        }
        return count;
    }
    size_t get_num_impropers() const { 
        size_t count = 0;
        for (const auto& dih : storage_.dihedrals) {
            if (dih.improper) count++;
        }
        return count;
    }
    size_t get_num_donors() const { return storage_.donors.size(); }
    size_t get_num_acceptors() const { return storage_.acceptors.size(); }
    size_t get_num_cmaps() const { return storage_.cmaps.size(); }
    size_t get_num_groups() const { return storage_.groups.size(); }

    // Advanced check methods for topology elements
    bool has_donor(int donor_atom) const {
        return std::any_of(storage_.donors.begin(), storage_.donors.end(),
            [donor_atom](const TopologyDonor& donor) {
                return donor.donor_atom == donor_atom;
            });
    }

    bool has_donor(int donor_atom, int hydrogen_atom) const {
        return std::any_of(storage_.donors.begin(), storage_.donors.end(),
            [donor_atom, hydrogen_atom](const TopologyDonor& donor) {
                return donor.donor_atom == donor_atom && donor.hydrogen_atom == hydrogen_atom;
            });
    }

    bool has_acceptor(int acceptor_atom) const {
        return std::any_of(storage_.acceptors.begin(), storage_.acceptors.end(),
            [acceptor_atom](const TopologyAcceptor& acceptor) {
                return acceptor.acceptor_atom == acceptor_atom;
            });
    }

    // Check methods
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
                       ((dihedral.atom1 == atom1 && dihedral.atom2 == atom2 && 
                         dihedral.atom3 == atom3 && dihedral.atom4 == atom4) ||
                        (dihedral.atom1 == atom4 && dihedral.atom2 == atom3 && 
                         dihedral.atom3 == atom2 && dihedral.atom4 == atom1));
            });
    }

    bool has_cmap() const {
        return !storage_.cmaps.empty();
    }

    bool has_group(int group_id) const {
        return std::any_of(storage_.groups.begin(), storage_.groups.end(),
            [group_id](const TopologyGroup& group) {
                return group.id == group_id;
            });
    }

    const TopologyGroup& get_group(int index) const {
        if (index < 0 || index >= static_cast<int>(storage_.groups.size())) {
            throw std::out_of_range("Invalid group index");
        }
        return storage_.groups[index];
    }

private:
    TopologyStorage& storage_;
};

} // namespace topology
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_TOPOLOGY_BONDS_HPP 