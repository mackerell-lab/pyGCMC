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
 * @brief Core bond, angle, and dihedral management operations
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

    // Getters for core bond structures
    const std::vector<TopologyBond>& get_bonds() const { return storage_.bonds; }
    const std::vector<TopologyAngle>& get_angles() const { return storage_.angles; }
    const std::vector<TopologyDihedral>& get_dihedrals() const { return storage_.dihedrals; }

    // Count methods for core structures
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

    // Check methods for core structures
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

private:
    TopologyStorage& storage_;
};

} // namespace topology
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_TOPOLOGY_BONDS_HPP 