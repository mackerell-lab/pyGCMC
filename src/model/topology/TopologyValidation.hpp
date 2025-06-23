#pragma once

#ifndef PYGCMC_MODEL_TOPOLOGY_VALIDATION_HPP
#define PYGCMC_MODEL_TOPOLOGY_VALIDATION_HPP

#include "TopologyCore.hpp"
#include "TopologyValidationMixin.hpp"

namespace pygcmc {
namespace model {
namespace topology {

/**
 * @brief Topology validation utilities
 */
class TopologyValidator : public ValidationMixin {
public:
    TopologyValidator(const TopologyStorage& storage) : ValidationMixin(storage) {}

    /**
     * @brief Validate topology consistency
     */
    bool validate_topology() const {
        return validate_atoms() && validate_bonds() && 
               validate_angles() && validate_dihedrals();
    }

    /**
     * @brief Validate atom references
     */
    bool validate_atoms() const {
        for (const auto& atom : get_storage().atoms) {
            if (!is_valid_residue(atom.residue_id) || !is_valid_segment(atom.segment_id)) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Validate bond references
     */
    bool validate_bonds() const {
        for (const auto& bond : get_storage().bonds) {
            if (!are_valid_bond_atoms(bond.atom1, bond.atom2)) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Validate angle references
     */
    bool validate_angles() const {
        for (const auto& angle : get_storage().angles) {
            if (!are_valid_angle_atoms(angle.atom1, angle.atom2, angle.atom3)) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Validate dihedral references
     */
    bool validate_dihedrals() const {
        for (const auto& dihedral : get_storage().dihedrals) {
            if (!are_valid_dihedral_atoms(dihedral.atom1, dihedral.atom2, 
                                         dihedral.atom3, dihedral.atom4)) {
                return false;
            }
        }
        return true;
    }

    // Public validation interface - these are the ONLY validation methods that should be used externally
    bool has_atom(int index) const { return is_valid_atom(index); }
    bool has_residue(int index) const { return is_valid_residue(index); }
    bool has_segment(int index) const { return is_valid_segment(index); }
};

} // namespace topology
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_TOPOLOGY_VALIDATION_HPP 