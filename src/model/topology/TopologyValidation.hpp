#pragma once

#ifndef PYGCMC_MODEL_TOPOLOGY_VALIDATION_HPP
#define PYGCMC_MODEL_TOPOLOGY_VALIDATION_HPP

#include "TopologyCore.hpp"

namespace pygcmc {
namespace model {
namespace topology {

/**
 * @brief Comprehensive validation mixin for topology classes
 * 
 * This class provides both basic validation methods (for inheritance) 
 * and complete topology validation capabilities.
 */
class ValidationMixin {
public:
    ValidationMixin(const TopologyStorage& storage) : storage_(storage) {}
    
protected:
    /**
     * @brief Check if atom index is valid
     */
    bool is_valid_atom(int index) const {
        return index >= 0 && index < static_cast<int>(storage_.atoms.size());
    }

    /**
     * @brief Check if residue index is valid
     */
    bool is_valid_residue(int index) const {
        return index >= 0 && index < static_cast<int>(storage_.residues.size());
    }

    /**
     * @brief Check if segment index is valid
     */
    bool is_valid_segment(int index) const {
        return index >= 0 && index < static_cast<int>(storage_.segments.size());
    }

    /**
     * @brief Check if bond atoms are valid
     */
    bool are_valid_bond_atoms(int atom1, int atom2) const {
        return is_valid_atom(atom1) && is_valid_atom(atom2);
    }

    /**
     * @brief Check if angle atoms are valid
     */
    bool are_valid_angle_atoms(int atom1, int atom2, int atom3) const {
        return is_valid_atom(atom1) && is_valid_atom(atom2) && is_valid_atom(atom3);
    }

    /**
     * @brief Check if dihedral atoms are valid
     */
    bool are_valid_dihedral_atoms(int atom1, int atom2, int atom3, int atom4) const {
        return is_valid_atom(atom1) && is_valid_atom(atom2) && 
               is_valid_atom(atom3) && is_valid_atom(atom4);
    }

    /**
     * @brief Get storage reference for derived classes
     */
    const TopologyStorage& get_storage() const { return storage_; }

public:
    // Complete topology validation interface (can be used directly or through inheritance)
    
    /**
     * @brief Validate complete topology consistency
     */
    bool validate_topology() const {
        return validate_atoms() && validate_bonds() && 
               validate_angles() && validate_dihedrals();
    }

    /**
     * @brief Validate atom references
     */
    bool validate_atoms() const {
        for (const auto& atom : storage_.atoms) {
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
        for (const auto& bond : storage_.bonds) {
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
        for (const auto& angle : storage_.angles) {
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
        for (const auto& dihedral : storage_.dihedrals) {
            if (!are_valid_dihedral_atoms(dihedral.atom1, dihedral.atom2, 
                                         dihedral.atom3, dihedral.atom4)) {
                return false;
            }
        }
        return true;
    }

    // Public validation interface - standard has_* methods
    bool has_atom(int index) const { return is_valid_atom(index); }
    bool has_residue(int index) const { return is_valid_residue(index); }
    bool has_segment(int index) const { return is_valid_segment(index); }

private:
    const TopologyStorage& storage_;
};

/**
 * @brief Standalone topology validator (for backward compatibility)
 * 
 * This is a typedef to ValidationMixin for cases where a dedicated validator is needed.
 */
using TopologyValidator = ValidationMixin;

} // namespace topology
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_TOPOLOGY_VALIDATION_HPP 