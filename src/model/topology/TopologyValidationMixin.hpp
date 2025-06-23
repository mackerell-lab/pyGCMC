#pragma once

#ifndef PYGCMC_MODEL_TOPOLOGY_VALIDATION_MIXIN_HPP
#define PYGCMC_MODEL_TOPOLOGY_VALIDATION_MIXIN_HPP

#include "TopologyCore.hpp"

namespace pygcmc {
namespace model {
namespace topology {

/**
 * @brief Validation mixin for topology classes
 * 
 * This class provides common validation methods that can be used by
 * multiple topology classes to avoid code duplication.
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

private:
    const TopologyStorage& storage_;
};

} // namespace topology
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_TOPOLOGY_VALIDATION_MIXIN_HPP 