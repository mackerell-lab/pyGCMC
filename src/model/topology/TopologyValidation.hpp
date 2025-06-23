#pragma once

#ifndef PYGCMC_MODEL_TOPOLOGY_VALIDATION_HPP
#define PYGCMC_MODEL_TOPOLOGY_VALIDATION_HPP

#include "TopologyCore.hpp"

namespace pygcmc {
namespace model {
namespace topology {

/**
 * @brief Topology validation utilities
 */
class TopologyValidator {
public:
    TopologyValidator(const TopologyStorage& storage) : storage_(storage) {}

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
        for (const auto& atom : storage_.atoms) {
            if (atom.residue_id < 0 || atom.residue_id >= static_cast<int>(storage_.residues.size()) ||
                atom.segment_id < 0 || atom.segment_id >= static_cast<int>(storage_.segments.size())) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Validate bond references
     */
    bool validate_bonds() const {
        const int num_atoms = static_cast<int>(storage_.atoms.size());
        for (const auto& bond : storage_.bonds) {
            if (bond.atom1 < 0 || bond.atom1 >= num_atoms ||
                bond.atom2 < 0 || bond.atom2 >= num_atoms) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Validate angle references
     */
    bool validate_angles() const {
        const int num_atoms = static_cast<int>(storage_.atoms.size());
        for (const auto& angle : storage_.angles) {
            if (angle.atom1 < 0 || angle.atom1 >= num_atoms ||
                angle.atom2 < 0 || angle.atom2 >= num_atoms ||
                angle.atom3 < 0 || angle.atom3 >= num_atoms) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Validate dihedral references
     */
    bool validate_dihedrals() const {
        const int num_atoms = static_cast<int>(storage_.atoms.size());
        for (const auto& dihedral : storage_.dihedrals) {
            if (dihedral.atom1 < 0 || dihedral.atom1 >= num_atoms ||
                dihedral.atom2 < 0 || dihedral.atom2 >= num_atoms ||
                dihedral.atom3 < 0 || dihedral.atom3 >= num_atoms ||
                dihedral.atom4 < 0 || dihedral.atom4 >= num_atoms) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Check if topology element exists
     */
    bool has_atom(int index) const {
        return index >= 0 && index < static_cast<int>(storage_.atoms.size());
    }

    bool has_residue(int index) const {
        return index >= 0 && index < static_cast<int>(storage_.residues.size());
    }

    bool has_segment(int index) const {
        return index >= 0 && index < static_cast<int>(storage_.segments.size());
    }

private:
    const TopologyStorage& storage_;
};

} // namespace topology
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_TOPOLOGY_VALIDATION_HPP 