#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_CHECKER_OPS_HPP
#define PYGCMC_MODEL_FORCEFIELD_CHECKER_OPS_HPP

#include "ForceFieldCore.hpp"

namespace pygcmc {
namespace model {
namespace forcefield {

// === Parameter Existence Check Method Implementations ===

inline bool ForceField::has_atom_mass(const std::string& type) const {
    return operations_.has_atom_mass(type);
}

inline bool ForceField::has_lj_params(const std::string& type) const {
    return operations_.has_lj_params(type);
}

inline bool ForceField::has_nbfix(const std::string& type1, const std::string& type2) const {
    return operations_.has_nbfix(type1, type2);
}

inline bool ForceField::has_bond_params(const std::string& type1, const std::string& type2) const {
    return operations_.has_bond_params(type1, type2);
}

inline bool ForceField::has_angle_params(const std::string& type1, const std::string& type2,
                                         const std::string& type3) const {
    return operations_.has_angle_params(type1, type2, type3);
}

inline bool ForceField::has_dihedral_params(const std::string& type1, const std::string& type2,
                                            const std::string& type3, const std::string& type4) const {
    return operations_.has_dihedral_params(type1, type2, type3, type4);
}

inline bool ForceField::has_improper_params(const std::string& type1, const std::string& type2,
                                            const std::string& type3, const std::string& type4) const {
    return operations_.has_improper_params(type1, type2, type3, type4);
}

// === Size Query Method Implementations ===

inline size_t ForceField::get_num_atom_types() const {
    return operations_.get_num_atom_types();
}

inline size_t ForceField::get_num_lj_params() const {
    return operations_.get_num_lj_params();
}

inline size_t ForceField::get_num_nbfix() const {
    return operations_.get_num_nbfix();
}

inline size_t ForceField::get_num_bond_types() const {
    return operations_.get_num_bond_types();
}

inline size_t ForceField::get_num_angle_types() const {
    return operations_.get_num_angle_types();
}

inline size_t ForceField::get_num_dihedral_types() const {
    return operations_.get_num_dihedral_types();
}

inline size_t ForceField::get_num_improper_types() const {
    return operations_.get_num_improper_types();
}

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_CHECKER_OPS_HPP 