#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_ACCESSORS_HPP
#define PYGCMC_MODEL_FORCEFIELD_ACCESSORS_HPP

#include "ForceFieldMain.hpp"

namespace pygcmc {
namespace model {
namespace forcefield {

// Inline implementations of all accessor methods for ForceField

// === Parameter Addition Methods ===

inline void ForceField::add_atom_mass(const std::string& type, double mass) {
    atom_masses_[type] = mass;
}

inline void ForceField::add_lj_params(const std::string& type, double epsilon, double rmin_half) {
    lj_params_[type] = LJParams{epsilon, rmin_half};
}

inline void ForceField::add_nbfix(const std::string& type1, const std::string& type2, 
                                 double epsilon, double rmin) {
    auto key = ForceField::makeTypePair(type1, type2);
    nbfix_[key] = NBFIXParams{epsilon, rmin};
}

inline void ForceField::add_bond_params(const std::string& type1, const std::string& type2, double kb, double b0) {
    auto key = ForceField::makeTypePair(type1, type2);
    bond_params_[key] = BondParams{kb, b0};
}

inline void ForceField::add_angle_params(const std::string& type1, const std::string& type2, const std::string& type3,
                                        double ktheta, double theta0, double kub, double s0) {
    auto key = ForceField::makeTypeTriple(type1, type2, type3);
    angle_params_[key] = AngleParams{ktheta, theta0, kub, s0};
}

inline void ForceField::add_dihedral_params(const std::string& type1, const std::string& type2,
                                           const std::string& type3, const std::string& type4,
                                           double kchi, int n, double delta) {
    auto key = ForceField::makeTypeQuad(type1, type2, type3, type4);
    dihedral_params_[key].emplace_back(DihedralParams{kchi, n, delta});
}

inline void ForceField::add_improper_params(const std::string& type1, const std::string& type2,
                                           const std::string& type3, const std::string& type4,
                                           double kpsi, double psi0) {
    auto key = ForceField::makeTypeQuad(type1, type2, type3, type4);
    improper_params_[key] = ImproperParams{kpsi, psi0};
}

// === Parameter Retrieval Methods ===

inline double ForceField::get_atom_mass(const std::string& type) const {
    auto it = atom_masses_.find(type);
    if (it == atom_masses_.end()) {
        throw std::runtime_error("Atom mass not found for type: " + type);
    }
    return it->second;
}

inline const LJParams& ForceField::get_lj_params(const std::string& type) const {
    auto it = lj_params_.find(type);
    if (it == lj_params_.end()) {
        throw std::runtime_error("LJ parameters not found for type: " + type);
    }
    return it->second;
}

inline std::pair<NBFIXParams, bool> ForceField::get_nbfix(const std::string& type1, 
                                                         const std::string& type2) const {
    auto key = ForceField::makeTypePair(type1, type2);
    auto it = nbfix_.find(key);
    if (it == nbfix_.end()) {
        return std::make_pair(NBFIXParams{}, false);
    }
    return std::make_pair(it->second, true);
}

inline const BondParams& ForceField::get_bond_params(const std::string& type1, const std::string& type2) const {
    auto key = ForceField::makeTypePair(type1, type2);
    auto it = bond_params_.find(key);
    if (it == bond_params_.end()) {
        throw std::runtime_error("Bond parameters not found for types: " + type1 + "-" + type2);
    }
    return it->second;
}

inline const AngleParams& ForceField::get_angle_params(const std::string& type1,
                                                      const std::string& type2,
                                                      const std::string& type3) const {
    auto key1 = std::make_tuple(type1, type2, type3);
    auto it = angle_params_.find(key1);
    if (it != angle_params_.end()) {
        return it->second;
    }
    
    auto key2 = std::make_tuple(type3, type2, type1);
    it = angle_params_.find(key2);
    if (it != angle_params_.end()) {
        return it->second;
    }
    
    throw std::runtime_error("Angle parameters not found for types: " + 
                           type1 + "-" + type2 + "-" + type3);
}

inline const std::vector<DihedralParams>& ForceField::get_dihedral_params(const std::string& type1,
                                                                         const std::string& type2,
                                                                         const std::string& type3,
                                                                         const std::string& type4) const {
    auto key = ForceField::makeTypeQuad(type1, type2, type3, type4);
    auto it = dihedral_params_.find(key);
    if (it == dihedral_params_.end()) {
        throw std::runtime_error("Dihedral parameters not found for types: " +
                               type1 + "-" + type2 + "-" + type3 + "-" + type4);
    }
    return it->second;
}

inline const ImproperParams& ForceField::get_improper_params(const std::string& type1, const std::string& type2,
                                                           const std::string& type3, const std::string& type4) const {
    auto key = ForceField::makeTypeQuad(type1, type2, type3, type4);
    auto it = improper_params_.find(key);
    if (it == improper_params_.end()) {
        throw std::runtime_error("Improper parameters not found for types: " +
                               type1 + "-" + type2 + "-" + type3 + "-" + type4);
    }
    return it->second;
}

// === Existence Check Methods ===

inline bool ForceField::has_atom_mass(const std::string& type) const {
    return atom_masses_.find(type) != atom_masses_.end();
}

inline bool ForceField::has_lj_params(const std::string& type) const {
    return lj_params_.find(type) != lj_params_.end();
}

inline bool ForceField::has_nbfix(const std::string& type1, const std::string& type2) const {
    auto key = ForceField::makeTypePair(type1, type2);
    return nbfix_.find(key) != nbfix_.end();
}

inline bool ForceField::has_bond_params(const std::string& type1, const std::string& type2) const {
    auto key = ForceField::makeTypePair(type1, type2);
    return bond_params_.find(key) != bond_params_.end();
}

inline bool ForceField::has_angle_params(const std::string& type1, const std::string& type2,
                                        const std::string& type3) const {
    auto key1 = std::make_tuple(type1, type2, type3);
    auto key2 = std::make_tuple(type3, type2, type1);
    return angle_params_.find(key1) != angle_params_.end() ||
           angle_params_.find(key2) != angle_params_.end();
}

inline bool ForceField::has_dihedral_params(const std::string& type1, const std::string& type2,
                                           const std::string& type3, const std::string& type4) const {
    auto key = ForceField::makeTypeQuad(type1, type2, type3, type4);
    return dihedral_params_.find(key) != dihedral_params_.end();
}

inline bool ForceField::has_improper_params(const std::string& type1, const std::string& type2,
                                           const std::string& type3, const std::string& type4) const {
    auto key = ForceField::makeTypeQuad(type1, type2, type3, type4);
    return improper_params_.find(key) != improper_params_.end();
}


} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_ACCESSORS_HPP 