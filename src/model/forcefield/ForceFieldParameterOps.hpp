#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_PARAMETER_OPS_HPP
#define PYGCMC_MODEL_FORCEFIELD_PARAMETER_OPS_HPP

#include "ForceFieldCore.hpp"

namespace pygcmc {
namespace model {
namespace forcefield {

// === Parameter Addition Method Implementations ===

inline void ForceField::add_atom_mass(const std::string& type, double mass) {
    manager_.add_atom_mass(type, mass);
}

inline void ForceField::add_lj_params(const std::string& type, double epsilon, double rmin_half) {
    manager_.add_lj_params(type, epsilon, rmin_half);
}

inline void ForceField::add_nbfix(const std::string& type1, const std::string& type2, 
                                 double epsilon, double rmin) {
    manager_.add_nbfix(type1, type2, epsilon, rmin);
}

inline void ForceField::add_bond_params(const std::string& type1, const std::string& type2, 
                                        double kb, double b0) {
    manager_.add_bond_params(type1, type2, kb, b0);
}

inline void ForceField::add_angle_params(const std::string& type1, const std::string& type2, 
                                         const std::string& type3, double ktheta, double theta0, 
                                         double kub, double s0) {
    manager_.add_angle_params(type1, type2, type3, ktheta, theta0, kub, s0);
}

inline void ForceField::add_dihedral_params(const std::string& type1, const std::string& type2,
                                            const std::string& type3, const std::string& type4,
                                            double kchi, int n, double delta) {
    manager_.add_dihedral_params(type1, type2, type3, type4, kchi, n, delta);
}

inline void ForceField::add_improper_params(const std::string& type1, const std::string& type2,
                                            const std::string& type3, const std::string& type4,
                                            double kpsi, double psi0) {
    manager_.add_improper_params(type1, type2, type3, type4, kpsi, psi0);
}

// === Parameter Retrieval Method Implementations ===

inline double ForceField::get_atom_mass(const std::string& type) const {
    return manager_.get_atom_mass(type);
}

inline const LJParams& ForceField::get_lj_params(const std::string& type) const {
    return manager_.get_lj_params(type);
}

inline std::pair<NBFIXParams, bool> ForceField::get_nbfix(const std::string& type1, 
                                                         const std::string& type2) const {
    return manager_.get_nbfix(type1, type2);
}

inline const BondParams& ForceField::get_bond_params(const std::string& type1, 
                                                     const std::string& type2) const {
    return manager_.get_bond_params(type1, type2);
}

inline const AngleParams& ForceField::get_angle_params(const std::string& type1,
                                                       const std::string& type2,
                                                       const std::string& type3) const {
    return manager_.get_angle_params(type1, type2, type3);
}

inline const std::vector<DihedralParams>& ForceField::get_dihedral_params(const std::string& type1,
                                                                          const std::string& type2,
                                                                          const std::string& type3,
                                                                          const std::string& type4) const {
    return manager_.get_dihedral_params(type1, type2, type3, type4);
}

inline const ImproperParams& ForceField::get_improper_params(const std::string& type1, 
                                                             const std::string& type2,
                                                             const std::string& type3, 
                                                             const std::string& type4) const {
    return manager_.get_improper_params(type1, type2, type3, type4);
}

// === Direct Map Access Method Implementations ===

inline const std::map<std::string, double>& ForceField::get_atom_masses() const {
    return analyzer_.get_atom_masses();
}

inline const std::map<std::string, LJParams>& ForceField::get_lj_params_map() const {
    return analyzer_.get_lj_params();
}

inline const std::map<std::pair<std::string, std::string>, NBFIXParams>& ForceField::get_nbfix_map() const {
    return analyzer_.get_nbfix();
}

inline const std::map<std::pair<std::string, std::string>, BondParams>& ForceField::get_bond_params_map() const {
    return analyzer_.get_bond_params();
}

inline const std::map<std::tuple<std::string, std::string, std::string>, AngleParams>& ForceField::get_angle_params_map() const {
    return analyzer_.get_angle_params();
}

inline const std::map<std::tuple<std::string, std::string, std::string, std::string>, std::vector<DihedralParams>>& ForceField::get_dihedral_params_map() const {
    return analyzer_.get_dihedral_params();
}

inline const std::map<std::tuple<std::string, std::string, std::string, std::string>, ImproperParams>& ForceField::get_improper_params_map() const {
    return analyzer_.get_improper_params();
}

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_PARAMETER_OPS_HPP 