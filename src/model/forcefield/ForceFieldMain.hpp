#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_MAIN_HPP
#define PYGCMC_MODEL_FORCEFIELD_MAIN_HPP

#include "ForceFieldTypes.hpp"
#include <stdexcept>
#include <algorithm>
#include <optional>
#include <memory>

namespace pygcmc {
namespace model {
namespace forcefield {

class ForceFieldAnalysis; // Forward declaration

/**
 * @brief Main force field class that holds all force field parameters
 */
class ForceField {
public:
    ForceField() = default;
    ~ForceField() = default;

    // === Parameter Addition Methods ===
    void add_atom_mass(const std::string& type, double mass);
    void add_lj_params(const std::string& type, double epsilon, double rmin_half);
    void add_nbfix(const std::string& type1, const std::string& type2, double epsilon, double rmin);
    void add_bond_params(const std::string& type1, const std::string& type2, double kb, double b0);
    void add_angle_params(const std::string& type1, const std::string& type2, const std::string& type3,
                         double ktheta, double theta0, double kub = 0.0, double s0 = 0.0);
    void add_dihedral_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4,
                            double kchi, int n, double delta);
    void add_improper_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4,
                            double kpsi, double psi0);
    
    // === Drude-specific Parameter Addition Methods ===
    void add_alpha_thole_params(const std::string& type, double alpha, double thole);
    void add_lonepair(const LonePairParams& params);
    void add_anisotropy(const AnisotropyParams& params);

    // === Parameter Retrieval Methods ===
    double get_atom_mass(const std::string& type) const;
    const LJParams& get_lj_params(const std::string& type) const;
    std::pair<NBFIXParams, bool> get_nbfix(const std::string& type1, const std::string& type2) const;
    const BondParams& get_bond_params(const std::string& type1, const std::string& type2) const;
    const AngleParams& get_angle_params(const std::string& type1, const std::string& type2,
                                      const std::string& type3) const;
    const std::vector<DihedralParams>& get_dihedral_params(const std::string& type1,
                                                          const std::string& type2,
                                                          const std::string& type3,
                                                          const std::string& type4) const;
    const ImproperParams& get_improper_params(const std::string& type1, const std::string& type2,
                                            const std::string& type3, const std::string& type4) const;
    
    // === Drude-specific Parameter Retrieval Methods ===
    const AlphaTHoleParams& get_alpha_params(const std::string& type) const;
    const std::vector<LonePairParams>& get_lonepairs() const { return lonepairs_; }
    const std::vector<AnisotropyParams>& get_anisotropies() const { return anisotropies_; }

    // === Existence Check Methods ===
    bool has_atom_mass(const std::string& type) const;
    bool has_lj_params(const std::string& type) const;
    bool has_nbfix(const std::string& type1, const std::string& type2) const;
    bool has_bond_params(const std::string& type1, const std::string& type2) const;
    bool has_angle_params(const std::string& type1, const std::string& type2,
                         const std::string& type3) const;
    bool has_dihedral_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4) const;
    bool has_improper_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4) const;
    
    // === Drude-specific Existence Check Methods ===
    bool has_alpha_params(const std::string& type) const;

    // === Size Methods ===
    size_t get_num_atom_types() const { return atom_masses_.size(); }
    size_t get_num_lj_params() const { return lj_params_.size(); }
    size_t get_num_nbfix() const { return nbfix_.size(); }
    size_t get_num_bond_types() const { return bond_params_.size(); }
    size_t get_num_angle_types() const { return angle_params_.size(); }
    size_t get_num_dihedral_types() const { return dihedral_params_.size(); }
    size_t get_num_improper_types() const { return improper_params_.size(); }
    
    // === Drude-specific Size Methods ===
    size_t get_num_alpha_params() const { return alpha_thole_params_.size(); }
    size_t get_num_lonepairs() const { return lonepairs_.size(); }
    size_t get_num_anisotropies() const { return anisotropies_.size(); }


    // === Nonbonded Parameters ===
    const NonbondedParams& get_nonbonded_params() const { return nonbonded_params_; }
    NonbondedParams& get_nonbonded_params() { return nonbonded_params_; }

    // === Direct access to parameter maps (for Python bindings) ===
    const std::map<std::string, double>& get_atom_masses() const { return atom_masses_; }
    const std::map<std::string, LJParams>& get_lj_params_map() const { return lj_params_; }
    const std::map<std::pair<std::string, std::string>, NBFIXParams>& get_nbfix_map() const { return nbfix_; }
    const std::map<std::pair<std::string, std::string>, BondParams>& get_bond_params_map() const { return bond_params_; }
    const std::map<std::tuple<std::string, std::string, std::string>, AngleParams>& get_angle_params_map() const { return angle_params_; }
    const std::map<std::tuple<std::string, std::string, std::string, std::string>, std::vector<DihedralParams>>& get_dihedral_params_map() const { return dihedral_params_; }
    const std::map<std::tuple<std::string, std::string, std::string, std::string>, ImproperParams>& get_improper_params_map() const { return improper_params_; }

    // === Static helper methods for making parameter keys ===
    static std::pair<std::string, std::string> makeTypePair(const std::string& type1, const std::string& type2) {
        return type1 < type2 ? std::make_pair(type1, type2) : std::make_pair(type2, type1);
    }

    static std::tuple<std::string, std::string, std::string> makeTypeTriple(
        const std::string& type1, const std::string& type2, const std::string& type3) {
        return std::make_tuple(type1, type2, type3);
    }

    static std::tuple<std::string, std::string, std::string, std::string> makeTypeQuad(
        const std::string& type1, const std::string& type2,
        const std::string& type3, const std::string& type4) {
        return std::make_tuple(type1, type2, type3, type4);
    }

private:
    std::map<std::string, double> atom_masses_;
    std::map<std::string, LJParams> lj_params_;
    std::map<std::pair<std::string, std::string>, NBFIXParams> nbfix_;
    std::map<std::pair<std::string, std::string>, BondParams> bond_params_;
    std::map<std::tuple<std::string, std::string, std::string>, AngleParams> angle_params_;
    std::map<std::tuple<std::string, std::string, std::string, std::string>, std::vector<DihedralParams>> dihedral_params_;
    std::map<std::tuple<std::string, std::string, std::string, std::string>, ImproperParams> improper_params_;
    NonbondedParams nonbonded_params_;
    
    // Drude-specific parameters
    std::map<std::string, AlphaTHoleParams> alpha_thole_params_;
    std::vector<LonePairParams> lonepairs_;
    std::vector<AnisotropyParams> anisotropies_;
};

} // namespace forcefield
} // namespace model
} // namespace pygcmc

// Include implementations
#include "ForceFieldAccessors.hpp"

#endif // PYGCMC_MODEL_FORCEFIELD_MAIN_HPP 