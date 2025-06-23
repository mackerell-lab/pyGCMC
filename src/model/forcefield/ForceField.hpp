#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_HPP
#define PYGCMC_MODEL_FORCEFIELD_HPP

#include "ForceFieldTypes.hpp"
#include <stdexcept>
#include <algorithm>
#include <optional>
#include <memory>
#include <sstream>
#include <iomanip>

namespace pygcmc {
namespace model {
namespace forcefield {

/**
 * @brief Main force field class that holds all force field parameters
 * 
 * This class provides a simple, efficient interface for managing force field
 * parameters. It combines the simplicity of the original single-file design
 * with some useful analysis and validation features.
 */
class ForceField {
public:
    ForceField() = default;
    ~ForceField() = default;

    // === Parameter Addition Methods ===

    void add_atom_mass(const std::string& type, double mass) {
        atom_masses_[type] = mass;
    }

    /**
     * @brief Add Lennard-Jones parameters for an atom type
     * @param type The atom type
     * @param epsilon Well depth (kcal/mole)
     * @param rmin_half Rmin/2: HALF of the distance at minimum energy (Angstroms)
     */
    void add_lj_params(const std::string& type, double epsilon, double rmin_half) {
        lj_params_[type] = LJParams{epsilon, rmin_half};
    }

    /**
     * @brief Add NBFIX parameters for a specific pair of atom types
     * @param type1 First atom type
     * @param type2 Second atom type
     * @param epsilon Well depth (kcal/mole)
     * @param rmin Distance at minimum energy (Angstroms)
     */
    void add_nbfix(const std::string& type1, const std::string& type2, 
                  double epsilon, double rmin) {
        auto key = ParamKeyUtils::makeTypePair(type1, type2);
        nbfix_[key] = NBFIXParams{epsilon, rmin};
    }

    void add_bond_params(const std::string& type1, const std::string& type2, double kb, double b0) {
        auto key = ParamKeyUtils::makeTypePair(type1, type2);
        bond_params_[key] = BondParams{kb, b0};
    }

    void add_angle_params(const std::string& type1, const std::string& type2, const std::string& type3,
                         double ktheta, double theta0, double kub = 0.0, double s0 = 0.0) {
        auto key = ParamKeyUtils::makeTypeTriple(type1, type2, type3);
        angle_params_[key] = AngleParams{ktheta, theta0, kub, s0};
    }

    void add_dihedral_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4,
                            double kchi, int n, double delta) {
        auto key = ParamKeyUtils::makeTypeQuad(type1, type2, type3, type4);
        dihedral_params_[key].emplace_back(kchi, n, delta);
    }

    void add_improper_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4,
                            double kpsi, double psi0) {
        auto key = ParamKeyUtils::makeTypeQuad(type1, type2, type3, type4);
        improper_params_[key] = ImproperParams{kpsi, psi0};
    }

    // === Parameter Retrieval Methods ===

    double get_atom_mass(const std::string& type) const {
        auto it = atom_masses_.find(type);
        if (it == atom_masses_.end()) {
            throw std::runtime_error("Atom mass not found for type: " + type);
        }
        return it->second;
    }

    const LJParams& get_lj_params(const std::string& type) const {
        auto it = lj_params_.find(type);
        if (it == lj_params_.end()) {
            throw std::runtime_error("LJ parameters not found for type: " + type);
        }
        return it->second;
    }

    /**
     * @brief Get NBFIX parameters for a pair of atom types
     * @param type1 First atom type
     * @param type2 Second atom type
     * @return Pair of (NBFIXParams, bool) where bool indicates if NBFIX exists
     */
    std::pair<NBFIXParams, bool> get_nbfix(const std::string& type1, 
                                          const std::string& type2) const {
        auto key = ParamKeyUtils::makeTypePair(type1, type2);
        auto it = nbfix_.find(key);
        if (it == nbfix_.end()) {
            return std::make_pair(NBFIXParams{}, false);
        }
        return std::make_pair(it->second, true);
    }

    const BondParams& get_bond_params(const std::string& type1, const std::string& type2) const {
        auto key = ParamKeyUtils::makeTypePair(type1, type2);
        auto it = bond_params_.find(key);
        if (it == bond_params_.end()) {
            throw std::runtime_error("Bond parameters not found for types: " + type1 + "-" + type2);
        }
        return it->second;
    }

    const AngleParams& get_angle_params(const std::string& type1,
                                      const std::string& type2,
                                      const std::string& type3) const {
        // Try both orientations of the outer atoms while keeping the middle atom fixed
        auto key1 = std::make_tuple(type1, type2, type3);
        auto it = angle_params_.find(key1);
        if (it != angle_params_.end()) {
            return it->second;
        }

        // Try the reverse orientation
        auto key2 = std::make_tuple(type3, type2, type1);
        it = angle_params_.find(key2);
        if (it != angle_params_.end()) {
            return it->second;
        }

        throw std::runtime_error("Angle parameters not found for types: " + 
                               type1 + "-" + type2 + "-" + type3);
    }

    const std::vector<DihedralParams>& get_dihedral_params(const std::string& type1,
                                                          const std::string& type2,
                                                          const std::string& type3,
                                                          const std::string& type4) const {
        auto key = ParamKeyUtils::makeTypeQuad(type1, type2, type3, type4);
        auto it = dihedral_params_.find(key);
        if (it == dihedral_params_.end()) {
            throw std::runtime_error("Dihedral parameters not found for types: " +
                                   type1 + "-" + type2 + "-" + type3 + "-" + type4);
        }
        return it->second;
    }

    const ImproperParams& get_improper_params(const std::string& type1, const std::string& type2,
                                            const std::string& type3, const std::string& type4) const {
        auto key = ParamKeyUtils::makeTypeQuad(type1, type2, type3, type4);
        auto it = improper_params_.find(key);
        if (it == improper_params_.end()) {
            throw std::runtime_error("Improper parameters not found for types: " +
                                   type1 + "-" + type2 + "-" + type3 + "-" + type4);
        }
        return it->second;
    }

    // === Existence Check Methods ===

    bool has_atom_mass(const std::string& type) const {
        return atom_masses_.find(type) != atom_masses_.end();
    }

    bool has_lj_params(const std::string& type) const {
        return lj_params_.find(type) != lj_params_.end();
    }

    bool has_nbfix(const std::string& type1, const std::string& type2) const {
        auto key = ParamKeyUtils::makeTypePair(type1, type2);
        return nbfix_.find(key) != nbfix_.end();
    }

    bool has_bond_params(const std::string& type1, const std::string& type2) const {
        auto key = ParamKeyUtils::makeTypePair(type1, type2);
        return bond_params_.find(key) != bond_params_.end();
    }

    bool has_angle_params(const std::string& type1, const std::string& type2,
                         const std::string& type3) const {
        auto key1 = std::make_tuple(type1, type2, type3);
        auto key2 = std::make_tuple(type3, type2, type1);
        return angle_params_.find(key1) != angle_params_.end() ||
               angle_params_.find(key2) != angle_params_.end();
    }

    bool has_dihedral_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4) const {
        auto key = ParamKeyUtils::makeTypeQuad(type1, type2, type3, type4);
        return dihedral_params_.find(key) != dihedral_params_.end();
    }

    bool has_improper_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4) const {
        auto key = ParamKeyUtils::makeTypeQuad(type1, type2, type3, type4);
        return improper_params_.find(key) != improper_params_.end();
    }

    // === Size Methods ===

    size_t get_num_atom_types() const { return atom_masses_.size(); }
    size_t get_num_lj_params() const { return lj_params_.size(); }
    size_t get_num_nbfix() const { return nbfix_.size(); }
    size_t get_num_bond_types() const { return bond_params_.size(); }
    size_t get_num_angle_types() const { return angle_params_.size(); }
    size_t get_num_dihedral_types() const { return dihedral_params_.size(); }
    size_t get_num_improper_types() const { return improper_params_.size(); }

    // === Utility Methods ===

    void clear() {
        atom_masses_.clear();
        lj_params_.clear();
        nbfix_.clear();
        bond_params_.clear();
        angle_params_.clear();
        dihedral_params_.clear();
        improper_params_.clear();
        nonbonded_params_ = NonbondedParams{};
    }

    std::set<std::string> get_atom_types() const {
        std::set<std::string> types;
        for (const auto& pair : atom_masses_) {
            types.insert(pair.first);
        }
        return types;
    }

    // === Analysis Methods (from new version) ===

    ForceFieldStats get_statistics() const {
        ForceFieldStats stats;
        stats.num_atom_types = atom_masses_.size();
        stats.num_lj_params = lj_params_.size();
        stats.num_nbfix = nbfix_.size();
        stats.num_bond_types = bond_params_.size();
        stats.num_angle_types = angle_params_.size();
        stats.num_dihedral_types = dihedral_params_.size();
        stats.num_improper_types = improper_params_.size();
        return stats;
    }

    std::string get_summary() const {
        auto stats = get_statistics();
        std::ostringstream oss;
        oss << "Force Field Summary:\n";
        oss << "  Atom types: " << stats.num_atom_types << "\n";
        oss << "  LJ parameters: " << stats.num_lj_params << "\n";
        oss << "  NBFIX entries: " << stats.num_nbfix << "\n";
        oss << "  Bond types: " << stats.num_bond_types << "\n";
        oss << "  Angle types: " << stats.num_angle_types << "\n";
        oss << "  Dihedral types: " << stats.num_dihedral_types << "\n";
        oss << "  Improper types: " << stats.num_improper_types;
        return oss.str();
    }

    std::set<std::string> find_missing_lj_params() const {
        std::set<std::string> missing;
        for (const auto& pair : atom_masses_) {
            if (lj_params_.find(pair.first) == lj_params_.end()) {
                missing.insert(pair.first);
            }
        }
        return missing;
    }

    std::set<std::string> find_missing_masses() const {
        std::set<std::string> missing;
        for (const auto& pair : lj_params_) {
            if (atom_masses_.find(pair.first) == atom_masses_.end()) {
                missing.insert(pair.first);
            }
        }
        return missing;
    }

    bool validate_force_field() const {
        auto missing_lj = find_missing_lj_params();
        auto missing_masses = find_missing_masses();
        return missing_lj.empty() && missing_masses.empty();
    }

    CompletenessResult check_completeness(const std::set<std::string>& atom_types) const {
        CompletenessResult result;
        
        for (const auto& type : atom_types) {
            if (!has_atom_mass(type)) {
                result.missing_atom_masses.insert(type);
            }
            if (!has_lj_params(type)) {
                result.missing_lj_params.insert(type);
            }
        }
        
        result.is_complete = result.missing_atom_masses.empty() && result.missing_lj_params.empty();
        
        std::ostringstream oss;
        if (result.is_complete) {
            oss << "Force field is complete for all " << atom_types.size() << " atom types.";
        } else {
            oss << "Force field is incomplete. Missing masses: " << result.missing_atom_masses.size()
                << ", missing LJ params: " << result.missing_lj_params.size();
        }
        result.summary = oss.str();
        
        return result;
    }

    std::vector<std::string> validate_consistency() const {
        std::vector<std::string> errors;
        
        auto missing_lj = find_missing_lj_params();
        auto missing_masses = find_missing_masses();
        
        for (const auto& type : missing_lj) {
            errors.push_back("Missing LJ parameters for atom type: " + type);
        }
        
        for (const auto& type : missing_masses) {
            errors.push_back("Missing mass for atom type: " + type);
        }
        
        return errors;
    }

    // === Access to nonbonded parameters ===

    const NonbondedParams& get_nonbonded_params() const { return nonbonded_params_; }
    NonbondedParams& get_nonbonded_params() { return nonbonded_params_; }

    // === Backward compatibility helper methods ===

    static std::pair<std::string, std::string> makeTypePair(const std::string& type1,
                                                           const std::string& type2) {
        return ParamKeyUtils::makeTypePair(type1, type2);
    }

    static std::tuple<std::string, std::string, std::string> makeTypeTriple(
        const std::string& type1, const std::string& type2, const std::string& type3) {
        return ParamKeyUtils::makeTypeTriple(type1, type2, type3);
    }

    static std::tuple<std::string, std::string, std::string, std::string> makeTypeQuad(
        const std::string& type1, const std::string& type2,
        const std::string& type3, const std::string& type4) {
        return ParamKeyUtils::makeTypeQuad(type1, type2, type3, type4);
    }

    // === Direct access to parameter maps (for Python bindings) ===

    const std::map<std::string, double>& get_atom_masses() const { return atom_masses_; }
    
    // Python binding compatible methods
    const std::map<std::string, LJParams>& get_lj_params_map() const { return lj_params_; }
    const std::map<std::pair<std::string, std::string>, NBFIXParams>& get_nbfix_map() const { return nbfix_; }
    const std::map<std::pair<std::string, std::string>, BondParams>& get_bond_params_map() const { return bond_params_; }
    const std::map<std::tuple<std::string, std::string, std::string>, AngleParams>& get_angle_params_map() const { return angle_params_; }
    const std::map<std::tuple<std::string, std::string, std::string, std::string>, std::vector<DihedralParams>>& get_dihedral_params_map() const { return dihedral_params_; }
    const std::map<std::tuple<std::string, std::string, std::string, std::string>, ImproperParams>& get_improper_params_map() const { return improper_params_; }
    
    // Backward compatibility methods (without _map suffix)
    const std::map<std::string, LJParams>& get_lj_params() const { return lj_params_; }
    const std::map<std::pair<std::string, std::string>, NBFIXParams>& get_nbfix() const { return nbfix_; }
    const std::map<std::pair<std::string, std::string>, BondParams>& get_bond_params() const { return bond_params_; }
    const std::map<std::tuple<std::string, std::string, std::string>, AngleParams>& get_angle_params() const { return angle_params_; }
    const std::map<std::tuple<std::string, std::string, std::string, std::string>, std::vector<DihedralParams>>& get_dihedral_params() const { return dihedral_params_; }
    const std::map<std::tuple<std::string, std::string, std::string, std::string>, ImproperParams>& get_improper_params() const { return improper_params_; }

private:
    // Parameter storage
    std::map<std::string, double> atom_masses_;
    std::map<std::string, LJParams> lj_params_;
    std::map<std::pair<std::string, std::string>, NBFIXParams> nbfix_;
    std::map<std::pair<std::string, std::string>, BondParams> bond_params_;
    std::map<std::tuple<std::string, std::string, std::string>, AngleParams> angle_params_;
    std::map<std::tuple<std::string, std::string, std::string, std::string>, std::vector<DihedralParams>> dihedral_params_;
    std::map<std::tuple<std::string, std::string, std::string, std::string>, ImproperParams> improper_params_;
    NonbondedParams nonbonded_params_;
};

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_HPP 