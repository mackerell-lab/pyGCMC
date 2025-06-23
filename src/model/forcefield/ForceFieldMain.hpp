#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_MAIN_HPP
#define PYGCMC_MODEL_FORCEFIELD_MAIN_HPP

#include "ForceFieldParams.hpp"
#include "ForceFieldInterface.hpp"
#include "ForceFieldManager.hpp"
#include "ForceFieldOperations.hpp"
#include "ForceFieldAnalysis.hpp"
#include "ForceFieldValidation.hpp"
#include "../common/ModelInterface.hpp"

namespace pygcmc {
namespace model {
namespace forcefield {

/**
 * @brief Main force field class with modular architecture
 * 
 * This class provides a complete force field interface using composition
 * of specialized components for better maintainability and AI friendliness.
 * Each component has a specific responsibility and is under 300 lines.
 * 
 * **Components:**
 * - ForceFieldParameterManager: Parameter addition and retrieval
 * - ForceFieldOperations: Existence checks and size queries
 * - ForceFieldAnalyzer: Statistics and analysis
 * - ForceFieldValidator: Validation and completeness checks
 * 
 * **Usage:**
 * ```cpp
 * ForceField ff;
 * ff.add_atom_mass("CA", 12.011);
 * ff.add_lj_params("CA", 0.070, 1.992);
 * auto mass = ff.get_atom_mass("CA");
 * auto stats = ff.get_summary();
 * ```
 */
class ForceField : public IForceField {
public:
    ForceField() : 
        manager_(storage_), 
        operations_(storage_), 
        analyzer_(storage_), 
        validator_(storage_) {}

    // === Copy/Move constructors ===
    ForceField(const ForceField& other) : 
        storage_(other.storage_), 
        manager_(storage_), 
        operations_(storage_), 
        analyzer_(storage_), 
        validator_(storage_) {}

    ForceField(ForceField&& other) noexcept : 
        storage_(std::move(other.storage_)),
        manager_(storage_), 
        operations_(storage_), 
        analyzer_(storage_), 
        validator_(storage_) {}

    // Delete assignment operators (cannot reassign reference members)
    ForceField& operator=(const ForceField&) = delete;
    ForceField& operator=(ForceField&&) = delete;

    ~ForceField() = default;

    // === IValidatable interface ===
    bool is_valid() const override {
        return validator_.validate_force_field();
    }

    // === IForceFieldOperations interface implementation ===

    void add_atom_mass(const std::string& type, double mass) override {
        manager_.add_atom_mass(type, mass);
    }

    void add_lj_params(const std::string& type, double epsilon, double rmin_half) override {
        manager_.add_lj_params(type, epsilon, rmin_half);
    }

    void add_nbfix(const std::string& type1, const std::string& type2, 
                  double epsilon, double rmin) override {
        manager_.add_nbfix(type1, type2, epsilon, rmin);
    }

    void add_bond_params(const std::string& type1, const std::string& type2, 
                        double kb, double b0) override {
        manager_.add_bond_params(type1, type2, kb, b0);
    }

    void add_angle_params(const std::string& type1, const std::string& type2, 
                         const std::string& type3, double ktheta, double theta0, 
                         double kub = 0.0, double s0 = 0.0) override {
        manager_.add_angle_params(type1, type2, type3, ktheta, theta0, kub, s0);
    }

    void add_dihedral_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4,
                            double kchi, int n, double delta) override {
        manager_.add_dihedral_params(type1, type2, type3, type4, kchi, n, delta);
    }

    void add_improper_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4,
                            double kpsi, double psi0) override {
        manager_.add_improper_params(type1, type2, type3, type4, kpsi, psi0);
    }

    double get_atom_mass(const std::string& type) const override {
        return manager_.get_atom_mass(type);
    }

    const LJParams& get_lj_params(const std::string& type) const override {
        return manager_.get_lj_params(type);
    }

    std::pair<NBFIXParams, bool> get_nbfix(const std::string& type1, 
                                          const std::string& type2) const override {
        return manager_.get_nbfix(type1, type2);
    }

    const BondParams& get_bond_params(const std::string& type1, 
                                     const std::string& type2) const override {
        return manager_.get_bond_params(type1, type2);
    }

    const AngleParams& get_angle_params(const std::string& type1,
                                       const std::string& type2,
                                       const std::string& type3) const override {
        return manager_.get_angle_params(type1, type2, type3);
    }

    const std::vector<DihedralParams>& get_dihedral_params(const std::string& type1,
                                                          const std::string& type2,
                                                          const std::string& type3,
                                                          const std::string& type4) const override {
        return manager_.get_dihedral_params(type1, type2, type3, type4);
    }

    const ImproperParams& get_improper_params(const std::string& type1, 
                                             const std::string& type2,
                                             const std::string& type3, 
                                             const std::string& type4) const override {
        return manager_.get_improper_params(type1, type2, type3, type4);
    }

    // === IForceFieldChecker interface implementation ===

    bool has_atom_mass(const std::string& type) const override {
        return operations_.has_atom_mass(type);
    }

    bool has_lj_params(const std::string& type) const override {
        return operations_.has_lj_params(type);
    }

    bool has_nbfix(const std::string& type1, const std::string& type2) const override {
        return operations_.has_nbfix(type1, type2);
    }

    bool has_bond_params(const std::string& type1, const std::string& type2) const override {
        return operations_.has_bond_params(type1, type2);
    }

    bool has_angle_params(const std::string& type1, const std::string& type2,
                         const std::string& type3) const override {
        return operations_.has_angle_params(type1, type2, type3);
    }

    bool has_dihedral_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4) const override {
        return operations_.has_dihedral_params(type1, type2, type3, type4);
    }

    bool has_improper_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4) const override {
        return operations_.has_improper_params(type1, type2, type3, type4);
    }

    size_t get_num_atom_types() const override { return operations_.get_num_atom_types(); }
    size_t get_num_lj_params() const override { return operations_.get_num_lj_params(); }
    size_t get_num_nbfix() const override { return operations_.get_num_nbfix(); }
    size_t get_num_bond_types() const override { return operations_.get_num_bond_types(); }
    size_t get_num_angle_types() const override { return operations_.get_num_angle_types(); }
    size_t get_num_dihedral_types() const override { return operations_.get_num_dihedral_types(); }
    size_t get_num_improper_types() const override { return operations_.get_num_improper_types(); }

    // === IForceFieldOperations interface implementation ===

    void clear() override {
        manager_.clear();
    }

    std::set<std::string> get_atom_types() const override {
        return manager_.get_atom_types();
    }

    const NonbondedParams& get_nonbonded_params() const override {
        return manager_.get_nonbonded_params();
    }

    NonbondedParams& get_nonbonded_params() override {
        return manager_.get_nonbonded_params();
    }

    // === IForceFieldAnalysis interface implementation ===

    bool validate_force_field() const override {
        return validator_.validate_force_field();
    }

    CompletenessResult check_completeness(const std::set<std::string>& atom_types) const override {
        return validator_.check_completeness(atom_types);
    }

    std::vector<std::string> validate_consistency() const override {
        return validator_.validate_consistency();
    }

    ForceFieldStats get_statistics() const override {
        return analyzer_.get_statistics();
    }

    std::string get_summary() const override {
        return analyzer_.get_summary();
    }

    std::string get_detailed_statistics() const override {
        return analyzer_.get_detailed_statistics();
    }

    std::set<std::string> find_missing_lj_params() const override {
        return analyzer_.find_missing_lj_params();
    }

    std::set<std::string> find_missing_masses() const override {
        return analyzer_.find_missing_masses();
    }

    std::set<std::string> get_all_referenced_types() const override {
        return analyzer_.get_all_referenced_types();
    }

    // === IForceField interface implementation ===

    const std::map<std::string, double>& get_atom_masses() const override {
        return analyzer_.get_atom_masses();
    }

    const std::map<std::string, LJParams>& get_lj_params_map() const override {
        return analyzer_.get_lj_params();
    }

    const std::map<std::pair<std::string, std::string>, NBFIXParams>& get_nbfix_map() const override {
        return analyzer_.get_nbfix();
    }

    const std::map<std::pair<std::string, std::string>, BondParams>& get_bond_params_map() const override {
        return analyzer_.get_bond_params();
    }

    const std::map<std::tuple<std::string, std::string, std::string>, AngleParams>& get_angle_params_map() const override {
        return analyzer_.get_angle_params();
    }

    const std::map<std::tuple<std::string, std::string, std::string, std::string>, std::vector<DihedralParams>>& get_dihedral_params_map() const override {
        return analyzer_.get_dihedral_params();
    }

    const std::map<std::tuple<std::string, std::string, std::string, std::string>, ImproperParams>& get_improper_params_map() const override {
        return analyzer_.get_improper_params();
    }

    std::string to_string() const override {
        return get_summary();
    }

    // === Backward Compatibility Helper Methods ===

    /**
     * @brief Static helper for making parameter keys (backward compatibility)
     */
    static std::pair<std::string, std::string> makeTypePair(const std::string& type1,
                                                           const std::string& type2) {
        return ParamKeyUtils::make_pair(type1, type2);
    }

    static std::tuple<std::string, std::string, std::string> makeTypeTriple(
        const std::string& type1, const std::string& type2, const std::string& type3) {
        return ParamKeyUtils::make_triple(type1, type2, type3);
    }

    static std::tuple<std::string, std::string, std::string, std::string> makeTypeQuad(
        const std::string& type1, const std::string& type2,
        const std::string& type3, const std::string& type4) {
        return ParamKeyUtils::make_quad(type1, type2, type3, type4);
    }

private:
    ForceFieldStorage storage_;                 ///< Core data storage
    ForceFieldParameterManager manager_;        ///< Parameter management
    ForceFieldOperations operations_;           ///< Existence checks and size queries
    ForceFieldAnalyzer analyzer_;               ///< Statistics and analysis
    ForceFieldValidator validator_;             ///< Validation and completeness
};

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_MAIN_HPP 