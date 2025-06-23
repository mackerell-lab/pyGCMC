#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_CORE_HPP
#define PYGCMC_MODEL_FORCEFIELD_CORE_HPP

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
    // === Constructors and Destructor ===
    
    ForceField() : 
        manager_(storage_), 
        operations_(storage_), 
        analyzer_(storage_), 
        validator_(storage_) {}

    // Copy constructor
    ForceField(const ForceField& other) : 
        storage_(other.storage_), 
        manager_(storage_), 
        operations_(storage_), 
        analyzer_(storage_), 
        validator_(storage_) {}

    // Move constructor
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

    // === Core Interface Implementation ===

    bool is_valid() const override {
        return validator_.validate_force_field();
    }

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

    // === Parameter Operations Interface (implemented in separate files) ===

    // Parameter addition methods
    void add_atom_mass(const std::string& type, double mass) override;
    void add_lj_params(const std::string& type, double epsilon, double rmin_half) override;
    void add_nbfix(const std::string& type1, const std::string& type2, 
                  double epsilon, double rmin) override;
    void add_bond_params(const std::string& type1, const std::string& type2, 
                        double kb, double b0) override;
    void add_angle_params(const std::string& type1, const std::string& type2, 
                         const std::string& type3, double ktheta, double theta0, 
                         double kub = 0.0, double s0 = 0.0) override;
    void add_dihedral_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4,
                            double kchi, int n, double delta) override;
    void add_improper_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4,
                            double kpsi, double psi0) override;

    // Parameter retrieval methods
    double get_atom_mass(const std::string& type) const override;
    const LJParams& get_lj_params(const std::string& type) const override;
    std::pair<NBFIXParams, bool> get_nbfix(const std::string& type1, 
                                          const std::string& type2) const override;
    const BondParams& get_bond_params(const std::string& type1, 
                                     const std::string& type2) const override;
    const AngleParams& get_angle_params(const std::string& type1,
                                       const std::string& type2,
                                       const std::string& type3) const override;
    const std::vector<DihedralParams>& get_dihedral_params(const std::string& type1,
                                                          const std::string& type2,
                                                          const std::string& type3,
                                                          const std::string& type4) const override;
    const ImproperParams& get_improper_params(const std::string& type1, 
                                             const std::string& type2,
                                             const std::string& type3, 
                                             const std::string& type4) const override;

    // Parameter checking methods
    bool has_atom_mass(const std::string& type) const override;
    bool has_lj_params(const std::string& type) const override;
    bool has_nbfix(const std::string& type1, const std::string& type2) const override;
    bool has_bond_params(const std::string& type1, const std::string& type2) const override;
    bool has_angle_params(const std::string& type1, const std::string& type2,
                         const std::string& type3) const override;
    bool has_dihedral_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4) const override;
    bool has_improper_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4) const override;

    // Size query methods
    size_t get_num_atom_types() const override;
    size_t get_num_lj_params() const override;
    size_t get_num_nbfix() const override;
    size_t get_num_bond_types() const override;
    size_t get_num_angle_types() const override;
    size_t get_num_dihedral_types() const override;
    size_t get_num_improper_types() const override;

    // Analysis and validation methods
    bool validate_force_field() const override;
    CompletenessResult check_completeness(const std::set<std::string>& atom_types) const override;
    std::vector<std::string> validate_consistency() const override;
    ForceFieldStats get_statistics() const override;
    std::string get_summary() const override;
    std::string get_detailed_statistics() const override;
    std::set<std::string> find_missing_lj_params() const override;
    std::set<std::string> find_missing_masses() const override;
    std::set<std::string> get_all_referenced_types() const override;

    // Direct map access methods
    const std::map<std::string, double>& get_atom_masses() const override;
    const std::map<std::string, LJParams>& get_lj_params_map() const override;
    const std::map<std::pair<std::string, std::string>, NBFIXParams>& get_nbfix_map() const override;
    const std::map<std::pair<std::string, std::string>, BondParams>& get_bond_params_map() const override;
    const std::map<std::tuple<std::string, std::string, std::string>, AngleParams>& get_angle_params_map() const override;
    const std::map<std::tuple<std::string, std::string, std::string, std::string>, std::vector<DihedralParams>>& get_dihedral_params_map() const override;
    const std::map<std::tuple<std::string, std::string, std::string, std::string>, ImproperParams>& get_improper_params_map() const override;

protected:
    // === Core Components ===
    ForceFieldStorage storage_;                 ///< Core data storage
    ForceFieldParameterManager manager_;        ///< Parameter management
    ForceFieldOperations operations_;           ///< Existence checks and size queries
    ForceFieldAnalyzer analyzer_;               ///< Statistics and analysis
    ForceFieldValidator validator_;             ///< Validation and completeness
};

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_CORE_HPP 