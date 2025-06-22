#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_MAIN_HPP
#define PYGCMC_MODEL_FORCEFIELD_MAIN_HPP

#include "ForceFieldParams.hpp"
#include "ForceFieldCore.hpp"
#include "ForceFieldUtils.hpp"
#include "../common/ModelInterface.hpp"

namespace pygcmc {
namespace model {
namespace forcefield {

/**
 * @brief Main force field class that holds all force field parameters
 * 
 * This class maintains backward compatibility with the original forcefield.hpp
 * while providing a modular and extensible architecture. It uses composition
 * with specialized managers for better organization and maintainability.
 * 
 * The class provides a complete interface for:
 * - Parameter management (adding, retrieving, validating parameters)
 * - Analysis and statistics (completeness checks, summaries)
 * - Backward compatibility (all original APIs preserved)
 * 
 * **Architecture:**
 * - ForceFieldStorage: Core data container
 * - ForceFieldManager: Parameter management operations
 * - ForceFieldUtils: Analysis and utility functions
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
class ForceField : public common::IValidatable {
public:
    ForceField() : manager_(storage_), utils_(storage_) {}

    // === Copy/Move constructors to ensure manager_ / utils_ correctly reference storage_ ===
    // Note: Assignment operators are deleted because manager_/utils_ contain reference members
    ForceField(const ForceField& other) : storage_(other.storage_), manager_(storage_), utils_(storage_) {}

    ForceField(ForceField&& other) noexcept : storage_(std::move(other.storage_)),
                                              manager_(storage_), utils_(storage_) {}

    // Delete assignment operators (cannot reassign reference members)
    ForceField& operator=(const ForceField&) = delete;
    ForceField& operator=(ForceField&&) = delete;

    ~ForceField() = default;

    // IValidatable interface
    bool is_valid() const override {
        return utils_.validate_force_field();
    }

    // === Parameter Addition Methods ===

    /**
     * @brief Add atom mass parameter
     * @param type Atom type name
     * @param mass Atomic mass (amu)
     */
    void add_atom_mass(const std::string& type, double mass) {
        manager_.add_atom_mass(type, mass);
    }

    /**
     * @brief Add Lennard-Jones parameters for an atom type
     * @param type The atom type
     * @param epsilon Well depth (kcal/mole)
     * @param rmin_half Rmin/2: HALF of the distance at minimum energy (Angstroms)
     */
    void add_lj_params(const std::string& type, double epsilon, double rmin_half) {
        manager_.add_lj_params(type, epsilon, rmin_half);
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
        manager_.add_nbfix(type1, type2, epsilon, rmin);
    }

    /**
     * @brief Add bond parameters
     */
    void add_bond_params(const std::string& type1, const std::string& type2, double kb, double b0) {
        manager_.add_bond_params(type1, type2, kb, b0);
    }

    /**
     * @brief Add angle parameters
     */
    void add_angle_params(const std::string& type1, const std::string& type2, const std::string& type3,
                         double ktheta, double theta0, double kub = 0.0, double s0 = 0.0) {
        manager_.add_angle_params(type1, type2, type3, ktheta, theta0, kub, s0);
    }

    /**
     * @brief Add dihedral parameters
     */
    void add_dihedral_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4,
                            double kchi, int n, double delta) {
        manager_.add_dihedral_params(type1, type2, type3, type4, kchi, n, delta);
    }

    /**
     * @brief Add improper parameters
     */
    void add_improper_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4,
                            double kpsi, double psi0) {
        manager_.add_improper_params(type1, type2, type3, type4, kpsi, psi0);
    }

    // === Parameter Retrieval Methods ===

    /**
     * @brief Get atom mass with error checking
     */
    double get_atom_mass(const std::string& type) const {
        return manager_.get_atom_mass(type);
    }

    /**
     * @brief Get LJ parameters with error checking
     */
    const LJParams& get_lj_params(const std::string& type) const {
        return manager_.get_lj_params(type);
    }

    /**
     * @brief Get NBFIX parameters for a pair of atom types
     * @param type1 First atom type
     * @param type2 Second atom type
     * @return Pair of (NBFIXParams, bool) where bool indicates if NBFIX exists
     */
    std::pair<NBFIXParams, bool> get_nbfix(const std::string& type1, 
                                          const std::string& type2) const {
        return manager_.get_nbfix(type1, type2);
    }

    /**
     * @brief Get bond parameters with error checking
     */
    const BondParams& get_bond_params(const std::string& type1, const std::string& type2) const {
        return manager_.get_bond_params(type1, type2);
    }

    /**
     * @brief Get angle parameters with error checking and bidirectional lookup
     */
    const AngleParams& get_angle_params(const std::string& type1,
                                      const std::string& type2,
                                      const std::string& type3) const {
        return manager_.get_angle_params(type1, type2, type3);
    }

    /**
     * @brief Get dihedral parameters with error checking
     */
    const std::vector<DihedralParams>& get_dihedral_params(const std::string& type1,
                                                          const std::string& type2,
                                                          const std::string& type3,
                                                          const std::string& type4) const {
        return manager_.get_dihedral_params(type1, type2, type3, type4);
    }

    /**
     * @brief Get improper parameters with error checking
     */
    const ImproperParams& get_improper_params(const std::string& type1, const std::string& type2,
                                            const std::string& type3, const std::string& type4) const {
        return manager_.get_improper_params(type1, type2, type3, type4);
    }

    // === Existence Check Methods ===

    bool has_atom_mass(const std::string& type) const {
        return manager_.has_atom_mass(type);
    }

    bool has_lj_params(const std::string& type) const {
        return manager_.has_lj_params(type);
    }

    bool has_nbfix(const std::string& type1, const std::string& type2) const {
        return manager_.has_nbfix(type1, type2);
    }

    bool has_bond_params(const std::string& type1, const std::string& type2) const {
        return manager_.has_bond_params(type1, type2);
    }

    bool has_angle_params(const std::string& type1, const std::string& type2,
                         const std::string& type3) const {
        return manager_.has_angle_params(type1, type2, type3);
    }

    bool has_dihedral_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4) const {
        return manager_.has_dihedral_params(type1, type2, type3, type4);
    }

    bool has_improper_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4) const {
        return manager_.has_improper_params(type1, type2, type3, type4);
    }

    // === Size Methods ===
    size_t get_num_atom_types() const { return manager_.get_num_atom_types(); }
    size_t get_num_lj_params() const { return manager_.get_num_lj_params(); }
    size_t get_num_nbfix() const { return manager_.get_num_nbfix(); }
    size_t get_num_bond_types() const { return manager_.get_num_bond_types(); }
    size_t get_num_angle_types() const { return manager_.get_num_angle_types(); }
    size_t get_num_dihedral_types() const { return manager_.get_num_dihedral_types(); }
    size_t get_num_improper_types() const { return manager_.get_num_improper_types(); }

    // === Access to Nonbonded Parameters ===
    const NonbondedParams& get_nonbonded_params() const { 
        return storage_.nonbonded_params; 
    }
    
    NonbondedParams& get_nonbonded_params() { 
        return storage_.nonbonded_params; 
    }

    // === Static Helper Methods for Making Parameter Keys ===
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

    // === Direct Access to Parameter Maps (for Python bindings) ===
    const std::map<std::string, double>& get_atom_masses() const { 
        return utils_.get_atom_masses(); 
    }
    
    const std::map<std::string, LJParams>& get_lj_params() const { 
        return utils_.get_lj_params(); 
    }
    
    const std::map<std::pair<std::string, std::string>, NBFIXParams>& get_nbfix() const { 
        return utils_.get_nbfix(); 
    }
    
    const std::map<std::pair<std::string, std::string>, BondParams>& get_bond_params() const { 
        return utils_.get_bond_params(); 
    }
    
    const std::map<std::tuple<std::string, std::string, std::string>, AngleParams>& get_angle_params() const { 
        return utils_.get_angle_params(); 
    }
    
    const std::map<std::tuple<std::string, std::string, std::string, std::string>, std::vector<DihedralParams>>& get_dihedral_params() const { 
        return utils_.get_dihedral_params(); 
    }
    
    const std::map<std::tuple<std::string, std::string, std::string, std::string>, ImproperParams>& get_improper_params() const { 
        return utils_.get_improper_params(); 
    }

    // === Additional Utility Methods ===

    /**
     * @brief Get all atom types
     */
    std::set<std::string> get_atom_types() const {
        return manager_.get_atom_types();
    }

    /**
     * @brief Clear all parameters
     */
    void clear() {
        manager_.clear();
    }

    // === String Representation ===
    std::string to_string() const {
        return utils_.get_summary();
    }

    // === Statistics and Analysis ===

    /**
     * @brief Get a human-readable summary of the force field
     */
    std::string get_summary() const {
        return utils_.get_summary();
    }

    /**
     * @brief Get comprehensive statistics about the force field
     */
    auto get_statistics() const {
        return utils_.get_statistics();
    }

    /**
     * @brief Get detailed statistics with parameter information
     */
    std::string get_detailed_statistics() const {
        return utils_.get_detailed_statistics();
    }

    /**
     * @brief Check completeness for a given set of atom types
     */
    auto check_completeness(const std::set<std::string>& atom_types) const {
        return utils_.check_completeness(atom_types);
    }

    /**
     * @brief Find all atom types that have masses but no LJ parameters
     */
    std::set<std::string> find_missing_lj_params() const {
        return utils_.find_missing_lj_params();
    }

    /**
     * @brief Find all atom types that have LJ parameters but no masses
     */
    std::set<std::string> find_missing_masses() const {
        return utils_.find_missing_masses();
    }

    /**
     * @brief Get all unique atom types referenced in the force field
     */
    std::set<std::string> get_all_referenced_types() const {
        return utils_.get_all_referenced_types();
    }

    /**
     * @brief Check if force field parameters are consistent
     */
    std::vector<std::string> validate_consistency() const {
        return utils_.validate_consistency();
    }

private:
    ForceFieldStorage storage_;   ///< Core data storage
    ForceFieldManager manager_;   ///< Parameter management
    ForceFieldUtils utils_;       ///< Analysis and utilities
};

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_MAIN_HPP 