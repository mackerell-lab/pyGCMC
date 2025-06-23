#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_MANAGER_HPP
#define PYGCMC_MODEL_FORCEFIELD_MANAGER_HPP

#include "ForceFieldParams.hpp"
#include "ForceFieldInterface.hpp"
#include "ForceFieldParameterAdder.hpp"
#include "ForceFieldParameterRetriever.hpp"
#include <stdexcept>
#include <sstream>
#include <algorithm>

namespace pygcmc {
namespace model {
namespace forcefield {

/**
 * @brief Lightweight force field parameter manager using composition
 * 
 * This class orchestrates parameter addition and retrieval operations
 * by delegating to specialized components. It provides a unified interface
 * while maintaining clear separation of concerns and optimal performance.
 */
class ForceFieldParameterManager {
public:
    /**
     * @brief Constructor
     * @param storage Reference to the force field storage container
     */
    explicit ForceFieldParameterManager(ForceFieldStorage& storage) : 
        storage_(storage),
        adder_(storage_),
        retriever_(storage_) {}
    ~ForceFieldParameterManager() = default;

    // === Parameter Addition Methods (Delegated) ===

    /**
     * @brief Add atom mass parameter with validation
     */
    void add_atom_mass(const std::string& type, double mass) {
        adder_.add_atom_mass(type, mass);
    }

    /**
     * @brief Add Lennard-Jones parameters for an atom type
     */
    void add_lj_params(const std::string& type, double epsilon, double rmin_half) {
        adder_.add_lj_params(type, epsilon, rmin_half);
    }

    /**
     * @brief Add NBFIX parameters for a specific pair of atom types
     */
    void add_nbfix(const std::string& type1, const std::string& type2, 
                  double epsilon, double rmin) {
        adder_.add_nbfix(type1, type2, epsilon, rmin);
    }

    /**
     * @brief Add bond parameters
     */
    void add_bond_params(const std::string& type1, const std::string& type2, 
                        double kb, double b0) {
        adder_.add_bond_params(type1, type2, kb, b0);
    }

    /**
     * @brief Add angle parameters
     */
    void add_angle_params(const std::string& type1, const std::string& type2, 
                         const std::string& type3, double ktheta, double theta0, 
                         double kub = 0.0, double s0 = 0.0) {
        adder_.add_angle_params(type1, type2, type3, ktheta, theta0, kub, s0);
    }

    /**
     * @brief Add dihedral parameters
     */
    void add_dihedral_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4,
                            double kchi, int n, double delta) {
        adder_.add_dihedral_params(type1, type2, type3, type4, kchi, n, delta);
    }

    /**
     * @brief Add improper parameters
     */
    void add_improper_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4,
                            double kpsi, double psi0) {
        adder_.add_improper_params(type1, type2, type3, type4, kpsi, psi0);
    }

    // === Parameter Retrieval Methods (Delegated) ===

    /**
     * @brief Get atom mass with enhanced error reporting
     */
    double get_atom_mass(const std::string& type) const {
        return retriever_.get_atom_mass(type);
    }

    /**
     * @brief Get LJ parameters with enhanced error reporting
     */
    const LJParams& get_lj_params(const std::string& type) const {
        return retriever_.get_lj_params(type);
    }

    /**
     * @brief Get NBFIX parameters for a pair of atom types
     */
    std::pair<NBFIXParams, bool> get_nbfix(const std::string& type1, 
                                          const std::string& type2) const {
        return retriever_.get_nbfix(type1, type2);
    }

    /**
     * @brief Get bond parameters with enhanced error reporting
     */
    const BondParams& get_bond_params(const std::string& type1, const std::string& type2) const {
        return retriever_.get_bond_params(type1, type2);
    }

    /**
     * @brief Get angle parameters with bidirectional lookup
     */
    const AngleParams& get_angle_params(const std::string& type1,
                                      const std::string& type2,
                                      const std::string& type3) const {
        return retriever_.get_angle_params(type1, type2, type3);
    }

    /**
     * @brief Get dihedral parameters with enhanced error reporting
     */
    const std::vector<DihedralParams>& get_dihedral_params(const std::string& type1,
                                                          const std::string& type2,
                                                          const std::string& type3,
                                                          const std::string& type4) const {
        return retriever_.get_dihedral_params(type1, type2, type3, type4);
    }

    /**
     * @brief Get improper parameters with enhanced error reporting
     */
    const ImproperParams& get_improper_params(const std::string& type1, const std::string& type2,
                                            const std::string& type3, const std::string& type4) const {
        return retriever_.get_improper_params(type1, type2, type3, type4);
    }

    // === Additional Methods ===

    /**
     * @brief Add multiple parameters at once using batch operations
     */
    template<typename Container>
    void add_atom_masses_batch(const Container& mass_data) {
        adder_.add_atom_masses_batch(mass_data);
    }

    template<typename Container>
    void add_lj_params_batch(const Container& lj_data) {
        adder_.add_lj_params_batch(lj_data);
    }

    /**
     * @brief Get multiple parameters at once with error collection
     */
    std::map<std::string, double> get_atom_masses_batch(const std::vector<std::string>& types) const {
        return retriever_.get_atom_masses_batch(types);
    }

    std::map<std::string, LJParams> get_lj_params_batch(const std::vector<std::string>& types) const {
        return retriever_.get_lj_params_batch(types);
    }

    /**
     * @brief Get parameter suggestions for missing types
     */
    std::vector<std::string> suggest_similar_types(const std::string& target_type) const {
        return retriever_.suggest_similar_types(target_type);
    }

    /**
     * @brief Clear all parameters
     */
    void clear() {
        storage_.clear();
    }

    /**
     * @brief Get all atom types with masses defined
     */
    std::set<std::string> get_atom_types() const {
        std::set<std::string> types;
        for (const auto& pair : storage_.atom_masses) {
            types.insert(pair.first);
        }
        return types;
    }

    /**
     * @brief Access to nonbonded parameters
     */
    const NonbondedParams& get_nonbonded_params() const { 
        return storage_.nonbonded_params; 
    }
    
    NonbondedParams& get_nonbonded_params() { 
        return storage_.nonbonded_params; 
    }

private:
    ForceFieldStorage& storage_;              ///< Reference to parameter storage
    ForceFieldParameterAdder adder_;          ///< Specialized parameter adder
    ForceFieldParameterRetriever retriever_;  ///< Specialized parameter retriever
};

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_MANAGER_HPP 