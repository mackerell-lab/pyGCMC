#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_INTERFACE_HPP
#define PYGCMC_MODEL_FORCEFIELD_INTERFACE_HPP

#include "ForceFieldParams.hpp"
#include "../common/ModelInterface.hpp"
#include <set>
#include <vector>

namespace pygcmc {
namespace model {
namespace forcefield {

/**
 * @brief Abstract interface for force field parameter operations
 * 
 * This interface defines the core operations that any force field
 * implementation must provide. It separates the interface from 
 * implementation details for better modularity.
 */
class IForceFieldOperations {
public:
    virtual ~IForceFieldOperations() = default;

    // === Parameter Addition Interface ===
    virtual void add_atom_mass(const std::string& type, double mass) = 0;
    virtual void add_lj_params(const std::string& type, double epsilon, double rmin_half) = 0;
    virtual void add_nbfix(const std::string& type1, const std::string& type2, 
                          double epsilon, double rmin) = 0;
    virtual void add_bond_params(const std::string& type1, const std::string& type2, 
                                double kb, double b0) = 0;
    virtual void add_angle_params(const std::string& type1, const std::string& type2, 
                                 const std::string& type3, double ktheta, double theta0, 
                                 double kub = 0.0, double s0 = 0.0) = 0;
    virtual void add_dihedral_params(const std::string& type1, const std::string& type2,
                                    const std::string& type3, const std::string& type4,
                                    double kchi, int n, double delta) = 0;
    virtual void add_improper_params(const std::string& type1, const std::string& type2,
                                    const std::string& type3, const std::string& type4,
                                    double kpsi, double psi0) = 0;

    // === Parameter Retrieval Interface ===
    virtual double get_atom_mass(const std::string& type) const = 0;
    virtual const LJParams& get_lj_params(const std::string& type) const = 0;
    virtual std::pair<NBFIXParams, bool> get_nbfix(const std::string& type1, 
                                                   const std::string& type2) const = 0;
    virtual const BondParams& get_bond_params(const std::string& type1, 
                                             const std::string& type2) const = 0;
    virtual const AngleParams& get_angle_params(const std::string& type1,
                                               const std::string& type2,
                                               const std::string& type3) const = 0;
    virtual const std::vector<DihedralParams>& get_dihedral_params(const std::string& type1,
                                                                  const std::string& type2,
                                                                  const std::string& type3,
                                                                  const std::string& type4) const = 0;
    virtual const ImproperParams& get_improper_params(const std::string& type1, 
                                                      const std::string& type2,
                                                      const std::string& type3, 
                                                      const std::string& type4) const = 0;

    // === Utility Interface ===
    virtual void clear() = 0;
    virtual std::set<std::string> get_atom_types() const = 0;
    virtual const NonbondedParams& get_nonbonded_params() const = 0;
    virtual NonbondedParams& get_nonbonded_params() = 0;
};

/**
 * @brief Abstract interface for force field analysis operations
 * 
 * This interface defines analysis and validation operations
 * that can be performed on force field data.
 */
class IForceFieldAnalysis {
public:
    virtual ~IForceFieldAnalysis() = default;

    // === Validation Interface ===
    virtual bool validate_force_field() const = 0;
    virtual CompletenessResult check_completeness(const std::set<std::string>& atom_types) const = 0;
    virtual std::vector<std::string> validate_consistency() const = 0;

    // === Statistics Interface ===
    virtual ForceFieldStats get_statistics() const = 0;
    virtual std::string get_summary() const = 0;
    virtual std::string get_detailed_statistics() const = 0;

    // === Analysis Interface ===
    virtual std::set<std::string> find_missing_lj_params() const = 0;
    virtual std::set<std::string> find_missing_masses() const = 0;
    virtual std::set<std::string> get_all_referenced_types() const = 0;
};

/**
 * @brief Abstract interface for checking parameter existence
 * 
 * This interface provides methods to check if specific
 * parameters exist without retrieving them.
 */
class IForceFieldChecker {
public:
    virtual ~IForceFieldChecker() = default;

    // === Existence Check Interface ===
    virtual bool has_atom_mass(const std::string& type) const = 0;
    virtual bool has_lj_params(const std::string& type) const = 0;
    virtual bool has_nbfix(const std::string& type1, const std::string& type2) const = 0;
    virtual bool has_bond_params(const std::string& type1, const std::string& type2) const = 0;
    virtual bool has_angle_params(const std::string& type1, const std::string& type2,
                                 const std::string& type3) const = 0;
    virtual bool has_dihedral_params(const std::string& type1, const std::string& type2,
                                    const std::string& type3, const std::string& type4) const = 0;
    virtual bool has_improper_params(const std::string& type1, const std::string& type2,
                                    const std::string& type3, const std::string& type4) const = 0;

    // === Size Interface ===
    virtual size_t get_num_atom_types() const = 0;
    virtual size_t get_num_lj_params() const = 0;
    virtual size_t get_num_nbfix() const = 0;
    virtual size_t get_num_bond_types() const = 0;
    virtual size_t get_num_angle_types() const = 0;
    virtual size_t get_num_dihedral_types() const = 0;
    virtual size_t get_num_improper_types() const = 0;
};

/**
 * @brief Complete force field interface combining all operations
 * 
 * This interface brings together all force field operations
 * for a complete implementation.
 */
class IForceField : public IForceFieldOperations, 
                   public IForceFieldAnalysis,
                   public IForceFieldChecker,
                   public common::IValidatable {
public:
    virtual ~IForceField() = default;

    // === Direct Map Access Interface (for Python bindings) ===
    virtual const std::map<std::string, double>& get_atom_masses() const = 0;
    virtual const std::map<std::string, LJParams>& get_lj_params_map() const = 0;
    virtual const std::map<std::pair<std::string, std::string>, NBFIXParams>& get_nbfix_map() const = 0;
    virtual const std::map<std::pair<std::string, std::string>, BondParams>& get_bond_params_map() const = 0;
    virtual const std::map<std::tuple<std::string, std::string, std::string>, AngleParams>& get_angle_params_map() const = 0;
    virtual const std::map<std::tuple<std::string, std::string, std::string, std::string>, std::vector<DihedralParams>>& get_dihedral_params_map() const = 0;
    virtual const std::map<std::tuple<std::string, std::string, std::string, std::string>, ImproperParams>& get_improper_params_map() const = 0;

    // === String Representation ===
    virtual std::string to_string() const = 0;
};

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_INTERFACE_HPP 