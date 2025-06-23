#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_MANAGER_HPP
#define PYGCMC_MODEL_FORCEFIELD_MANAGER_HPP

#include "ForceFieldParams.hpp"
#include "ForceFieldInterface.hpp"
#include <stdexcept>
#include <sstream>
#include <algorithm>

namespace pygcmc {
namespace model {
namespace forcefield {

/**
 * @brief Core force field parameter manager
 * 
 * This class handles all parameter addition and retrieval operations.
 * It provides error checking and validation for all parameter types.
 * All methods perform consistency checks and provide meaningful error messages.
 */
class ForceFieldParameterManager {
public:
    /**
     * @brief Constructor
     * @param storage Reference to the force field storage container
     */
    explicit ForceFieldParameterManager(ForceFieldStorage& storage) : storage_(storage) {}
    ~ForceFieldParameterManager() = default;

    // === Parameter Addition Methods ===

    /**
     * @brief Add atom mass parameter with validation
     * @param type Atom type name (cannot be empty)
     * @param mass Atomic mass (amu, must be positive)
     */
    void add_atom_mass(const std::string& type, double mass) {
        if (type.empty()) {
            throw std::invalid_argument("Atom type cannot be empty");
        }
        if (mass <= 0.0) {
            throw std::invalid_argument("Atom mass must be positive for type: " + type);
        }
        storage_.atom_masses[type] = mass;
    }

    /**
     * @brief Add Lennard-Jones parameters for an atom type
     * @param type The atom type (cannot be empty)
     * @param epsilon Well depth (kcal/mole)
     * @param rmin_half Rmin/2: HALF of the distance at minimum energy (Angstroms)
     * 
     * Note: In CHARMM force fields, epsilon can be negative and rmin_half can be zero
     * We follow the original implementation without strict validation
     */
    void add_lj_params(const std::string& type, double epsilon, double rmin_half) {
        if (type.empty()) {
            throw std::invalid_argument("Atom type cannot be empty");
        }
        LJParams params{epsilon, rmin_half};
        storage_.lj_params[type] = params;
    }

    /**
     * @brief Add NBFIX parameters for a specific pair of atom types
     * @param type1 First atom type (cannot be empty)
     * @param type2 Second atom type (cannot be empty)
     * @param epsilon Well depth (kcal/mole)
     * @param rmin Distance at minimum energy (Angstroms)
     */
    void add_nbfix(const std::string& type1, const std::string& type2, 
                  double epsilon, double rmin) {
        if (type1.empty() || type2.empty()) {
            throw std::invalid_argument("Atom types cannot be empty for NBFIX");
        }
        auto key = ParamKeyUtils::make_pair(type1, type2);
        NBFIXParams params{epsilon, rmin};
        storage_.nbfix[key] = params;
    }

    /**
     * @brief Add bond parameters
     * @param type1 First atom type (cannot be empty)
     * @param type2 Second atom type (cannot be empty)
     * @param kb Force constant
     * @param b0 Equilibrium bond length
     */
    void add_bond_params(const std::string& type1, const std::string& type2, 
                        double kb, double b0) {
        if (type1.empty() || type2.empty()) {
            throw std::invalid_argument("Atom types cannot be empty for bond parameters");
        }
        auto key = ParamKeyUtils::make_pair(type1, type2);
        BondParams params{kb, b0};
        storage_.bond_params[key] = params;
    }

    /**
     * @brief Add angle parameters
     * @param type1 First atom type (cannot be empty)
     * @param type2 Second atom type (center atom, cannot be empty)
     * @param type3 Third atom type (cannot be empty)
     * @param ktheta Angle force constant
     * @param theta0 Equilibrium angle
     * @param kub Urey-Bradley force constant (optional)
     * @param s0 Urey-Bradley equilibrium distance (optional)
     */
    void add_angle_params(const std::string& type1, const std::string& type2, 
                         const std::string& type3, double ktheta, double theta0, 
                         double kub = 0.0, double s0 = 0.0) {
        if (type1.empty() || type2.empty() || type3.empty()) {
            throw std::invalid_argument("Atom types cannot be empty for angle parameters");
        }
        auto key = ParamKeyUtils::make_triple(type1, type2, type3);
        AngleParams params{ktheta, theta0, kub, s0};
        storage_.angle_params[key] = params;
    }

    /**
     * @brief Add dihedral parameters
     * @param type1 First atom type (cannot be empty)
     * @param type2 Second atom type (cannot be empty)
     * @param type3 Third atom type (cannot be empty)
     * @param type4 Fourth atom type (cannot be empty)
     * @param kchi Force constant
     * @param n Multiplicity
     * @param delta Phase shift
     */
    void add_dihedral_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4,
                            double kchi, int n, double delta) {
        if (type1.empty() || type2.empty() || type3.empty() || type4.empty()) {
            throw std::invalid_argument("Atom types cannot be empty for dihedral parameters");
        }
        auto key = ParamKeyUtils::make_quad(type1, type2, type3, type4);
        DihedralParams params{kchi, n, delta};
        storage_.dihedral_params[key].push_back(params);
    }

    /**
     * @brief Add improper parameters
     * @param type1 First atom type (cannot be empty)
     * @param type2 Second atom type (cannot be empty)
     * @param type3 Third atom type (cannot be empty)
     * @param type4 Fourth atom type (cannot be empty)
     * @param kpsi Force constant
     * @param psi0 Equilibrium improper angle
     */
    void add_improper_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4,
                            double kpsi, double psi0) {
        if (type1.empty() || type2.empty() || type3.empty() || type4.empty()) {
            throw std::invalid_argument("Atom types cannot be empty for improper parameters");
        }
        auto key = ParamKeyUtils::make_quad(type1, type2, type3, type4);
        ImproperParams params{kpsi, psi0};
        storage_.improper_params[key] = params;
    }

    // === Parameter Retrieval Methods ===

    /**
     * @brief Get atom mass with error checking
     * @param type Atom type name
     * @return Atomic mass (amu)
     * @throws std::runtime_error if type not found
     */
    double get_atom_mass(const std::string& type) const {
        auto it = storage_.atom_masses.find(type);
        if (it == storage_.atom_masses.end()) {
            throw std::runtime_error("Atom mass not found for type: " + type);
        }
        return it->second;
    }

    /**
     * @brief Get LJ parameters with error checking
     * @param type Atom type name
     * @return LJ parameters
     * @throws std::runtime_error if type not found
     */
    const LJParams& get_lj_params(const std::string& type) const {
        auto it = storage_.lj_params.find(type);
        if (it == storage_.lj_params.end()) {
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
        auto key = ParamKeyUtils::make_pair(type1, type2);
        auto it = storage_.nbfix.find(key);
        if (it == storage_.nbfix.end()) {
            return std::make_pair(NBFIXParams{}, false);
        }
        return std::make_pair(it->second, true);
    }

    /**
     * @brief Get bond parameters with error checking
     * @param type1 First atom type
     * @param type2 Second atom type
     * @return Bond parameters
     * @throws std::runtime_error if parameters not found
     */
    const BondParams& get_bond_params(const std::string& type1, const std::string& type2) const {
        auto key = ParamKeyUtils::make_pair(type1, type2);
        auto it = storage_.bond_params.find(key);
        if (it == storage_.bond_params.end()) {
            throw std::runtime_error("Bond parameters not found for types: " + type1 + "-" + type2);
        }
        return it->second;
    }

    /**
     * @brief Get angle parameters with error checking and bidirectional lookup
     * @param type1 First atom type
     * @param type2 Second atom type (center)
     * @param type3 Third atom type
     * @return Angle parameters
     * @throws std::runtime_error if parameters not found
     */
    const AngleParams& get_angle_params(const std::string& type1,
                                      const std::string& type2,
                                      const std::string& type3) const {
        // Try both orientations of the outer atoms while keeping the middle atom fixed
        auto key1 = std::make_tuple(type1, type2, type3);
        auto it = storage_.angle_params.find(key1);
        if (it != storage_.angle_params.end()) {
            return it->second;
        }

        // Try the reverse orientation
        auto key2 = std::make_tuple(type3, type2, type1);
        it = storage_.angle_params.find(key2);
        if (it != storage_.angle_params.end()) {
            return it->second;
        }

        throw std::runtime_error("Angle parameters not found for types: " + 
                               type1 + "-" + type2 + "-" + type3);
    }

    /**
     * @brief Get dihedral parameters with error checking
     * @param type1 First atom type
     * @param type2 Second atom type
     * @param type3 Third atom type
     * @param type4 Fourth atom type
     * @return Vector of dihedral parameters (multiple terms possible)
     * @throws std::runtime_error if parameters not found
     */
    const std::vector<DihedralParams>& get_dihedral_params(const std::string& type1,
                                                          const std::string& type2,
                                                          const std::string& type3,
                                                          const std::string& type4) const {
        auto key = ParamKeyUtils::make_quad(type1, type2, type3, type4);
        auto it = storage_.dihedral_params.find(key);
        if (it == storage_.dihedral_params.end()) {
            throw std::runtime_error("Dihedral parameters not found for types: " +
                                   type1 + "-" + type2 + "-" + type3 + "-" + type4);
        }
        return it->second;
    }

    /**
     * @brief Get improper parameters with error checking
     * @param type1 First atom type
     * @param type2 Second atom type
     * @param type3 Third atom type
     * @param type4 Fourth atom type
     * @return Improper parameters
     * @throws std::runtime_error if parameters not found
     */
    const ImproperParams& get_improper_params(const std::string& type1, const std::string& type2,
                                            const std::string& type3, const std::string& type4) const {
        auto key = ParamKeyUtils::make_quad(type1, type2, type3, type4);
        auto it = storage_.improper_params.find(key);
        if (it == storage_.improper_params.end()) {
            throw std::runtime_error("Improper parameters not found for types: " +
                                   type1 + "-" + type2 + "-" + type3 + "-" + type4);
        }
        return it->second;
    }

    /**
     * @brief Clear all parameters
     */
    void clear() {
        storage_.clear();
    }

    /**
     * @brief Get all atom types with masses defined
     * @return Set of atom types
     */
    std::set<std::string> get_atom_types() const {
        std::set<std::string> types;
        for (const auto& pair : storage_.atom_masses) {
            types.insert(pair.first);
        }
        return types;
    }

    /**
     * @brief Access to nonbonded parameters (const)
     * @return Const reference to nonbonded parameters
     */
    const NonbondedParams& get_nonbonded_params() const { 
        return storage_.nonbonded_params; 
    }
    
    /**
     * @brief Access to nonbonded parameters (mutable)
     * @return Mutable reference to nonbonded parameters
     */
    NonbondedParams& get_nonbonded_params() { 
        return storage_.nonbonded_params; 
    }

private:
    ForceFieldStorage& storage_;  ///< Reference to parameter storage
};

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_MANAGER_HPP 