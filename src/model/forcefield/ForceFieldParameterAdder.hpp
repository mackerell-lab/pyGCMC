#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_PARAMETER_ADDER_HPP
#define PYGCMC_MODEL_FORCEFIELD_PARAMETER_ADDER_HPP

#include "ForceFieldParams.hpp"
#include <stdexcept>

namespace pygcmc {
namespace model {
namespace forcefield {

/**
 * @brief Template-based parameter adder for type safety and code reuse
 * 
 * This class provides templated methods for adding different types of
 * force field parameters with consistent validation and error handling.
 */
class ForceFieldParameterAdder {
public:
    explicit ForceFieldParameterAdder(ForceFieldStorage& storage) : storage_(storage) {}
    ~ForceFieldParameterAdder() = default;

    // === Single Parameter Addition ===

    /**
     * @brief Add atom mass parameter with validation
     */
    void add_atom_mass(const std::string& type, double mass) {
        validate_non_empty_type(type, "atom mass");
        if (mass <= 0.0) {
            throw std::invalid_argument("Atom mass must be positive for type: " + type);
        }
        storage_.atom_masses[type] = mass;
    }

    /**
     * @brief Add Lennard-Jones parameters
     */
    void add_lj_params(const std::string& type, double epsilon, double rmin_half) {
        validate_non_empty_type(type, "LJ parameters");
        storage_.lj_params[type] = LJParams{epsilon, rmin_half};
    }

    /**
     * @brief Add NBFIX parameters for specific atom type pairs
     */
    void add_nbfix(const std::string& type1, const std::string& type2, 
                  double epsilon, double rmin) {
        validate_pair_types(type1, type2, "NBFIX");
        auto key = ParamKeyUtils::make_pair(type1, type2);
        storage_.nbfix[key] = NBFIXParams{epsilon, rmin};
    }

    // === Pair Parameter Addition ===

    /**
     * @brief Add bond parameters
     */
    void add_bond_params(const std::string& type1, const std::string& type2, 
                        double kb, double b0) {
        validate_pair_types(type1, type2, "bond parameters");
        auto key = ParamKeyUtils::make_pair(type1, type2);
        storage_.bond_params[key] = BondParams{kb, b0};
    }

    // === Triple Parameter Addition ===

    /**
     * @brief Add angle parameters
     */
    void add_angle_params(const std::string& type1, const std::string& type2, 
                         const std::string& type3, double ktheta, double theta0, 
                         double kub = 0.0, double s0 = 0.0) {
        validate_triple_types(type1, type2, type3, "angle parameters");
        auto key = ParamKeyUtils::make_triple(type1, type2, type3);
        storage_.angle_params[key] = AngleParams{ktheta, theta0, kub, s0};
    }

    // === Quadruple Parameter Addition ===

    /**
     * @brief Add dihedral parameters
     */
    void add_dihedral_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4,
                            double kchi, int n, double delta) {
        validate_quad_types(type1, type2, type3, type4, "dihedral parameters");
        auto key = ParamKeyUtils::make_quad(type1, type2, type3, type4);
        storage_.dihedral_params[key].emplace_back(kchi, n, delta);
    }

    /**
     * @brief Add improper parameters
     */
    void add_improper_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4,
                            double kpsi, double psi0) {
        validate_quad_types(type1, type2, type3, type4, "improper parameters");
        auto key = ParamKeyUtils::make_quad(type1, type2, type3, type4);
        storage_.improper_params[key] = ImproperParams{kpsi, psi0};
    }

    // === Batch Parameter Addition ===

    /**
     * @brief Add multiple atom masses at once
     */
    template<typename Container>
    void add_atom_masses_batch(const Container& mass_data) {
        for (const auto& item : mass_data) {
            add_atom_mass(item.first, item.second);
        }
    }

    /**
     * @brief Add multiple LJ parameters at once
     */
    template<typename Container>
    void add_lj_params_batch(const Container& lj_data) {
        for (const auto& item : lj_data) {
            const auto& type = item.first;
            const auto& params = item.second;
            add_lj_params(type, params.epsilon, params.rmin_half);
        }
    }

private:
    ForceFieldStorage& storage_;

    // === Validation Helpers ===

    void validate_non_empty_type(const std::string& type, const std::string& param_name) const {
        if (type.empty()) {
            throw std::invalid_argument("Atom type cannot be empty for " + param_name);
        }
    }

    void validate_pair_types(const std::string& type1, const std::string& type2, 
                           const std::string& param_name) const {
        if (type1.empty() || type2.empty()) {
            throw std::invalid_argument("Atom types cannot be empty for " + param_name);
        }
    }

    void validate_triple_types(const std::string& type1, const std::string& type2, 
                             const std::string& type3, const std::string& param_name) const {
        if (type1.empty() || type2.empty() || type3.empty()) {
            throw std::invalid_argument("Atom types cannot be empty for " + param_name);
        }
    }

    void validate_quad_types(const std::string& type1, const std::string& type2,
                           const std::string& type3, const std::string& type4,
                           const std::string& param_name) const {
        if (type1.empty() || type2.empty() || type3.empty() || type4.empty()) {
            throw std::invalid_argument("Atom types cannot be empty for " + param_name);
        }
    }
};

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_PARAMETER_ADDER_HPP 