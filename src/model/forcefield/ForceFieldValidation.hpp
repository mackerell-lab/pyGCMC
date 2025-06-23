#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_VALIDATION_HPP
#define PYGCMC_MODEL_FORCEFIELD_VALIDATION_HPP

#include "ForceFieldParams.hpp"
#include "ForceFieldInterface.hpp"
#include <sstream>

namespace pygcmc {
namespace model {
namespace forcefield {

/**
 * @brief Force field validation and completeness checking class
 * 
 * This class provides validation functionality for force fields including
 * completeness checks, consistency validation, and parameter verification.
 * All validation methods are read-only and safe for concurrent access.
 */
class ForceFieldValidator {
public:
    /**
     * @brief Constructor
     * @param storage Reference to the force field storage container
     */
    explicit ForceFieldValidator(const ForceFieldStorage& storage) : storage_(storage) {}
    ~ForceFieldValidator() = default;

    // === Basic Validation Methods ===

    /**
     * @brief Validate the force field has basic required parameters
     * @return true if force field has minimum required data
     */
    bool validate_force_field() const {
        // Basic validation: check if we have any atom masses and LJ parameters
        return !storage_.atom_masses.empty() && !storage_.lj_params.empty();
    }

    /**
     * @brief Check if force field is ready for simulations
     * @return true if force field has sufficient parameters for basic simulations
     */
    bool is_simulation_ready() const {
        // More stringent check: need masses and LJ params for the same types
        if (storage_.atom_masses.empty() || storage_.lj_params.empty()) {
            return false;
        }

        // Check if all atom types with masses have LJ parameters
        for (const auto& mass_pair : storage_.atom_masses) {
            if (storage_.lj_params.find(mass_pair.first) == storage_.lj_params.end()) {
                return false;
            }
        }

        return true;
    }

    // === Completeness Check Methods ===

    /**
     * @brief Check completeness for a given set of atom types
     * @param atom_types Set of atom types to check
     * @return CompletenessResult with detailed analysis
     */
    CompletenessResult check_completeness(const std::set<std::string>& atom_types) const {
        CompletenessResult result;
        result.is_complete = true;

        // Check atom masses
        for (const auto& type : atom_types) {
            if (storage_.atom_masses.find(type) == storage_.atom_masses.end()) {
                result.missing_atom_masses.insert(type);
                result.is_complete = false;
            }
        }

        // Check LJ parameters
        for (const auto& type : atom_types) {
            if (storage_.lj_params.find(type) == storage_.lj_params.end()) {
                result.missing_lj_params.insert(type);
                result.is_complete = false;
            }
        }

        // Generate summary
        std::ostringstream oss;
        if (result.is_complete) {
            oss << "Force field is complete for all " << atom_types.size() << " atom types.";
        } else {
            oss << "Force field is incomplete:\n";
            if (!result.missing_atom_masses.empty()) {
                oss << "  Missing atom masses: " << result.missing_atom_masses.size() << " types\n";
            }
            if (!result.missing_lj_params.empty()) {
                oss << "  Missing LJ parameters: " << result.missing_lj_params.size() << " types\n";
            }
        }
        result.summary = oss.str();

        return result;
    }

    /**
     * @brief Check completeness for nonbonded interactions only
     * @param atom_types Set of atom types to check
     * @return CompletenessResult focusing on nonbonded parameters
     */
    CompletenessResult check_nonbonded_completeness(const std::set<std::string>& atom_types) const {
        CompletenessResult result;
        result.is_complete = true;

        // Only check masses and LJ parameters for nonbonded interactions
        for (const auto& type : atom_types) {
            if (storage_.atom_masses.find(type) == storage_.atom_masses.end()) {
                result.missing_atom_masses.insert(type);
                result.is_complete = false;
            }
            if (storage_.lj_params.find(type) == storage_.lj_params.end()) {
                result.missing_lj_params.insert(type);
                result.is_complete = false;
            }
        }

        // Generate summary
        std::ostringstream oss;
        if (result.is_complete) {
            oss << "Nonbonded force field is complete for all " << atom_types.size() << " atom types.";
        } else {
            oss << "Nonbonded force field is incomplete:\n";
            if (!result.missing_atom_masses.empty()) {
                oss << "  Missing atom masses: " << result.missing_atom_masses.size() << " types\n";
            }
            if (!result.missing_lj_params.empty()) {
                oss << "  Missing LJ parameters: " << result.missing_lj_params.size() << " types\n";
            }
        }
        result.summary = oss.str();

        return result;
    }

    // === Consistency Validation Methods ===

    /**
     * @brief Check if force field parameters are internally consistent
     * @return Vector of issues found (empty if consistent)
     */
    std::vector<std::string> validate_consistency() const {
        std::vector<std::string> issues;

        // Check for missing masses
        auto missing_masses = find_types_with_lj_but_no_mass();
        if (!missing_masses.empty()) {
            std::ostringstream oss;
            oss << "Missing atom masses for " << missing_masses.size() << " types with LJ parameters";
            issues.push_back(oss.str());
        }

        // Check for missing LJ parameters
        auto missing_lj = find_types_with_mass_but_no_lj();
        if (!missing_lj.empty()) {
            std::ostringstream oss;
            oss << "Missing LJ parameters for " << missing_lj.size() << " types with masses";
            issues.push_back(oss.str());
        }

        // Check nonbonded parameter consistency
        const auto& nb = storage_.nonbonded_params;
        if (nb.cutnb <= nb.ctofnb) {
            issues.push_back("cutnb should be greater than ctofnb");
        }
        if (nb.ctofnb <= nb.ctonnb) {
            issues.push_back("ctofnb should be greater than ctonnb");
        }
        if (nb.cutnb <= 0.0) {
            issues.push_back("cutnb must be positive");
        }
        if (nb.eps <= 0.0) {
            issues.push_back("dielectric constant (eps) should be positive");
        }

        return issues;
    }

    /**
     * @brief Validate nonbonded parameter ranges
     * @return Vector of issues with nonbonded parameters
     */
    std::vector<std::string> validate_nonbonded_params() const {
        std::vector<std::string> issues;
        const auto& nb = storage_.nonbonded_params;

        // Check cutoff distances
        if (nb.cutnb < 5.0 || nb.cutnb > 50.0) {
            issues.push_back("cutnb outside typical range (5-50 Å)");
        }
        if (nb.ctofnb < 0.0) {
            issues.push_back("ctofnb cannot be negative");
        }
        if (nb.ctonnb < 0.0) {
            issues.push_back("ctonnb cannot be negative");
        }

        // Check dielectric and scaling factors
        if (nb.eps < 0.1 || nb.eps > 100.0) {
            issues.push_back("eps outside typical range (0.1-100)");
        }
        if (nb.e14fac < 0.0 || nb.e14fac > 2.0) {
            issues.push_back("e14fac outside typical range (0-2)");
        }

        return issues;
    }

    // === Helper Methods ===

    /**
     * @brief Find atom types that have LJ parameters but no masses
     * @return Set of problematic atom types
     */
    std::set<std::string> find_types_with_lj_but_no_mass() const {
        std::set<std::string> missing;
        for (const auto& pair : storage_.lj_params) {
            if (storage_.atom_masses.find(pair.first) == storage_.atom_masses.end()) {
                missing.insert(pair.first);
            }
        }
        return missing;
    }

    /**
     * @brief Find atom types that have masses but no LJ parameters
     * @return Set of problematic atom types
     */
    std::set<std::string> find_types_with_mass_but_no_lj() const {
        std::set<std::string> missing;
        for (const auto& pair : storage_.atom_masses) {
            if (storage_.lj_params.find(pair.first) == storage_.lj_params.end()) {
                missing.insert(pair.first);
            }
        }
        return missing;
    }

    /**
     * @brief Get summary of validation issues
     * @return Human-readable validation report
     */
    std::string get_validation_report() const {
        std::ostringstream oss;
        oss << "=== Force Field Validation Report ===\n\n";

        // Basic validation
        if (validate_force_field()) {
            oss << "✓ Basic validation: PASSED\n";
        } else {
            oss << "✗ Basic validation: FAILED\n";
        }

        // Simulation readiness
        if (is_simulation_ready()) {
            oss << "✓ Simulation readiness: PASSED\n";
        } else {
            oss << "✗ Simulation readiness: FAILED\n";
        }

        // Consistency check
        auto consistency_issues = validate_consistency();
        if (consistency_issues.empty()) {
            oss << "✓ Consistency check: PASSED\n";
        } else {
            oss << "✗ Consistency check: FAILED\n";
            for (const auto& issue : consistency_issues) {
                oss << "  - " << issue << "\n";
            }
        }

        // Nonbonded parameter validation
        auto nb_issues = validate_nonbonded_params();
        if (nb_issues.empty()) {
            oss << "✓ Nonbonded parameters: PASSED\n";
        } else {
            oss << "⚠ Nonbonded parameters: WARNINGS\n";
            for (const auto& issue : nb_issues) {
                oss << "  - " << issue << "\n";
            }
        }

        return oss.str();
    }

private:
    const ForceFieldStorage& storage_;  ///< Reference to parameter storage
};

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_VALIDATION_HPP 