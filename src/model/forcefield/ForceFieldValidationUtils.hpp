#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_VALIDATION_UTILS_HPP
#define PYGCMC_MODEL_FORCEFIELD_VALIDATION_UTILS_HPP

#include "ForceFieldValidationCore.hpp"
#include <limits>
#include <algorithm>

namespace pygcmc {
namespace model {
namespace forcefield {

// === Implementation of ForceFieldValidator utility methods ===

inline std::vector<std::string> ForceFieldValidator::validate_consistency() const {
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

inline std::vector<std::string> ForceFieldValidator::validate_nonbonded_params() const {
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

inline std::string ForceFieldValidator::get_validation_report() const {
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

/**
 * @brief Utility class for advanced ForceFieldValidator operations
 * 
 * This class provides additional utility functions for force field validation
 * that are separate from the core validation operations.
 */
class ForceFieldValidatorUtils {
public:
    /**
     * @brief Perform comprehensive validation with detailed reporting
     */
    static std::string get_comprehensive_validation_report(const ForceFieldValidator& validator) {
        std::ostringstream oss;
        oss << "=== Comprehensive Force Field Validation ===\n\n";

        // Basic statistics
        oss << "Storage Statistics:\n";
        oss << "  Atom types with masses: " << validator.storage_.atom_masses.size() << "\n";
        oss << "  Atom types with LJ params: " << validator.storage_.lj_params.size() << "\n";
        oss << "  NBFIX entries: " << validator.storage_.nbfix.size() << "\n";
        oss << "  Bond parameter types: " << validator.storage_.bond_params.size() << "\n";
        oss << "  Angle parameter types: " << validator.storage_.angle_params.size() << "\n";
        oss << "  Dihedral parameter types: " << validator.storage_.dihedral_params.size() << "\n";
        oss << "  Improper parameter types: " << validator.storage_.improper_params.size() << "\n\n";

        // Type consistency analysis
        auto missing_masses = validator.find_types_with_lj_but_no_mass();
        auto missing_lj = validator.find_types_with_mass_but_no_lj();

        if (!missing_masses.empty()) {
            oss << "Types with LJ params but no masses (" << missing_masses.size() << "):\n";
            for (const auto& type : missing_masses) {
                oss << "  - " << type << "\n";
            }
            oss << "\n";
        }

        if (!missing_lj.empty()) {
            oss << "Types with masses but no LJ params (" << missing_lj.size() << "):\n";
            for (const auto& type : missing_lj) {
                oss << "  - " << type << "\n";
            }
            oss << "\n";
        }

        // Parameter range analysis
        oss << analyze_parameter_ranges(validator);

        // Add basic validation report
        oss << validator.get_validation_report();

        return oss.str();
    }

    /**
     * @brief Validate specific molecule types
     */
    static CompletenessResult validate_molecule_types(const ForceFieldValidator& validator,
                                                     const std::vector<std::string>& molecule_types) {
        std::set<std::string> unique_types(molecule_types.begin(), molecule_types.end());
        return validator.check_completeness(unique_types);
    }

    /**
     * @brief Check force field compatibility with simulation requirements
     */
    static std::vector<std::string> check_simulation_compatibility(const ForceFieldValidator& validator,
                                                                  bool requires_bonded_params = false,
                                                                  bool requires_electrostatics = true) {
        std::vector<std::string> issues;

        // Basic requirements
        if (!validator.validate_force_field()) {
            issues.push_back("Force field lacks basic required parameters");
        }

        // Bonded parameter requirements
        if (requires_bonded_params) {
            if (validator.storage_.bond_params.empty()) {
                issues.push_back("Simulation requires bonded parameters but none found");
            }
        }

        // Electrostatics requirements
        if (requires_electrostatics) {
            const auto& nb = validator.storage_.nonbonded_params;
            if (nb.eps <= 0.0) {
                issues.push_back("Electrostatics enabled but dielectric constant is invalid");
            }
        }

        return issues;
    }

private:
    static std::string analyze_parameter_ranges(const ForceFieldValidator& validator) {
        std::ostringstream oss;
        oss << "Parameter Range Analysis:\n";

        // Analyze mass ranges
        if (!validator.storage_.atom_masses.empty()) {
            double min_mass = std::numeric_limits<double>::max();
            double max_mass = 0.0;
            for (const auto& [type, mass] : validator.storage_.atom_masses) {
                min_mass = std::min(min_mass, mass);
                max_mass = std::max(max_mass, mass);
            }
            oss << "  Mass range: " << min_mass << " - " << max_mass << " amu\n";
        }

        // Analyze LJ parameter ranges
        if (!validator.storage_.lj_params.empty()) {
            double min_eps = std::numeric_limits<double>::max();
            double max_eps = -std::numeric_limits<double>::max();
            double min_rmin = std::numeric_limits<double>::max();
            double max_rmin = 0.0;

            for (const auto& [type, params] : validator.storage_.lj_params) {
                min_eps = std::min(min_eps, params.epsilon);
                max_eps = std::max(max_eps, params.epsilon);
                min_rmin = std::min(min_rmin, params.rmin_half);
                max_rmin = std::max(max_rmin, params.rmin_half);
            }

            oss << "  LJ epsilon range: " << min_eps << " - " << max_eps << " kcal/mol\n";
            oss << "  LJ rmin/2 range: " << min_rmin << " - " << max_rmin << " Å\n";
        }

        oss << "\n";
        return oss.str();
    }
};

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_VALIDATION_UTILS_HPP 