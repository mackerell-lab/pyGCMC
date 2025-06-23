#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_ANALYSIS_HPP
#define PYGCMC_MODEL_FORCEFIELD_ANALYSIS_HPP

#include "ForceFieldTypes.hpp"
#include <set>
#include <string>
#include <vector>
#include <sstream>

namespace pygcmc {
namespace model {
namespace forcefield {

// Forward declaration
class ForceField;

/**
 * @brief Analysis and validation utilities for ForceField
 * 
 * This class provides analysis, validation, and statistical methods
 * for ForceField objects. It's separated from the main ForceField
 * class to keep the core functionality focused and maintainable.
 */
class ForceFieldAnalysis {
public:
    /**
     * @brief Get statistics about the force field parameters
     */
    static ForceFieldStats get_statistics(const ForceField& ff);

    /**
     * @brief Get a human-readable summary of the force field
     */
    static std::string get_summary(const ForceField& ff);

    /**
     * @brief Find atom types that have masses but no LJ parameters
     */
    static std::set<std::string> find_missing_lj_params(const ForceField& ff);

    /**
     * @brief Find atom types that have LJ parameters but no masses
     */
    static std::set<std::string> find_missing_masses(const ForceField& ff);

    /**
     * @brief Validate that the force field is consistent
     * @return true if all atom types have both masses and LJ parameters
     */
    static bool validate_force_field(const ForceField& ff);

    /**
     * @brief Check completeness for a specific set of atom types
     * @param ff The force field to check
     * @param atom_types Set of atom types to validate
     * @return CompletenessResult with detailed information
     */
    static CompletenessResult check_completeness(const ForceField& ff, 
                                                const std::set<std::string>& atom_types);

    /**
     * @brief Validate consistency and return detailed error messages
     * @return Vector of error messages (empty if no errors)
     */
    static std::vector<std::string> validate_consistency(const ForceField& ff);
};

} // namespace forcefield
} // namespace model
} // namespace pygcmc

// Include implementation after ForceField is defined
#include "ForceField.hpp"

namespace pygcmc {
namespace model {
namespace forcefield {

inline ForceFieldStats ForceFieldAnalysis::get_statistics(const ForceField& ff) {
    ForceFieldStats stats;
    stats.num_atom_types = ff.atom_masses_.size();
    stats.num_lj_params = ff.lj_params_.size();
    stats.num_nbfix = ff.nbfix_.size();
    stats.num_bond_types = ff.bond_params_.size();
    stats.num_angle_types = ff.angle_params_.size();
    stats.num_dihedral_types = ff.dihedral_params_.size();
    stats.num_improper_types = ff.improper_params_.size();
    return stats;
}

inline std::string ForceFieldAnalysis::get_summary(const ForceField& ff) {
    auto stats = get_statistics(ff);
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

inline std::set<std::string> ForceFieldAnalysis::find_missing_lj_params(const ForceField& ff) {
    std::set<std::string> missing;
    for (const auto& pair : ff.atom_masses_) {
        if (ff.lj_params_.find(pair.first) == ff.lj_params_.end()) {
            missing.insert(pair.first);
        }
    }
    return missing;
}

inline std::set<std::string> ForceFieldAnalysis::find_missing_masses(const ForceField& ff) {
    std::set<std::string> missing;
    for (const auto& pair : ff.lj_params_) {
        if (ff.atom_masses_.find(pair.first) == ff.atom_masses_.end()) {
            missing.insert(pair.first);
        }
    }
    return missing;
}

inline bool ForceFieldAnalysis::validate_force_field(const ForceField& ff) {
    auto missing_lj = find_missing_lj_params(ff);
    auto missing_masses = find_missing_masses(ff);
    return missing_lj.empty() && missing_masses.empty();
}

inline CompletenessResult ForceFieldAnalysis::check_completeness(const ForceField& ff, 
                                                                const std::set<std::string>& atom_types) {
    CompletenessResult result;
    
    for (const auto& type : atom_types) {
        if (!ff.has_atom_mass(type)) {
            result.missing_atom_masses.insert(type);
        }
        if (!ff.has_lj_params(type)) {
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

inline std::vector<std::string> ForceFieldAnalysis::validate_consistency(const ForceField& ff) {
    std::vector<std::string> errors;
    
    auto missing_lj = find_missing_lj_params(ff);
    auto missing_masses = find_missing_masses(ff);
    
    for (const auto& type : missing_lj) {
        errors.push_back("Missing LJ parameters for atom type: " + type);
    }
    
    for (const auto& type : missing_masses) {
        errors.push_back("Missing mass for atom type: " + type);
    }
    
    return errors;
}

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_ANALYSIS_HPP 