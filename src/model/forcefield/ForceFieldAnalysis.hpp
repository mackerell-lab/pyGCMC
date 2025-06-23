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

class ForceField; // Forward declaration

/**
 * @brief Analysis utilities for ForceField
 */
class ForceFieldAnalysis {
public:
    static ForceFieldStats get_statistics(const ForceField& ff);
    static std::string get_summary(const ForceField& ff);
    static std::set<std::string> find_missing_lj_params(const ForceField& ff);
    static std::set<std::string> find_missing_masses(const ForceField& ff);
    static bool validate_force_field(const ForceField& ff);
    static CompletenessResult check_completeness(const ForceField& ff, const std::set<std::string>& atom_types);
    static std::vector<std::string> validate_consistency(const ForceField& ff);
};

} // namespace forcefield
} // namespace model
} // namespace pygcmc

// Implementation after ForceField is defined
#include "ForceFieldMain.hpp"

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
    oss << "Force Field Summary:\n  Atom types: " << stats.num_atom_types 
        << "\n  LJ parameters: " << stats.num_lj_params 
        << "\n  NBFIX entries: " << stats.num_nbfix
        << "\n  Bond types: " << stats.num_bond_types 
        << "\n  Angle types: " << stats.num_angle_types
        << "\n  Dihedral types: " << stats.num_dihedral_types 
        << "\n  Improper types: " << stats.num_improper_types;
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
    return find_missing_lj_params(ff).empty() && find_missing_masses(ff).empty();
}

inline CompletenessResult ForceFieldAnalysis::check_completeness(const ForceField& ff, const std::set<std::string>& atom_types) {
    CompletenessResult result;
    for (const auto& type : atom_types) {
        if (!ff.has_atom_mass(type)) result.missing_atom_masses.insert(type);
        if (!ff.has_lj_params(type)) result.missing_lj_params.insert(type);
    }
    result.is_complete = result.missing_atom_masses.empty() && result.missing_lj_params.empty();
    
    std::ostringstream oss;
    if (result.is_complete) {
        oss << "Force field is complete for all " << atom_types.size() << " atom types.";
    } else {
        oss << "Force field incomplete. Missing masses: " << result.missing_atom_masses.size()
            << ", missing LJ params: " << result.missing_lj_params.size();
    }
    result.summary = oss.str();
    return result;
}

inline std::vector<std::string> ForceFieldAnalysis::validate_consistency(const ForceField& ff) {
    std::vector<std::string> errors;
    for (const auto& type : find_missing_lj_params(ff)) {
        errors.push_back("Missing LJ parameters for atom type: " + type);
    }
    for (const auto& type : find_missing_masses(ff)) {
        errors.push_back("Missing mass for atom type: " + type);
    }
    return errors;
}

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_ANALYSIS_HPP 